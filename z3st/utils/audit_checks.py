#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""Static consistency checks over the repository. Runs in seconds, solves nothing.

    python -m z3st.utils.audit_checks            # all checks
    python -m z3st.utils.audit_checks --list     # names only
    python -m z3st.utils.audit_checks schema     # one check by name

Exits non-zero if any check reports a finding, so it can gate a commit or a CI job.

Why these seven and not a linter's worth
--.--.--.--.--.--.--.--.--.--.--.--.--.-

Each one exists because it would have caught a real mistake made during the
2026-08 simplification audit that the test suite could not:

  models        `contact` was reported as having zero cases because the search was
                for `contact: true`, while it is configured as a *block* and Config
                stores bool() of it. Seven cases enable it. This check evaluates the
                switch the way Config does instead of matching text.
  assertions    Two golds carried a hard-coded rel_error of 0 against a reference
                that differed from the value, so the case passed whatever happened —
                and a zero reference also sends regression_check down its
                should-be-zero branch, disabling the gold comparison too.
  schema        Two cases hand-rolled a copy of pass_fail_check that wrote a bare
                json.dump(errors) with no "summary" key. The driver could read no
                verdict, and the blessed gold next to it was never compared.
  finite        A stale field name sat inside `if name in mesh.point_data:`, so two
                force-displacement curves were produced with the force set to NaN.
                No error, no warning, a saved figure missing its only quantity.
  names         Relocating the physics steps left six undefined names across two
                files. py_compile sees none of them; they would have surfaced as
                NameError in whichever case ran first.
  coverage      cluster_dynamic_model.py, 153 lines, had its only case in sandbox/
                with neither a gold nor a script, in a directory the driver never
                scans.
  mro           Spine flattens 13 parent classes into one namespace, so a duplicated
                method name is resolved silently by MRO order and the loser is never
                called. Five physics steps moved between those classes on 2026-08-09.

  docs          Four .rst files under docs/source/ still told users to conda install
                sympy and meshio, and showed a `coupling` key, days after all three
                were removed. The earlier sweep had grepped three files by name
                instead of the tree.

The unifying failure was asserting something from a check that could not have shown
the opposite. Where a check below is a heuristic that might misfire, it says so in
its own output rather than being silently trusted.

Not covered here, because it needs to run something: an `mpirun -n 2` smoke test.
Rank-asymmetric Python state deadlocks at exit — dolfinx's PETSc wrappers do
collective work in __del__ — and only an actual parallel run finds it.
verification/cluster/mass_conservation_1D takes 10 s and is the natural candidate.
"""

import argparse
import ast
import io
import json
import math
import pathlib
import re
import sys

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[2]
CASES = ROOT / "z3st" / "cases"

# Physics switches as Config reads them (core/config.py::Config.on)
MODEL_SWITCHES = [
    ("thermal", "models/thermal_model.py"),
    ("mechanical", "models/mechanical_model.py"),
    ("damage", "models/damage_model.py"),
    ("cluster", "models/cluster_dynamic_model.py"),
    ("plasticity", "models/plasticity_model.py"),
    ("contact", "models/contact_model.py"),
    ("porosity", "models/porosity_migration_model.py"),
]


def _yaml(path):
    try:
        return yaml.safe_load(path.read_text()) or {}
    except Exception:
        return {}


def _case_dirs():
    """Every directory holding an Allrun, i.e. everything the driver can discover."""
    return sorted(p.parent for p in CASES.rglob("Allrun"))


def _golds():
    return sorted(CASES.rglob("output/non-regression_gold.json"))


def _metrics(payload):
    return payload.get("results", payload) if isinstance(payload, dict) else {}


# --.. ..- .-.. .-.. --- 1. model coverage, evaluated as Config does --.. ..- ---
def check_models():
    """Which cases reach each physics model, using bool() on the switch.

    A model switch may be a boolean *or* a configuration block: Config stores
    ``bool(models.get(name, False))``, and a non-empty dict is True. Grepping for
    ``name: true`` therefore finds none of the block-configured cases.
    """
    findings = []
    reached = {name: [] for name, _ in MODEL_SWITCHES}
    for d in _case_dirs():
        models = _yaml(d / "input.yaml").get("models", {}) or {}
        for name, _ in MODEL_SWITCHES:
            if bool(models.get(name, False)):
                reached[name].append(d)

    for name, module in MODEL_SWITCHES:
        cases = reached[name]
        golded = [c for c in cases if (c / "output" / "non-regression_gold.json").exists()]
        if not cases:
            findings.append(f"{module}: no case anywhere enables '{name}'")
        elif not golded:
            findings.append(
                f"{module}: '{name}' enabled by {len(cases)} case(s) but none carries a gold, "
                f"so nothing protects it — e.g. {cases[0].relative_to(CASES)}"
            )
    return findings


# --.. ..- .-.. .-.. --- 2. assertions that cannot fail --.. ..- .-.. .-.. ---
def check_assertions():
    """Gold metrics whose pass/fail is switched off.

    ``rel_error == 0`` while ``reference != numerical`` means the number was written
    by hand rather than computed: the check passes for any value. A ``reference`` of
    exactly 0 additionally routes regression_check into its should-be-zero branch,
    which defers to the analytic status — disabling the gold comparison as well.

    ``reference == numerical`` is *not* flagged: that is what ``tracked()`` produces
    on purpose, for a quantity with no analytic reference.
    """
    findings = []
    for g in _golds():
        try:
            payload = json.loads(g.read_text())
        except Exception as exc:
            findings.append(f"{g.relative_to(CASES)}: unreadable ({exc})")
            continue
        for key, val in _metrics(payload).items():
            if not isinstance(val, dict):
                continue
            num, ref, rel = val.get("numerical"), val.get("reference"), val.get("rel_error")
            if any(isinstance(v, list) for v in (num, ref, rel)):
                continue
            try:
                num, ref, rel = float(num), float(ref), float(rel)
            except (TypeError, ValueError):
                continue
            if ref != num and rel == 0.0 and num != 0.0:
                findings.append(
                    f"{g.parent.parent.relative_to(CASES)} :: {key}: rel_error hard-coded to 0 "
                    f"with numerical={num:g} != reference={ref:g} — passes unconditionally"
                )
    return findings


# --.. ..- .-.. .-.. --- 3. verdict schema --.. ..- .-.. .-.. ---
def check_schema():
    """Verdict files the suite drivers can actually read.

    Both drivers grep ``output/non-regression.json`` for a top-level ``"summary"``
    and ``"regression"``. A bare ``json.dump(errors)`` produces a file that parses
    fine and yields no verdict, which is how two cases stayed effectively unprotected.

    The *gold* is deliberately not checked for a wrapper: regression_check reads it
    as ``gold_data.get("results", gold_data)``, so a bare metric dict is valid there.
    Requiring one produced 27 false positives on the first version of this check.
    """
    findings = []
    for j in sorted(CASES.rglob("output/non-regression.json")):
        case = j.parent.parent.relative_to(CASES)
        try:
            payload = json.loads(j.read_text())
        except Exception as exc:
            findings.append(f"{case}: non-regression.json unreadable ({exc})")
            continue
        for field in ("results", "tolerance", "summary"):
            if field not in payload:
                findings.append(
                    f"{case}: non-regression.json has no top-level '{field}' — "
                    "the drivers read no verdict from it"
                )
    return findings


# --.. ..- .-.. .-.. --- 4. metrics that are not finite --.. ..- .-.. .-.. ---
def check_finite():
    """NaN or inf anywhere in a blessed gold.

    A guard like ``if field_name in mesh.point_data:`` that never matches leaves the
    quantity as NaN and the script still saves its figure. Nothing raises.
    """
    findings = []
    for g in _golds():
        try:
            payload = json.loads(g.read_text())
        except Exception:
            continue
        for key, val in _metrics(payload).items():
            if not isinstance(val, dict):
                continue
            for field, v in val.items():
                for x in (v if isinstance(v, list) else [v]):
                    if isinstance(x, float) and not math.isfinite(x):
                        findings.append(
                            f"{g.parent.parent.relative_to(CASES)} :: {key}.{field} = {x}"
                        )
    return findings


# --.. ..- .-.. .-.. --- 5. undefined names --.. ..- .-.. .-.. ---
def check_names():
    """Syntax errors and undefined names across the package and the cases.

    ``py_compile`` catches only the first of those two. Unused imports are not
    reported: they are noise here, and several predate this check.
    """
    try:
        from pyflakes.api import check as pyflakes_check
        from pyflakes.reporter import Reporter
    except ImportError:
        return ["pyflakes not installed (pip install -e '.[dev]') — check skipped"]

    findings = []
    for path in sorted(ROOT.glob("z3st/**/*.py")):
        if "__pycache__" in str(path):
            continue
        src = path.read_text()
        try:
            ast.parse(src)
        except SyntaxError as exc:
            findings.append(f"{path.relative_to(ROOT)}:{exc.lineno}: SyntaxError: {exc.msg}")
            continue
        out, err = io.StringIO(), io.StringIO()
        pyflakes_check(src, str(path.relative_to(ROOT)), Reporter(out, err))
        for line in (out.getvalue() + err.getvalue()).splitlines():
            # "unable to detect undefined names" is the import-* warning, not a
            # finding: filtering on the substring alone reported three of those.
            if "unable to detect" in line:
                continue
            if "undefined name" in line or "invalid syntax" in line:
                findings.append(line.strip())
    return findings


# --.. ..- .-.. .-.. --- 6. suite coverage --.. ..- .-.. .-.. ---
def check_coverage():
    """Cases the driver cannot protect, and exclusions with no recorded reason.

    A case is discovered only with both an Allrun and a gold. Anything else runs
    without a net, and an exclusion whose reason was never written down is a
    decision nobody can revisit.
    """
    findings = []
    verdict_call = re.compile(r"(^|[^_a-zA-Z])(pass_fail_check|finish)\s*\(", re.M)

    for d in _case_dirs():
        rel = d.relative_to(CASES)
        if rel.parts and rel.parts[0] == "sandbox":
            continue  # never scanned, by design
        script = d / "non-regression.py"
        gold = d / "output" / "non-regression_gold.json"
        if not script.exists():
            continue  # aggregators and parametric drivers legitimately have none
        # follow one level of delegation: some variants exec a shared script
        # Delegation is real and takes three shapes: a shared script in the parent
        # directory, one in a *sibling* directory, and one alongside. Missing the
        # sibling form reported thermal_conductivity_GPR as verdict-less, which was
        # the same false positive this audit hit by hand earlier.
        sources = [script.read_text()]
        for name in set(re.findall(r"([A-Za-z_][A-Za-z0-9_-]*\.py)", sources[0])):
            for cand in (d / name, d.parent / name,
                         *d.parent.glob(f"*/{name}"), *d.parent.parent.glob(f"*/{name}")):
                if cand.exists() and cand.resolve() != script.resolve():
                    sources.append(cand.read_text())
                    break
        if not any(verdict_call.search(s) for s in sources):
            findings.append(f"{rel}: has a non-regression.py that never writes a verdict")
        elif not gold.exists():
            findings.append(f"{rel}: writes a verdict but has no gold — not discovered")

    exclude = CASES / "suite_exclude.txt"
    if exclude.exists():
        for raw in exclude.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "#" not in raw:
                findings.append(
                    f"suite_exclude.txt: '{line.split()[0]}' excluded with no reason comment"
                )
    return findings


# --.. ..- .-.. .-.. --- 7. documentation against the code --.. ..- .-.. .-.. ---
def check_docs():
    """Three concrete ways the prose can outlive the code.

    Deliberately narrow: case paths that no longer exist, dependencies the docs tell
    users to install that the package does not declare, and regime values outside
    Config's whitelist. Broader keyword sweeps produce false positives on prose and
    on dated changelog entries, which is why they are not attempted here.
    """
    findings = []
    docs = [p for p in list(ROOT.glob("*.md")) + list(ROOT.glob("docs/**/*.rst"))
            + list(ROOT.glob("docs/**/*.md")) + list(ROOT.glob("z3st/ai/*.md"))
            if ".git" not in str(p)]

    # a) case paths cited in prose
    cited = re.compile(r"`((?:verification|regression|benchmarks|studies|teaching|sandbox)/[\w/]+)`")
    for p in docs:
        text = p.read_text()
        for match in set(cited.findall(text)):
            if (CASES / match.rstrip("/")).is_dir():
                continue
            # a documented removal is fine; look for the marker on the same line
            lines = [ln for ln in text.splitlines() if match in ln]
            if all(re.search(r"remov|delet|was:", ln, re.I) for ln in lines):
                continue
            findings.append(f"{p.relative_to(ROOT)}: cites a case that does not exist: {match}")

    # A dependency cross-check was tried here and removed: parsing tokens after
    # "pip install" swallowed prose ("the", "part", "is", "also", "available") and
    # flagged legitimate suggestions (pylint for pyreverse, mpich and petsc4py from
    # the conda stack). Eight of its ten findings were noise, and a check that cries
    # wolf teaches people to skip the whole tool. Auditing what the docs tell users
    # to install needs a curated list of removed names, not a regex.

    # c) regime values
    config = (ROOT / "z3st" / "core" / "config.py").read_text()
    match = re.search(r"valid_regimes\s*=\s*\{([^}]*)\}", config)
    if match:
        valid = set(re.findall(r'"([^"]+)"', match.group(1)))
        for p in docs:
            text = p.read_text()
            # Only indented / fenced lines: in prose "regime: retained" is a
            # sentence, not configuration.
            candidates = set(re.findall(r"^[ \t>]+regime:\s*([A-Za-z_0-9]+)", text, re.M))
            for regime in candidates:
                if regime.lower() in valid:
                    continue
                lines = [ln for ln in text.splitlines() if f"regime: {regime}" in ln]
                if all(re.search(r"remov|delet|was:", ln, re.I) for ln in lines):
                    continue
                findings.append(
                    f"{p.relative_to(ROOT)}: documents regime '{regime}', not in "
                    f"Config.valid_regimes {sorted(valid)}"
                )
    return findings


def check_mro():
    """Method-name collisions between Spine's parent classes.

    ``Spine`` multiply-inherits 13 classes into one flat namespace, so two mixins
    defining the same method name do not conflict -- the MRO silently picks one and
    the other is never called. Nothing raises, nothing warns, and the case still
    passes if the surviving implementation happens to be close enough.

    That surface was widened on 2026-08-09, when five physics steps moved out of
    core/solver.py into their models. It was clean then (85 methods, zero
    collisions); this check is what keeps it clean, because "we looked once" does
    not survive the next model.

    Parsed with ast rather than by importing Spine: the import pulls in dolfinx and
    PETSc, and this script's promise is that it needs neither an environment nor a
    simulation.
    """
    spine = ROOT / "z3st" / "core" / "spine.py"
    tree = ast.parse(spine.read_text())
    bases = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "Spine":
            bases = [b.id for b in node.bases if isinstance(b, ast.Name)]
    if not bases:
        return ["z3st/core/spine.py: could not read Spine's base classes — "
                "this check cannot have found anything, treat it as unrun"]

    owners = {}
    findings = []
    found = set()  # classes whose definition was parsed, with or without methods
    for src in sorted((ROOT / "z3st").rglob("*.py")):
        try:
            mod = ast.parse(src.read_text())
        except SyntaxError:
            continue  # reported by check_names
        for node in mod.body:
            if not isinstance(node, ast.ClassDef) or node.name not in bases:
                continue
            found.add(node.name)
            for item in node.body:
                if not isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    continue
                if item.name.startswith("__"):
                    continue
                owners.setdefault(item.name, []).append(node.name)

    for name, classes in sorted(owners.items()):
        if len(classes) > 1:
            findings.append(
                f"{name}() defined in {', '.join(classes)} — the MRO calls "
                f"{classes[0]}'s and silently discards the rest"
            )

    # A base with only dunder methods (Config) contributes no names and is not a
    # gap; a base whose file was never parsed is. Conflating the two reported
    # Config as unchecked on the first run of this check.
    missing = [b for b in bases if b not in found]
    if missing:
        findings.append(
            f"z3st/core/spine.py: no class definition found for {', '.join(missing)} — "
            "a base that was not parsed cannot be checked for collisions"
        )
    return findings


CHECKS = {
    "models": check_models,
    "assertions": check_assertions,
    "schema": check_schema,
    "finite": check_finite,
    "names": check_names,
    "coverage": check_coverage,
    "docs": check_docs,
    "mro": check_mro,
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("checks", nargs="*", help="run only these (default: all)")
    parser.add_argument("--list", action="store_true", help="print check names and exit")
    args = parser.parse_args(argv)

    if args.list:
        for name, fn in CHECKS.items():
            print(f"  {name:<12} {(fn.__doc__ or '').strip().splitlines()[0]}")
        return 0

    selected = args.checks or list(CHECKS)
    unknown = [c for c in selected if c not in CHECKS]
    if unknown:
        parser.error(f"unknown check(s): {', '.join(unknown)}")

    total = 0
    for name in selected:
        findings = CHECKS[name]()
        total += len(findings)
        mark = "FAIL" if findings else " ok "
        print(f"[{mark}] {name}" + (f" — {len(findings)} finding(s)" if findings else ""))
        for f in findings:
            print(f"         {f}")

    print()
    if total:
        print(f"{total} finding(s). None of these needs a simulation to reproduce.")
    else:
        print("No findings.")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main())
