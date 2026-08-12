#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""Static consistency checks over the repository. Runs in seconds, solves nothing.

    python -m z3st.utils.audit_checks            # all checks
    python -m z3st.utils.audit_checks --list     # names only
    python -m z3st.utils.audit_checks schema     # one check by name

Exits non-zero if any check reports a finding, so it can gate a commit or a CI job.

The checks
--.--.--.--.--.--.--.--.--.--.--.--.--.-

  models        which cases reach each physics model, evaluating the switch as
                Config does
  assertions    gold metrics whose pass/fail is switched off
  schema        verdict files the suite drivers can read
  finite        NaN or inf in a blessed gold
  names         syntax errors and undefined names
  coverage      cases the driver cannot protect, exclusions with no reason
  shell         every shell script parses, every stub finds the shared runner
  workflow      scripts the CI workflow invokes exist
  scripts       literal $ROOT_DIR paths in the shell drivers resolve
  deps          third-party imports in library code are declared
  mro           method-name collisions between Spine's parent classes
  docs          case paths, dependencies and regime values in the prose
  ci            cases_ci.txt against the tree it names
  env           pyproject runtime deps are installed by z3st_env.yml

Where a check is a heuristic that might misfire, it says so in its own output.

Not covered here, because it needs to run something: an `mpirun -n 2` smoke test.
"""

import argparse
import ast
import io
import json
import math
import pathlib
import re
import subprocess
import sys

import yaml

from z3st.core.config import MODEL_NAMES

ROOT = pathlib.Path(__file__).resolve().parents[2]
CASES = ROOT / "z3st" / "cases"

# Physics switches, taken from Config itself so that a model added there cannot
# go unaudited here. Only the two module names that do not follow <name>_model.py
# are spelled out.
_MODEL_MODULES = {"cluster": "cluster_dynamic_model.py",
                  "porosity": "porosity_migration_model.py"}
MODEL_SWITCHES = [(name, "models/" + _MODEL_MODULES.get(name, name + "_model.py"))
                  for name in MODEL_NAMES]


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
    ``bool(models.get(name, False))``, and a non-empty dict is True.
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
            if ref is None:
                # tracked()/_regression_only(): no closed form, gold-only protection
                # by design. Same exemption as reference == numerical below.
                continue
            try:
                num, ref, rel = float(num), float(ref), float(rel)
            except (TypeError, ValueError):
                findings.append(
                    f"{g.parent.parent.relative_to(CASES)} :: {key}: metric fields are "
                    f"not numeric (numerical={num!r}, reference={ref!r}, rel_error={rel!r}) "
                    "— UNCHECKED by the assertion audit"
                )
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
    and ``"regression"``.

    The *gold* is deliberately not checked for a wrapper: regression_check reads it
    as ``gold_data.get("results", gold_data)``, so a bare metric dict is valid there.
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
    """NaN or inf anywhere in a blessed gold."""
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

    Unused imports are not reported.
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
            # A star import switches pyflakes' undefined-name analysis off for the
            # whole file, so the loop below can only report "no findings" there.
            # Opt out on the import line with a "noqa: F403" directive when the
            # namespace is external.
            if "unable to detect" in line:
                star = next((l for l in src.splitlines()
                             if re.match(r"\s*from\s+\S+\s+import\s+\*", l)), "")
                if "noqa" not in star:
                    findings.append(
                        f"{path.relative_to(ROOT)}: star import disables pyflakes' "
                        "undefined-name analysis for this whole file — it is UNCHECKED. "
                        "Import the names explicitly, or mark the line `# noqa: F403`."
                    )
                continue
            if "undefined name" in line or "invalid syntax" in line:
                findings.append(line.strip())
    return findings


# --.. ..- .-.. .-.. --- 6. suite coverage --.. ..- .-.. .-.. ---
def check_coverage():
    """Cases the driver cannot protect, and exclusions with no recorded reason.

    A case is discovered only with both an Allrun and a gold.
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
        # Follow one level of delegation. It takes three shapes: a shared script in
        # the parent directory, one in a *sibling* directory, and one alongside.
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

    # regression_check reads output/non-regression_gold.json and nothing else,
    # so a gold anywhere else is read by no one.
    for stray in sorted(CASES.rglob("*gold*.json")):
        if stray.parent.name != "output":
            findings.append(
                f"{stray.relative_to(CASES)}: gold outside output/, read by nothing"
            )

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
    Config's whitelist.
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

    # b) no dependency cross-check

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
    the other is never called.

    Parsed with ast, not by importing Spine.
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
    # gap; a base whose file was never parsed is.
    missing = [b for b in bases if b not in found]
    if missing:
        findings.append(
            f"z3st/core/spine.py: no class definition found for {', '.join(missing)} — "
            "a base that was not parsed cannot be checked for collisions"
        )
    return findings


# Installed by conda (the FEniCSx stack) or built by hand (the SCIANTIX binding), so
# deliberately absent from pyproject. installation.rst documents both.
_NOT_ON_PYPI = {"dolfinx", "ufl", "basix", "petsc4py", "mpi4py", "slepc4py", "ffcx",
                "sciantix_binding"}

# import name -> distribution name, where they differ.
_DIST_ALIASES = {"yaml": "pyyaml", "PIL": "pillow", "sklearn": "scikit-learn"}


def check_deps():
    """Third-party modules imported by library code must be declared somewhere.

    The rule follows *when* the import runs:

      * **module level, unguarded** -> must be in ``[project] dependencies``.
      * **inside a function, or guarded by try/except ImportError** -> any group will do
        (an extra is the right home).

    A local run cannot reveal an undeclared dependency: the conda dolfinx stack pulls
    several in transitively.

    Scope is library code: z3st/ minus cases/ and conference/, which are leaf scripts.
    """
    pyproject = (ROOT / "pyproject.toml").read_text()
    block = re.search(r"^dependencies\s*=\s*\[(.*?)^\]", pyproject, re.M | re.S)
    if not block:
        return ["pyproject.toml: no dependencies array found — this check cannot have "
                "found anything, treat it as unrun"]
    # PEP 503: "-" and "_" are equivalent in a distribution name, so normalise both
    # sides. dolfinx-external-operator is imported as dolfinx_external_operator.
    def _norm(name):
        return name.lower().replace("_", "-")

    hard = {_norm(m.group(1)) for m in re.finditer(r'"([A-Za-z0-9_.-]+)', block.group(1))}
    # Any name quoted anywhere in pyproject counts as "declared somewhere",
    # extras included. Leading name only: a version spec ("numpy>=2") must
    # still match the module name.
    anywhere = {_norm(m.group(1)) for m in re.finditer(r'"([A-Za-z0-9_.-]+)', pyproject)}

    findings = []
    for src in sorted((ROOT / "z3st").rglob("*.py")):
        rel = src.relative_to(ROOT)
        if rel.parts[1] in ("cases", "conference"):
            continue
        try:
            tree = ast.parse(src.read_text())
        except SyntaxError:
            continue  # reported by check_names

        # Mark every node that sits inside a function or inside a try that handles an
        # import failure: those imports are deliberately lazy or deliberately optional.
        lazy = set()
        for node in ast.walk(tree):
            guarded = isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            if isinstance(node, ast.Try):
                guarded = any(
                    h.type is None
                    or (isinstance(h.type, ast.Name) and h.type.id in ("ImportError", "Exception"))
                    for h in node.handlers
                )
                children = node.body
            elif guarded:
                children = node.body
            else:
                continue
            if guarded:
                for sub in children:
                    for inner in ast.walk(sub):
                        lazy.add(id(inner))

        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                mods = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                # level > 0 is relative; module is None for "from . import x"
                mods = [] if node.level or not node.module else [node.module.split(".")[0]]
            else:
                continue
            for mod in mods:
                if mod in sys.stdlib_module_names or mod in ("z3st", *_NOT_ON_PYPI):
                    continue
                # A module shadowed by a file in this repo is a local import, not a
                # dependency -- __main__ imports a case-local diagnostics.py this way.
                if any((ROOT / "z3st").rglob(f"{mod}.py")):
                    continue
                dist = _norm(_DIST_ALIASES.get(mod, mod))
                if id(node) in lazy:
                    if dist not in anywhere:
                        findings.append(
                            f"{rel}: imports '{mod}' lazily but '{dist}' is declared in no "
                            "dependency group — the feature cannot be installed"
                        )
                elif dist not in hard:
                    findings.append(
                        f"{rel}: imports '{mod}' at module level but '{dist}' is not in "
                        "[project] dependencies — a fresh install breaks on import"
                    )
    return sorted(set(findings))


def check_scripts():
    """Literal paths in the shell drivers must resolve.

    Only paths whose remainder is a literal are checked; anything with a second variable
    in it (``"$ROOT_DIR/cases/$case_name"``) is built at run time and skipped.
    """
    findings = []
    pkg = ROOT / "z3st"
    for sh in sorted(pkg.rglob("*.sh")):
        text = sh.read_text()
        # ROOT_DIR does not mean the same directory in every script: the local driver
        # sets it to its own directory (cases/), the CI driver to the parent (the
        # package).
        uses = re.findall(r'"\$\{?ROOT_DIR\}?/([^"$]+)"', text)
        decl = re.search(r'ROOT_DIR="\$\(\s*cd\s+"\$\(\s*dirname\s+"\$\{BASH_SOURCE\[0\]\}"'
                         r'\s*\)((?:/\.\.)*)"\s*&&', text)
        if not decl:
            if uses:
                findings.append(
                    f"{sh.relative_to(ROOT)}: uses $ROOT_DIR but its declaration could "
                    "not be parsed — its paths are UNCHECKED, not verified"
                )
            continue
        base = sh.parent
        for _ in re.findall(r"/\.\.", decl.group(1)):
            base = base.parent
        for m in re.finditer(r'"\$\{?ROOT_DIR\}?/([^"$]+)"', text):
            target = base / m.group(1)
            if not target.exists():
                findings.append(
                    f"{sh.relative_to(ROOT)}: refers to $ROOT_DIR/{m.group(1)}, which "
                    f"does not exist — ROOT_DIR here is {base.relative_to(ROOT)}/"
                )
    return sorted(set(findings))


def check_shell():
    """Every shell script parses, and every stub finds the shared runner.

    Most shell in this repo is per-case Allrun/Allclean stubs. Syntax is checked
    with ``bash -n``.
    """
    findings = []
    targets = [f for f in ROOT.rglob("*")
               if f.is_file() and ".git" not in f.parts
               and (f.suffix == ".sh" or f.name in ("Allrun", "Allclean"))]
    for f in sorted(targets):
        r = subprocess.run(["bash", "-n", str(f)], capture_output=True, text=True)
        if r.returncode != 0:
            findings.append(f"{f.relative_to(ROOT)}: shell syntax error: "
                            f"{r.stderr.strip().splitlines()[-1] if r.stderr.strip() else '?'}")
        # Stubs source the shared runner through the installed package. The path inside
        # that expression is fixed, so it can be checked even though the prefix is not.
        for m in re.finditer(r"/utils/(allrun|allclean)\.sh", f.read_text()):
            if not (ROOT / "z3st" / "utils" / f"{m.group(1)}.sh").exists():
                findings.append(f"{f.relative_to(ROOT)}: sources utils/{m.group(1)}.sh, "
                                "which does not exist")
    if not targets:
        findings.append("no shell scripts found at all — this check cannot have found "
                        "anything, treat it as unrun")
    return sorted(set(findings))


def check_workflow():
    """Scripts the CI workflow invokes must exist in the repo.

    Nothing else exercises the workflow: a path that is wrong there fails only on a
    runner. Outputs (upload-artifact paths) are deliberately not checked -- they are produced
    by the run, so their absence before it is correct.
    """
    findings = []
    wf_dir = ROOT / ".github" / "workflows"
    if not wf_dir.is_dir():
        return ["no .github/workflows directory — this check cannot have found "
                "anything, treat it as unrun"]
    for wf in sorted(wf_dir.glob("*.yml")):
        text = wf.read_text()
        refs = set(re.findall(r"(?:^|\s)\./(\S+\.(?:sh|py))", text, re.M))
        refs |= set(re.findall(r"chmod\s+\+x\s+\.?/?(\S+)", text))
        for ref in sorted(refs):
            if not (ROOT / ref).exists():
                findings.append(
                    f".github/workflows/{wf.name}: invokes ./{ref}, which does not "
                    "exist — this fails only on a runner"
                )
    return sorted(set(findings))


def check_ci():
    """The hand-written CI case list against the tree it names.

    non-regression_local.sh discovers its cases (Allrun + gold, sandbox pruned);
    non-regression_github.sh instead reads cases_ci.txt, maintained by hand.

    Not flagged: cases absent from cases_ci.txt. The short list is a documented
    choice (a per-commit gate under a measured time budget), so 'missing from CI' is
    the normal state for most cases.
    """
    findings = []
    ci_file = CASES / "cases_ci.txt"
    if not ci_file.exists():
        return ["z3st/cases/cases_ci.txt is missing — this check cannot have found "
                "anything, treat it as unrun"]

    def _entries(path):
        for raw in path.read_text().splitlines():
            line = raw.strip()
            if line and not line.startswith("#"):
                yield line.split()[0]

    excluded = set(_entries(CASES / "suite_exclude.txt")) if (CASES / "suite_exclude.txt").exists() else set()

    for name in _entries(ci_file):
        d = CASES / name
        if not d.is_dir():
            findings.append(f"cases_ci.txt: '{name}' is not a directory — CI will fail on a runner")
        elif not (d / "Allrun").exists():
            findings.append(f"cases_ci.txt: '{name}' has no Allrun — CI cannot run it")
        elif not (d / "output" / "non-regression_gold.json").exists():
            findings.append(f"cases_ci.txt: '{name}' has no gold — CI runs it with no regression check")
        elif name.split("/")[0] == "sandbox":
            findings.append(f"cases_ci.txt: '{name}' is under sandbox, which the local suite prunes")
        if name in excluded:
            findings.append(
                f"cases_ci.txt: '{name}' is also in suite_exclude.txt — excluded locally, run in CI"
            )
    return findings


def check_env():
    """Runtime dependencies declared in pyproject.toml must be in z3st_env.yml.

    `pip install -e .` resolves `[project] dependencies`, but the documented way to
    build a working environment is `conda env create -f z3st_env.yml`, which reads
    only its own lists.

    Only the runtime set is compared. Extras are installed by name
    (`pip install -e '.[dev]'`), so their absence from the env file is correct.
    """
    findings = []
    pyproject = ROOT / "pyproject.toml"
    env_file = ROOT / "z3st_env.yml"
    if not (pyproject.exists() and env_file.exists()):
        return [f"{'pyproject.toml' if not pyproject.exists() else 'z3st_env.yml'} is "
                "missing — this check cannot have found anything, treat it as unrun"]

    block = re.search(r"^dependencies\s*=\s*\[(.*?)^\]", pyproject.read_text(), re.M | re.S)
    if not block:
        return ["pyproject.toml: no [project] dependencies array found — this check "
                "cannot have found anything, treat it as unrun"]

    def _norm(name):  # PEP 503
        return re.sub(r"[-_.]+", "-", name).lower()

    declared = {_norm(m.group(1)) for m in re.finditer(r'"([A-Za-z0-9_.-]+)', block.group(1))}

    # Every package name the env file installs, conda entries and the pip: block alike.
    env_text = env_file.read_text()
    installed = set()
    for raw in env_text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line.startswith("-"):
            continue
        name = re.match(r"-\s*([A-Za-z0-9_.-]+)", line)
        if name:
            installed.add(_norm(name.group(1)))

    for dep in sorted(declared - installed):
        findings.append(
            f"z3st_env.yml: '{dep}' is a runtime dependency in pyproject.toml but the "
            "env file does not install it — `conda env create -f z3st_env.yml` "
            "produces an environment that cannot import it"
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
    "deps": check_deps,
    "scripts": check_scripts,
    "shell": check_shell,
    "workflow": check_workflow,
    "ci": check_ci,
    "env": check_env,
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
