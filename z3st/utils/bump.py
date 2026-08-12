#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""Propagate the version across the repository, or check that it already agrees.

    python -m z3st.utils.bump --check     # every tracked file agrees with pyproject
    python -m z3st.utils.bump 0.3.3       # rewrite it everywhere
    python -m z3st.utils.bump 0.3.3 -n    # show what would change, touch nothing

The version lives in pyproject.toml and is repeated in ~130 tracked files, mostly
the header banner. --check is the guard: it fails if a tag is cut while the tree
still names the previous release.

Only tracked text files are visited, so build directories and outputs stay out.
CITATION.cff also carries date-released, which is set to today on a rewrite.
"""

import argparse
import datetime
import pathlib
import re
import subprocess
import sys
import tomllib

ROOT = pathlib.Path(__file__).resolve().parents[2]
SEMVER = re.compile(r"^\d+\.\d+\.\d+$")


def current_version():
    """The version in pyproject.toml, which is the one every other file follows."""
    with open(ROOT / "pyproject.toml", "rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def tracked_files():
    """Tracked paths, from git, so untracked output and build trees stay out."""
    out = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "-z"],
        capture_output=True, text=True, check=True,
    ).stdout
    return [ROOT / name for name in out.split("\0") if name]


def occurrences(text, version):
    """Spans of `version` not embedded in a longer number: 0.3.2 must not match
    inside 0.3.20, and 10.3.2 must not match at its tail."""
    pattern = re.compile(rf"(?<![\d.]){re.escape(version)}(?![\d.])")
    return list(pattern.finditer(text))


def visit(version, new=None):
    """Find `version` in every tracked text file, rewriting it to `new` if given.
    Returns {path: count}. Binary and unreadable files are skipped."""
    hits = {}
    for path in tracked_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        found = occurrences(text, version)
        if not found:
            continue
        hits[path] = len(found)
        if new is not None:
            path.write_text(
                re.sub(rf"(?<![\d.]){re.escape(version)}(?![\d.])", new, text),
                encoding="utf-8",
            )
    return hits


# Where a version is declared rather than merely mentioned. Prose is free to
# cite an older release on purpose; these forms are not.
# The lookbehind keeps a qualified neighbour out: cff-version and mermaid_version
# declare someone else's version, not ours.
DECLARATIONS = re.compile(
    r"""^.*?(?<![\w-])(?:
          Version:\s*                    # header banner
        | version\s*=\s*"                # pyproject
        | version:\s*                    # CITATION.cff
        | release\s*=\s*"                # Sphinx conf.py
        | __version__\s*=\s*"            # package fallback
    )(\d+\.\d+\.\d+)""",
    re.M | re.X,
)


def disagreements(version):
    """Version declarations that name something other than `version`.
    Returns [(path, line number, what it says)]."""
    bad = []
    for path in tracked_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        for match in DECLARATIONS.finditer(text):
            if match.group(1) != version:
                line = text.count("\n", 0, match.start()) + 1
                bad.append((path, line, match.group(1)))
    return bad


def stamp_citation(today):
    """CITATION.cff carries a release date alongside the version."""
    path = ROOT / "CITATION.cff"
    text = path.read_text(encoding="utf-8")
    text, n = re.subn(r"^date-released: .*$", f"date-released: {today}",
                      text, flags=re.M)
    path.write_text(text, encoding="utf-8")
    return n


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("version", nargs="?", help="the new version, e.g. 0.3.3")
    parser.add_argument("--check", action="store_true",
                        help="report disagreement with pyproject, change nothing")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="list the files a rewrite would touch")
    args = parser.parse_args(argv)

    old = current_version()

    if args.check:
        stale = disagreements(old)
        print(f"pyproject: {old}")
        for path, line, found in stale:
            print(f"  [FAIL] {path.relative_to(ROOT)}:{line} names {found}")
        if stale:
            print(f"{len(stale)} declaration(s) disagree — run: "
                  f"python -m z3st.utils.bump {old}", file=sys.stderr)
            return 1
        print("every version declaration agrees")
        return 0

    if not args.version:
        parser.error("give a version to set, or --check")
    if not SEMVER.match(args.version):
        parser.error(f"'{args.version}' is not X.Y.Z")
    if args.version == old:
        print(f"already at {old}, nothing to do")
        return 0

    if args.dry_run:
        hits = visit(old)
        for path, count in sorted(hits.items()):
            print(f"  {count:>2}  {path.relative_to(ROOT)}")
        print(f"{old} -> {args.version}: "
              f"{sum(hits.values())} occurrences in {len(hits)} files")
        return 0

    hits = visit(old, args.version)
    today = datetime.date.today().isoformat()
    stamp_citation(today)
    print(f"{old} -> {args.version}: "
          f"{sum(hits.values())} occurrences in {len(hits)} files")
    print(f"CITATION.cff date-released: {today}")

    left = visit(old)
    if left:
        print(f"[FAIL] {len(left)} files still name {old}", file=sys.stderr)
        return 1
    return 0


def demo():
    """Self-check on the substring edge cases, which is where a naive replace
    silently corrupts a file. Run: python -m z3st.utils.bump --demo"""
    assert len(occurrences("Version: 0.3.2 (2026)", "0.3.2")) == 1
    assert len(occurrences("0.3.20 is not it", "0.3.2")) == 0
    assert len(occurrences("see 10.3.2 above", "0.3.2")) == 0
    assert len(occurrences("10.5281/zenodo.17748028", "0.3.2")) == 0
    assert len(occurrences("v0.3.2 and 0.3.2", "0.3.2")) == 2
    assert len(occurrences("", "0.3.2")) == 0
    print("bump: self-check ok")


if __name__ == "__main__":
    if "--demo" in sys.argv:
        demo()
    else:
        sys.exit(main())
