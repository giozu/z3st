# SPDX-License-Identifier: Apache-2.0
"""
Isolated parametric sweeps for z3st cases.

Parametric study scripts (``study_theta.py``, ``study_jiang.py``,
``study_pressure_sweep.py``, ...) historically swept a parameter by rewriting
the case's own tracked ``mesh.geo`` / ``geometry.yaml`` / ``input.yaml`` in
place, once per sample, and left the case in whatever state the final sample
wrote. That has two failure modes:

* it destroys parameterisation -- a writer that only knows how to emit
  ``ay = <number>`` will happily overwrite the line that *computed* ``ay`` from
  a swept angle, leaving a case that still meshes but silently ignores the
  parameter being swept;
* it corrupts silently on interruption -- a sweep that dies at sample 2 of 3
  leaves the case configured for sample 2, unmarked, so the next ``./Allrun``
  produces a plausible result for the wrong geometry.

This module runs each sample in a throwaway copy instead. The tracked case
files are read, never written.

Variants are created as ``<parent>/.sweep/<case>__<label>/`` -- one directory
level deeper than the case itself -- so relative references inside the copied
files (``../../../../materials/uo2_jiang.yaml``) are re-targeted by exactly one
level. ``.sweep/`` is gitignored.

Typical use::

    from z3st.utils.case_sweep import make_variant, run_variant

    for theta in (30, 45, 60, 90):
        v = make_variant(CASE_DIR, f"theta_{theta}",
                         edits={"mesh.geo": lambda s: set_theta(s, theta)})
        run_variant(v)
        harvest(os.path.join(v, "output"))
"""

import os
import re
import shutil
import subprocess

#: Files copied into a variant. Anything not listed (logs, meshes, previous
#: results) is deliberately left behind so each variant starts clean.
CASE_FILES = (
    "input.yaml",
    "geometry.yaml",
    "mesh.geo",
    "boundary_conditions.yaml",
    "non-regression.py",
    "make_bcs.py",
    "diagnostics.py",
    "generate_yaml.py",
    "Allclean",
    "Allrun",
)

#: How many directory levels a variant sits below its case dir.
_EXTRA_DEPTH = 1

# Relative path references inside case files, e.g.
#   materials: ../../../../materials/uo2_jiang.yaml
#   python3 ../../../../utils/plot_convergence.py
_REL_PATH = re.compile(r'(?<![\w./])((?:\.\./)+)')


def _retarget(text):
    """Deepen every relative ``../`` chain by :data:`_EXTRA_DEPTH` levels."""
    return _REL_PATH.sub(lambda m: "../" * _EXTRA_DEPTH + m.group(1), text)


def sweep_root(case_dir):
    """Directory holding all variants of *case_dir*."""
    return os.path.join(os.path.dirname(os.path.abspath(case_dir)), ".sweep")


def make_variant(case_dir, label, edits=None):
    """Materialise an isolated copy of *case_dir* and return its path.

    Parameters
    ----------
    case_dir : str
        The case to copy. Never modified.
    label : str
        Identifies this sample; becomes part of the variant directory name.
    edits : dict, optional
        ``{filename: callable}``. Each callable takes the file's text and
        returns the modified text. Applied to the *copy*, after re-targeting.
    """
    case_dir = os.path.abspath(case_dir)
    name = os.path.basename(case_dir)
    variant = os.path.join(sweep_root(case_dir), f"{name}__{label}")

    if os.path.isdir(variant):
        shutil.rmtree(variant)
    os.makedirs(os.path.join(variant, "output"))

    for fname in CASE_FILES:
        src = os.path.join(case_dir, fname)
        if not os.path.isfile(src):
            continue
        dst = os.path.join(variant, fname)
        with open(src, encoding="utf-8") as fh:
            text = _retarget(fh.read())
        with open(dst, "w", encoding="utf-8") as fh:
            fh.write(text)
        shutil.copymode(src, dst)

    if edits:
        for fname, fn in edits.items():
            path = os.path.join(variant, fname)
            if not os.path.isfile(path):
                raise FileNotFoundError(
                    f"edit target {fname!r} is not among the files copied into "
                    f"the variant (see CASE_FILES): {path}"
                )
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(fn(text))

    return variant


def run_variant(variant, dim=2, log_name="log_z3st.md", check=True):
    """Mesh and solve a variant in place. Returns the z3st return code."""
    env = os.environ.copy()
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
    )
    env["PYTHONPATH"] = repo_root + os.pathsep + env.get("PYTHONPATH", "")

    make_bcs = os.path.join(variant, "make_bcs.py")
    if os.path.isfile(make_bcs):
        subprocess.run(["python3", "make_bcs.py"], cwd=variant, check=check, env=env)

    subprocess.run(
        ["gmsh", "mesh.geo", f"-{dim}"],
        cwd=variant, check=check, env=env,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

    with open(os.path.join(variant, log_name), "w", encoding="utf-8") as log:
        proc = subprocess.run(
            ["python3", "-m", "z3st"],
            cwd=variant, check=check, env=env, stdout=log, stderr=log,
        )
    return proc.returncode


def clean_sweeps(case_dir):
    """Delete every variant of *case_dir*."""
    root = sweep_root(case_dir)
    name = os.path.basename(os.path.abspath(case_dir))
    if not os.path.isdir(root):
        return
    for entry in os.listdir(root):
        if entry.startswith(f"{name}__"):
            shutil.rmtree(os.path.join(root, entry))
