# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""Shared skeleton for the per-case ``non-regression.py`` scripts.

The common prologue (paths, YAML loading), the extract-mask-sort idiom, the
abs/rel arithmetic in the ``errors`` dict and the epilogue live here. Each case
keeps only what is specific to it: which fields it reads, its analytic
references, its plots.

Metrics come in three shapes, and the distinction matters because pass/fail
reads ``rel_error`` only:

    metric(value, reference)   a real comparison: abs = |v-r|, rel = |v-r|/r
    error_metric(value, rel)   the value is itself an error (an L2 norm, a
                               non-uniformity): reference 0, abs = value,
                               rel = the supplied relative norm, or value
    tracked(value)             no analytic reference; recorded for the gold
                               baseline but always passes pass_fail_check
"""

import os

import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import extract_field
from z3st.utils.utils_verification import pass_fail_check, regression_check


def case_paths(case_file):
    """``CASE_DIR, VTU_FILE, OUT_JSON`` for a case's non-regression script.

    Call as ``case_paths(__file__)``.
    """
    case_dir = os.path.dirname(os.path.abspath(case_file))
    out = os.path.join(case_dir, "output")
    return case_dir, os.path.join(out, "fields.vtu"), os.path.join(out, "non-regression.json")


def load_yaml(case_dir, name):
    """Parse one of the case's YAML files by name (e.g. ``"geometry.yaml"``)."""
    with open(os.path.join(case_dir, name)) as f:
        return yaml.safe_load(f)


def load_case(case_dir, material=None):
    """Return ``(geometry, input, material)`` for the case.

    ``material`` selects an entry of ``input.yaml::materials`` by key; the
    default takes the first entry.
    """
    geom = load_yaml(case_dir, "geometry.yaml")
    inp = load_yaml(case_dir, "input.yaml")
    mats = inp["materials"]
    mat_file = mats[material] if material is not None else next(iter(mats.values()))
    return geom, inp, load_yaml(case_dir, mat_file)


def metric(value, reference, rel=None):
    """A comparison against an analytic or reference value.

    ``rel`` overrides the relative error when the case computes it its own way
    (a relative L2 norm over a profile rather than a pointwise ratio).
    """
    value, reference = float(value), float(reference)
    abs_error = abs(value - reference)
    if rel is None:
        rel = abs_error / abs(reference) if reference else abs_error
    return {
        "numerical": value,
        "reference": reference,
        "abs_error": abs_error,
        "rel_error": float(rel),
    }


def error_metric(value, rel=None):
    """A value that is itself an error, compared against an ideal zero."""
    value = float(value)
    return {
        "numerical": value,
        "reference": 0.0,
        "abs_error": value,
        "rel_error": float(rel) if rel is not None else value,
    }


def tracked(value):
    """A quantity with no analytic reference: recorded, never failed.

    Its errors are zero by construction, so pass_fail_check always passes it,
    while regression_check still compares it against the blessed gold.
    """
    value = float(value)
    return {"numerical": value, "reference": value, "abs_error": 0.0, "rel_error": 0.0}


def line(vtu, field, component=None, **planes):
    """Extract a field along a coordinate line and return ``(axis, values)``.

    The plane is given as keyword arguments naming the *fixed* coordinates plus
    a tolerance, e.g. ``line(vtu, "Temperature", y=0.5, z=0.0, tol=1e-5)``. The
    remaining coordinate is the one the profile runs along, and both arrays come
    back sorted by it.
    """
    tol = planes.pop("tol")
    coords = dict(zip("xyz", extract_field(vtu, field_name=field)[:3]))
    values = extract_field(vtu, field_name=field)[3]

    mask = np.ones(len(values), dtype=bool)
    for axis, target in planes.items():
        mask &= np.abs(coords[axis] - target) < tol

    if not np.any(mask):
        raise RuntimeError(f"No points found for '{field}' at {planes} (tol={tol})")

    (along,) = [a for a in "xyz" if a not in planes]
    order = np.argsort(coords[along][mask])
    if component is not None:
        values = values[:, component]
    return coords[along][mask][order], values[mask][order]


def finish(errors, tolerance, out_json, case_dir):
    """Run the pass/fail check and the gold regression check, then report."""
    pass_fail_check(errors, tolerance, out_json, case_dir)
    regression_check(errors, case_dir)
    print("\n[INFO] non-regression completed.\n")
