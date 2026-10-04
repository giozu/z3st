#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/thermal/two_heated_materials_2D

non-regression script
---------------------
A gamma-heated cylindrical wall split in two at r = Rm, the two parts being
the SAME material under two names. Splitting a material in two must change
nothing, so the temperature has to match, node by node, the run in single/
that meshes the same wall as one material.

What it guards: the heat source q_third is nodal, and a node on the interface
belongs to both materials. set_power() used to add each material's source
there, so the interface node carried twice the source and the split wall ran
hotter (2026-09-24, found on a liner bonded to a gamma-heated vessel).
"""

import os

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

from z3st.utils.non_regression import case_paths, error_metric, finish, tracked

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)
VTU_SINGLE = os.path.join(CASE_DIR, "single", "output", "fields.vtu")

TOLERANCE = 1.0e-6  # -   both runs solve the same linear system with MUMPS

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
split = pv.read(VTU_FILE)
single = pv.read(VTU_SINGLE)

# same geometry and transfinite divisions, so the node sets coincide; the
# ordering need not, hence the nearest-neighbour match
dist, idx = cKDTree(single.points).query(split.points)
if dist.max() > 1e-9:
    raise RuntimeError(f"node sets differ (max distance {dist.max():.2e} m)")

T_split = np.asarray(split["Temperature"], float)
T_single = np.asarray(single["Temperature"], float)[idx]

dT = np.abs(T_split - T_single)
span = T_single.max() - T_single.min()
print(f"[INFO] T range (single material): {T_single.min():.3f} .. {T_single.max():.3f} K")
print(f"[INFO] max |T_split - T_single| = {dT.max():.3e} K  (relative to the range: {dT.max()/span:.3e})")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "Linf_T_split_vs_single": error_metric(dT.max(), rel=dT.max() / span),
    "T_max": tracked(T_split.max()),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
