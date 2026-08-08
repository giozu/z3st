#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: benchmarks/damage/plate_thermal_shock_2D

Kamagate et al. (2025) quenched-plate replica: crack nucleation under thermal
shock, no pre-crack and no damage seed.

What this script deliberately does NOT check
--.--.--.--.--.--.--.--.--.--.--.--.--.--.-

The **number** of nucleated cracks. README.md says it plainly: the count is
sensitive to mesh, lc and heterogeneity -- a known feature of this benchmark --
and the shipped mesh is lc/3, "for a first look", not the lc/4.5 at which counts
stop moving under refinement. A crack-count metric here would be a coin toss
dressed up as a check. The sibling pellet_quench_2D_xy shows what that leads to:
a reference of 3 against a measured 5, kept green by a forced "pass": True.

What it checks instead: quantities that do not depend on how a damage band
happens to split. The temperature field is effectively one-way (constant k and c,
damage heat neglected), so it is the most reproducible thing the case produces.
Everything is tracked() -- recorded in the gold and guarded against regression --
because the paper gives a pattern to match, not numbers to hit.
"""

import os

import numpy as np

from z3st.utils.non_regression import finish, tracked
from z3st.utils.utils_extract_xdmf import extract_field_xdmf, list_fields_xdmf

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
XDMF_FILE = os.path.join(OUT_DIR, "fields.xdmf")
OUT_JSON = os.path.join(OUT_DIR, "non-regression.json")

TOLERANCE = 1e-2

if not os.path.exists(XDMF_FILE):
    raise FileNotFoundError(f"XDMF output not found: {XDMF_FILE} (run ./Allrun first)")

print(f"[INFO] Reading {os.path.basename(XDMF_FILE)}")
list_fields_xdmf(XDMF_FILE)

errors = {}

# --.. ..- .-.. .-.. --- temperature, final step --.. ..- .-.. .-.. ---
_, _, _, T = extract_field_xdmf(XDMF_FILE, "Temperature", step_index=-1)
T = np.asarray(T).reshape(-1)
print(f"[INFO] T final: min = {T.min():.2f} K, max = {T.max():.2f} K, mean = {T.mean():.2f} K")
errors["T_min_final"] = tracked(T.min())
errors["T_max_final"] = tracked(T.max())
errors["T_mean_final"] = tracked(T.mean())

# --.. ..- .-.. .-.. --- damage, final step --.. ..- .-.. .-.. ---
# D_max is a threshold, not a count: it says whether nucleation happened at all.
# The damaged-node fraction integrates "how much" damage there is without caring
# how it splits into bands. It is mesh-dependent in absolute terms, which is
# harmless here: the regression check compares runs on the same mesh, and a mesh
# change is exactly the event that should force a re-blessing.
_, _, _, D = extract_field_xdmf(XDMF_FILE, "Damage", step_index=-1)
D = np.asarray(D).reshape(-1)
frac = float(np.count_nonzero(D > 0.5)) / D.size
print(f"[INFO] D final: max = {D.max():.4f}, nodes with D > 0.5 = {frac*100:.3f} %")
errors["D_max_final"] = tracked(D.max())
errors["D_mean_final"] = tracked(D.mean())
errors["damaged_node_fraction"] = tracked(frac)

# --.. ..- .-.. .-.. --- global energy balance, if the solver wrote it --..
energy_file = os.path.join(CASE_DIR, "energies.txt")
if os.path.exists(energy_file):
    data = np.genfromtxt(energy_file, names=True)
    print(f"[INFO] energies: E_el = {data['E_el'][-1]:.4e} J, E_frac = {data['E_frac'][-1]:.4e} J")
    errors["E_el_final"] = tracked(data["E_el"][-1])
    errors["E_frac_final"] = tracked(data["E_frac"][-1])
    errors["E_tot_final"] = tracked(data["E_tot"][-1])
else:
    print(f"[WARN] {os.path.basename(energy_file)} not found; skipping energy metrics.")

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
