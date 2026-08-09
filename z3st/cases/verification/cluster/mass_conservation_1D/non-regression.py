#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/cluster/mass_conservation_1D

1-D cluster-dynamics transport (implicit-Euler DG1, upwind interior-facet flux,
SIPG diffusion). Checks the one property the model asserts by construction: the
total defect mass C_tot = integral of c(n) * n dn is conserved, defects changing
size rather than disappearing.

Source of the numbers
--.--.--.--.--.--.--.

C_tot is a weighted integral of a discontinuous P1 field: two dofs per cell, so a
trapezoid over sorted nodal coordinates double-counts every interface. The values
read here are the ones the solver assembles in UFL and prints each step.

What is checked
--.--.--.--.--.

The model enforces conservation: after each step it rescales c by target /
current. The quantity carrying information is the rescale factor: 1 - factor is
the mass the discretisation loses in one step before the correction. On the
shipped configuration the worst step needs a 3.6 % correction, a mixture of DG
transport error and outflow through the open boundary (velocity 1.0, no BC on the
cluster field). The script pins that value; it does not decide whether the split
is acceptable.
"""

import os
import re

import numpy as np

from z3st.utils.non_regression import case_paths, finish, tracked
from z3st.utils.utils_extract_xdmf import extract_field_xdmf

CASE_DIR, _, OUT_JSON = case_paths(__file__)
XDMF_FILE = os.path.join(CASE_DIR, "output", "results.xdmf")
SOLVER_LOG = os.path.join(CASE_DIR, "log_z3st.md")

TOLERANCE = 1e-2

# --.. ..- .-.. .-.. --- mass bookkeeping from the solver log --.. ..- .-.. .-.. ---
if not os.path.exists(SOLVER_LOG):
    raise FileNotFoundError(f"solver log not found: {SOLVER_LOG} (run ./Allrun first)")
log = open(SOLVER_LOG).read()

targets = [float(v) for v in re.findall(r"Mass conservation: target = ([0-9.eE+-]+)", log)]
before = [float(v) for v in re.findall(r"Mass conservation: before = ([0-9.eE+-]+)", log)]
factors = [float(v) for v in re.findall(r"Mass conservation: factor = ([0-9.eE+-]+)", log)]

if not (targets and before and factors):
    raise RuntimeError("no mass-conservation lines in the solver log; did the cluster model run?")

C_tot_target = targets[-1]
worst_deviation = float(np.max(np.abs(np.array(factors) - 1.0)))
worst_step = int(np.argmax(np.abs(np.array(factors) - 1.0)))

print(f"[INFO] steps with mass bookkeeping : {len(factors)}")
print(f"[INFO] conserved C_tot             : {C_tot_target:.6e}")
print(f"[INFO] pre-rescale mass, last step : {before[-1]:.6e}")
print(f"[INFO] worst rescale factor        : {factors[worst_step]:.8f} at step {worst_step} "
      f"({worst_deviation * 100:.2f} % correction)")

errors = {
    # The target is the physical invariant: it must not drift between runs.
    "C_tot_target": tracked(C_tot_target),
    # Pre-rescale mass on the final step: 1.88 % from the invariant, worst single
    # step 3.64 %. tracked(), not compared: the acceptable mass loss is a physics
    # decision about the DG transport and the open boundary.
    "C_tot_before_rescale_final": tracked(before[-1]),
    "rel_gap_before_rescale": tracked(abs(before[-1] - C_tot_target) / C_tot_target),
    # Worst single-step correction over the run.
    "worst_rescale_deviation": tracked(worst_deviation),
}

# --.. ..- .-.. .-.. --- peak density, from the field itself --.. ..- .-.. .-.. ---
# Cross-checks that the written output is consistent with what the solver reported.
if os.path.exists(XDMF_FILE):
    _, _, _, c = extract_field_xdmf(XDMF_FILE, "ClusterDensity", step_index=-1)
    errors["c_max_final"] = tracked(np.max(np.asarray(c)))
    print(f"[INFO] peak density, final step    : {np.max(np.asarray(c)):.4f}")
else:
    print(f"[WARN] {os.path.basename(XDMF_FILE)} not found; skipping the field cross-check.")

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
