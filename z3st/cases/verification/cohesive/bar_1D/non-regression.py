#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Non-regression for verification/cohesive/bar_1D.

Verifies the cohesive phase-field model against the closed-form solution of
the 1D bar in tension, Vicentini et al. (2026) Sec. 3.6.

Every equilibrium state of the localized solution satisfies Eq. (85),

    U_t E / (2 L sigma_c) = a0 / ((1 - a0) B) + (1 - a0)^2,
    B = 2 L / ell_ch,      ell_ch = Gc E / sigma_c^2,

with the control parameter a0 = alpha(0) fixed by the stress through Eq. (69),
sigma / sigma_c = (1 - a0)^2. This holds branch by branch, so it is a valid
check even where the structural response snaps back (B = 5 here) and the
displacement-controlled solution jumps between branches.

The comparison therefore reads a0 off the computed stress and asks whether the
computed displacement is the one the theory assigns to it.
"""

import os
import re

import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric, tracked

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, _VTU, OUT_JSON = case_paths(__file__)
CSV = os.path.join(CASE_DIR, "output", "response.csv")

# Discretizing a displacement jump with continuous linear elements overestimates
# the dissipated energy. Eq. (149) gives the overestimate in closed form,
#
#     D_h - D = (Gc / c_w) w(a0) (h / ell),
#
# which against D = Gc a0^2 is the fixed relative factor h / (2 ell) for AT2.
# The comparison below therefore uses the fracture energy the discretization
# actually delivers rather than the continuum Gc; the raw deviation is kept as
# a tracked quantity so the bias itself stays visible and is checked to be the
# size Eq. (149) predicts.
TOLERANCE = 0.02

geom, inp, mat = load_case(CASE_DIR)

L2 = float(geom["Lx"])              # 2 L, the full bar length
E = float(mat["E"])
Gc = float(mat["Gc"])
sigma_c = float(mat["p_c"])

ell = float(inp["models"]["cohesive"]["ell"])
n_elements = int(re.search(r"nx\s*=\s*(\d+)", open(os.path.join(CASE_DIR, "mesh.geo")).read()).group(1)) - 1
h = L2 / n_elements

Gc_h = Gc * (1.0 + h / (2.0 * ell))   # Eq. (149), AT2
ell_ch = Gc_h * E / sigma_c**2
B = L2 / ell_ch
U_peak = L2 * sigma_c / E             # displacement at the elastic limit

print(f"[INFO] h = {h:.4e} m, h/ell = {h / ell:.4f}, Gc_h/Gc = {Gc_h / Gc:.4f}")
print(f"[INFO] ell_ch = {ell_ch:.4e} m, B = {B:.4f}, U_peak = {U_peak:.4e} m")

# --.. ..- .-.. .-.. --- numerical response --.. ..- .-.. .-.. ---
data = np.genfromtxt(CSV, delimiter=",", names=True)
U = np.atleast_1d(data["U_t_m"])
sigma = np.atleast_1d(data["sigma_xx_Pa"])
alpha_max = np.atleast_1d(data["alpha_max"])
E_frac = np.atleast_1d(data["E_frac_J"])

# --.. ..- .-.. .-.. --- peak stress --.. ..- .-.. .-.. ---
# Nucleation is at sigma_c independently of ell: that independence is the point
# of the model, so the peak is a direct check on the prescribed strength.
sigma_peak = float(sigma.max())

# --.. ..- .-.. .-.. --- softening branch --.. ..- .-.. .-.. ---
# Post-peak states only: while the bar is still elastic a0 = 0 and Eq. (85)
# degenerates to the elastic line, which tests nothing.
localized = alpha_max > 1e-6
if not np.any(localized):
    raise RuntimeError("No localized state reached; the bar never cracked.")

a0 = 1.0 - np.sqrt(np.clip(sigma[localized] / sigma_c, 0.0, 1.0))
U_theory = U_peak * (a0 / ((1.0 - a0) * B) + (1.0 - a0) ** 2)
U_num = U[localized]

# The same curve against the continuum Gc, to keep the discretization bias
# itself under observation rather than merely corrected away.
B_raw = L2 * sigma_c**2 / (Gc * E)
U_raw = U_peak * (a0 / ((1.0 - a0) * B_raw) + (1.0 - a0) ** 2)
bias = float(np.mean(np.abs(U_num - U_raw) / U_raw))
print(f"[INFO] mean deviation from the continuum curve = {bias:.4%} "
      f"(Eq. 149 predicts a fracture-energy bias of {h / (2 * ell):.4%})")

rel_dev = np.abs(U_num - U_theory) / np.abs(U_theory)
print(f"[INFO] {localized.sum()} localized states, "
      f"max |dU|/U = {rel_dev.max():.4e}, mean = {rel_dev.mean():.4e}")

errors = {
    "peak_stress_Pa": metric(sigma_peak, sigma_c),
    "softening_branch_U_rel_L2": error_metric(
        float(np.sqrt(np.mean(rel_dev**2))), rel=float(np.sqrt(np.mean(rel_dev**2)))
    ),
    "softening_branch_U_rel_max": error_metric(float(rel_dev.max()), rel=float(rel_dev.max())),
    # Eq. (69): the stress on the softening branch is a(alpha(0)) sigma_c, so
    # the damage read off the stress must be the damage the solver reports.
    "alpha_from_stress": metric(
        float(1.0 - np.sqrt(np.clip(sigma[-1] / sigma_c, 0.0, 1.0))),
        float(alpha_max[-1]),
    ),
    # The bias is not merely corrected away: it must be the size Eq. (149)
    # predicts, which is what confirms the correction is the right one.
    "discretization_bias_vs_eq149": metric(bias, h / (2.0 * ell)),
    "alpha_max_final": tracked(float(alpha_max[-1])),
    "fracture_energy_final_J": tracked(float(E_frac[-1])),
    "stress_final_Pa": tracked(float(sigma[-1])),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
