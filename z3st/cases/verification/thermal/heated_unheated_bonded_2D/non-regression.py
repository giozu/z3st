#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/thermal/heated_unheated_bonded_2D

non-regression script
---------------------
A heated layer (uniform q''', 0 < x < a) bonded to an unheated one
(a < x < Lx) on one conforming mesh. xmin is adiabatic, xmax is held at T0,
so all the power qa leaves through the unheated layer:

    0 < x < a:   T = T(a) + q (a² - x²) / (2 k_h)
    a < x < Lx:  T = T0 + q a (Lx - x) / k_c

The solution does not depend on y and the source is constant on every heated
element, so bilinear elements reproduce it at the nodes to round-off.

What it guards: q_third is nodal, and the interface nodes carry the heated
value. The source must be integrated over the heated material only; over the
unheated one it would add a linearly decaying source across the first element
(q h / 2 per unit area, +12.5 % of the power with this mesh).
"""

import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_vtu import extract_field

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

geom, inp, mat_h = load_case(CASE_DIR, material="heated")
_, _, mat_c = load_case(CASE_DIR, material="cold")
a, Lx = float(geom["a"]), float(geom["Lx"])  # m
k_h, k_c = float(mat_h["k"]), float(mat_c["k"])  # W/m·K
q = float(mat_h["gamma_heating"])  # W/m³
T0 = 500.0  # K, Dirichlet at xmax

TOLERANCE = 1.0e-6  # - nodally exact up to the MUMPS solve


# --.. ..- .-.. .-.. --- analytic solution --.. ..- .-.. .-.. ---
def analytic_T(x):
    T_a = T0 + q * a * (Lx - a) / k_c
    return np.where(
        x >= a,
        T0 + q * a * (Lx - x) / k_c,
        T_a + q * (a**2 - x**2) / (2.0 * k_h),
    )


# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
x, y, z, T = extract_field(VTU_FILE, field_name="Temperature")
T_ref = analytic_T(x)
span = T_ref.max() - T0

dT = np.abs(T - T_ref)
print(f"[INFO] T range (analytic): {T0:.3f} .. {T_ref.max():.3f} K")
print(f"[INFO] max |T - T_ref| = {dT.max():.3e} K  (relative to the rise: {dT.max()/span:.3e})")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "Linf_error_T": error_metric(dT.max(), rel=dT.max() / span),
    # rise above T0 rather than T: a spurious source moves it by ~10 %
    "T_rise_max": metric(T.max() - T0, span),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
