#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 1_thin_slab_neumann_3D

non-regression script
---------------------
Steady-state 3D slab (Neumann-Dirichlet).

"""

import os

import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plotter_sigma_temperature_slab
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Lx, Ly, Lz = 0.100, 2.0, 2.0  # m (geometry dimensions)
k, E, nu, alpha = (
    48.1,
    1.77e11,
    0.3,
    1.7e-5,
)  # W/m·K, Pa, -, 1/K (thermal conductivity, Young's modulus, Poisson's ratio, thermal expansion)
Ti, To = 573.0, 583.0  # K (boundary temperature)
q0, mu = 0.00, 24.0  # W/m³, 1/m (volumetric heat source, attenuation coefficient)
y_target, z_target, mask_tol = Ly / 2, Lz / 2, 0.65  # m, m, m (plane selection and tolerance)

TOLERANCE = 3e-3  # - (relative tolerance for non-regression tests)


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile (slab)."""
    term1 = Ti + (To - Ti) * (x / Lx)
    term2 = (q0 / (mu**2 * k)) * ((x / Lx) * (np.exp(-mu * Lx) - 1) - (np.exp(-mu * x) - 1))
    return term1 + term2


def sigma_th(x, T_num, c=1.0):
    """Thermal stress profile (linear)."""
    T_mean = np.trapezoid(T_num, x) / (x.max() - x.min())
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
# Temperature
x_T, y_T, z_T, T_all = extract_field(VTU_FILE, field_name="Temperature")
mask = (np.abs(y_T - y_target) < mask_tol) & (np.abs(z_T - z_target) < mask_tol)
sort_idx = np.argsort(x_T[mask])

x_T = x_T[mask][sort_idx]
T = T_all[mask][sort_idx]

# Stress
x_S, y_S, z_S, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
mask = (np.abs(y_S - y_target) < mask_tol) & (np.abs(z_S - z_target) < mask_tol)
sort_idx = np.argsort(x_S[mask])

x_s = x_S[mask][sort_idx]

sigma_xx = S_all[mask, 0][sort_idx]
sigma_yy = S_all[mask, 4][sort_idx]

# Analytical results
T_ref = analytic_T(x_T)
sigma_th_ref = sigma_th(x_T, T_ref, c=1.0)
max_sigma_T = np.max(sigma_yy)

# Plot
plotter_sigma_temperature_slab(
    x_s=x_s,
    sigma=sigma_yy,
    x_s_ref=x_T,
    sigma_ref=sigma_th_ref,
    T_ref=T_ref,
    x_T=x_T,
    T=T,
    max_sigma_T=max_sigma_T,
    Ti=Ti,
    To=To,
    CASE_DIR=CASE_DIR,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref) ** 2)))
Linf_T = float(np.max(np.abs((T - T_ref))))
RelL2_T = float(L2_T / np.mean(np.abs(T_ref)))

Tmax_num = float(np.max(T))
Tmax_ref = float(np.max(T_ref))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / Tmax_ref

errors = {
    "L2_error_T": {
        "numerical": L2_T,
        "reference": 0.0,
        "abs_error": L2_T,
        "rel_error": RelL2_T,
    },
    "Linf_error_T": {
        "numerical": Linf_T,
        "reference": 0.0,
        "abs_error": Linf_T,
        "rel_error": Linf_T / np.mean(np.abs(T_ref)),
    },
    "T_max": {
        "numerical": Tmax_num,
        "reference": Tmax_ref,
        "abs_error": abs(Tmax_num - Tmax_ref),
        "rel_error": RelErr_Tmax,
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
