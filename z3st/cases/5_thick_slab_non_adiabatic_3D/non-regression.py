#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 5_thick_slab_non_adiabatic_3D

non-regression script
---------------------
Steady-state 3D slab (Dirichlet-Dirichlet).

"""

import os

import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plotter_sigma_temperature_slab
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material (loaded from YAML)
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])  # m (geometry dimensions)

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
mat_path = os.path.join(CASE_DIR, next(iter(inp["materials"].values())))
with open(mat_path) as f:
    mat = yaml.safe_load(f)
k, E, nu, alpha = (
    float(mat["k"]),
    float(mat["E"]),
    float(mat["nu"]),
    float(mat["alpha"]),
)  # W/m·K, Pa, -, 1/K (thermal conductivity, Young's modulus, Poisson's ratio, thermal expansion)
q0, mu = float(mat["gamma_heating"]), float(
    mat["mu_gamma"]
)  # W/m³, 1/m (volumetric heat source, attenuation coefficient)
Ti, To = 490, 500.0  # K (boundary temperature)
y_target, z_target, mask_tol = Ly / 2, Lz / 2, 0.1  # m, m, m (plane selection and tolerance)

TOLERANCE = 3e-2  # - (relative tolerance for non-regression tests)


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
x_T, y_T, z_T, T_all = extract_temperature(VTU_FILE)
x_T, T = average_section(x_T, y_T, z_T, T_all, y_target, z_target, mask_tol, label="T", decimals=5)

x_s, y_s, z_s, s = extract_stress(VTU_FILE, component="all", return_coords=True, prefer="cells")
x_s, sigma_yy = average_section(
    x_s, y_s, z_s, s["yy"], y_target, z_target, mask_tol, decimals=5, label="sigma_yy"
)

# Analytical results
T_ref = analytic_T(x_T)
sigma_th_ref = sigma_th(x_s, analytic_T(x_s), c=1.0)
max_sigma_T = np.max(sigma_yy)

# Plot
plotter_sigma_temperature_slab(
    x_s=x_s,
    sigma=sigma_yy,
    x_s_ref=x_s,
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

err_sigma = np.sqrt(np.mean((sigma_yy - sigma_th_ref) ** 2)) / np.max(np.abs(sigma_th_ref))

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
    "L2_error_sigma_yy": {
        "numerical": float(err_sigma),
        "reference": 0.0,
        "abs_error": float(err_sigma),
        "rel_error": float(err_sigma),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
