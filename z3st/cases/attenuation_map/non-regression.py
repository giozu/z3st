#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: attenuation_map

non-regression script
---------------------

"""

import os
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *
from z3st.utils.utils_plot import plotter_sigma_temperature_cylinder

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Ri, Ro, Lz = 1.0, 1.10, 10            # m          inner and outer radius
Pi, Po = 0.0, 0.0                     # Pa         internal and external pressure
k, E, nu, alpha = 18.5, 1.8e11, 0.28, 1.4e-5   # W/m·K, Pa, -, 1/K
Ti = 490.0                            # K          inner surface temperature
q0, mu = 9.0e4, 27.0                  # W/m³, 1/m  heat source, attenuation
Lx = Ro - Ri                          # m          wall thickness
slenderness = Ri / Lx                 # -          slenderness ratio
z_target, z_tol = Lz/2, 1               # m          z-plane for data extraction

TOLERANCE = 5.0e-2                    # -          tolerance for non-regression


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile (slab), Dirichlet-Neumann."""
    return Ti + q0 / (k * mu**2) * (1 - np.exp(-mu * x) - mu * x * np.exp(-mu * Lx))

def sigma_th(r, T_num, c=1.0):
    """Thermal stress profile (axisymmetric)."""
    T_mean = 2 / (Ro**2 - Ri**2) * np.trapezoid(T_num * r, r)
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
# Numerical results
x_T, y_T, z_T, T_all = extract_temperature(VTU_FILE)
r_T, T = average_section_radial(x_T, y_T, z_T, T_all, z_target=z_target, tol=z_tol, decimals=6)

r_s, sigma_rr, sigma_tt, sigma_zz = extract_cylindrical_stresses(
    filename=VTU_FILE,
    z_fixed=z_target,
    tol=z_tol,
    case_dir=CASE_DIR,
    stress_field_hint="Stress",
    data_source="auto",
    average=True,
    decimals=6,
)

# Analytical results
sigma_th_ref = sigma_th(r_T, T, c=1.0)
T_ref = analytic_T(r_T - Ri)
max_sigma_T = np.max(sigma_tt)

# map
print(f"Ro/Ri = {Ro/Ri:.2f}")
print(f"mu*Ri = {mu*Ri:.2f}") 
sigma_T_map = 0.52
sigma_th_max_map = alpha * E * q0 / ((1-nu)*k*mu**2) * sigma_T_map
print(f"From attenuation map, maximum thermal stress = {sigma_th_max_map/1e6:.2f} MPa")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T - T_ref)**2)))
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
    "sigma_th_tt": {
        "numerical": max_sigma_T,
        "rel_error": 0.0
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
