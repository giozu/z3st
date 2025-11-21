#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 14_thick_cylindrical_thermal_shield_adiabatic

non-regression script
---------------------
Steady-state 1D axisymmetric cylindrical wall (Dirichlet-Neumann).

"""

import os
import numpy as np

from utils.utils_extract_vtu import *
from utils.utils_verification import *
from utils.utils_plot import plotter_sigma_temperature_cylinder

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Ri, Ro = 2.0, 2.40                    # m          inner and outer radius
Pi, Po = 0.0, 0.0                     # Pa         internal and external pressure
k, E, nu, alpha = 48.1, 1.77e11, 0.3, 1.7e-5   # W/m·K, Pa, -, 1/K
Ti = 490.0                            # K          inner surface temperature
q0, mu = 2.0e6, 24.0                  # W/m³, 1/m  heat source, attenuation
Lx = Ro - Ri                          # m          wall thickness
slenderness = Ri / Lx                 # -          slenderness ratio
z_target, z_tol = 10, 1               # m          z-plane for data extraction

print(0.5 * alpha * E * q0 / (k * mu**2 * (1-nu))/1e6)

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
    data_source="cell",
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
sigma_T_map = 0.91
sigma_th_max_map = alpha * E * q0 / ((1-nu)*k*mu**2) * sigma_T_map
print(f"From attenuation map, maximum thermal stress = {sigma_th_max_map/1e6:.2f} MPa")

# Plot
plotter_sigma_temperature_cylinder(
    r_s=r_s,
    sigma_rr=sigma_rr,
    sigma_tt=sigma_tt,
    sigma_zz=sigma_zz,
    r_T=r_T,
    T=T,
    sigma_th_ref=sigma_th_ref,
    T_ref=T_ref,
    max_sigma_T=max_sigma_T,
    Ti=Ti,
    To=None,
    Ri=Ri,
    CASE_DIR=CASE_DIR,
    slenderness=slenderness,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
L2 = float(np.sqrt(np.mean((T - T_ref)**2)))
Linf = float(np.max(np.abs((T - T_ref))))
RelL2 = float(L2 / np.mean(np.abs(T_ref)))

Tmax_num = float(np.max(T))
Tmax_ref = float(np.max(T_ref))
RelErr_Tmax = abs(Tmax_num - Tmax_ref) / Tmax_ref

errors = {
    "L2_error_T": {
        "numerical": L2,
        "reference": 0.0,
        "abs_error": L2,
        "rel_error": RelL2,
    },
    "Linf_error_T": {
        "numerical": Linf,
        "reference": 0.0,
        "abs_error": Linf,
        "rel_error": Linf / np.mean(np.abs(T_ref)),
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

print("\n[INFO] Non-regression completed.\n")
