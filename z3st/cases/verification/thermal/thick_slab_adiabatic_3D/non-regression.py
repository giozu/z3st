#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/thermal/thick_slab_adiabatic_3D

non-regression script
---------------------
Steady-state 3D slab (Dirichlet-Neumann).

"""

import os

import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric
from z3st.utils.utils_extract_vtu import average_section, extract_stress, extract_temperature, list_fields
from z3st.utils.utils_plot import plotter_sigma_temperature_slab

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, VTU_FILE, OUT_JSON = case_paths(__file__)

# Geometry and material
geom, inp, mat = load_case(CASE_DIR)
Lx, Ly, Lz = float(geom["Lx"]), float(geom["Ly"]), float(geom["Lz"])  # m (geometry dimensions)

k, E, nu, alpha = (
    float(mat["k"]),
    float(mat["E"]),
    float(mat["nu"]),
    float(mat["alpha"]),
)  # W/m·K, Pa, -, 1/K (thermal conductivity, Young's modulus, Poisson's ratio, thermal expansion)
Ti = 490.0  # K (boundary temperature)
q0, mu = float(mat["gamma_heating"]), float(
    mat["mu_gamma"]
)  # W/m³, 1/m (volumetric heat source, attenuation coefficient)
y_target, z_target, mask_tol = Ly / 2, Lz / 2, 0.1  # m, m, m (plane selection and tolerance)

TOLERANCE = 3e-3  # - (relative tolerance for non-regression tests)


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile (slab)."""
    return Ti + q0 / (k * mu**2) * (1 - np.exp(-mu * x) - mu * x * np.exp(-mu * Lx))


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
sigma_th_ref = sigma_th(x_T, T_ref, c=1.0)
max_sigma_T = np.max(sigma_yy)

# attenuation map estimate of peak thermal stress
print(f"Ro/Ri = {2.4/2:.2f}")
print(f"mu*Ri = {mu*2:.2f}")
sigma_T_map = 0.91
sigma_th_max_map = alpha * E * q0 / ((1 - nu) * k * mu**2) * sigma_T_map
print(f"From attenuation map, maximum thermal stress = {sigma_th_max_map/1e6:.2f} MPa")


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
    To=None,
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
    "L2_error_T": error_metric(L2_T, rel=RelL2_T),
    "Linf_error_T": error_metric(Linf_T, rel=Linf_T / np.mean(np.abs(T_ref))),
    "T_max": metric(Tmax_num, Tmax_ref, rel=RelErr_Tmax),
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
