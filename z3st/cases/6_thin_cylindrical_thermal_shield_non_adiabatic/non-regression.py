#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 6_thin_cylindrical_thermal_shield_non_adiabatic

non-regression script
---------------------
Steady-state 1D axisymmetric cylindrical wall (Dirichlet-Dirichlet).

"""

import os

import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plotter_sigma_temperature_cylinder
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Ri, Ro, Lz = 2.0, 2.1, 10  # m          inner and outer radius
Pi, Po = 0.0, 0.0  # Pa         internal and external pressure
k, E, nu, alpha = (
    48.1,
    1.77e11,
    0.3,
    1.7e-5,
)  # W/m·K, Pa, -, 1/K (thermal conductivity, Young's modulus, Poisson's ratio, thermal expansion)
Ti, To = 490.0, 500  # K          inner and outer surface temperature
q0, mu = 2.0e6, 24.0  # W/m³, 1/m  heat source, attenuation
Lx = Ro - Ri  # m          wall thickness
slenderness = Ri / Lx  # -          slenderness ratio
z_target, z_tol = Lz / 2, Lz / 10  # m          z-plane for data extraction

TOLERANCE = 7.0e-1  # -          tolerance for non-regression


# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(x):
    """Analytical temperature profile (slab), Dirichlet-Dirichlet."""
    x = x - Ri
    term1 = Ti + (To - Ti) * (x / Lx)
    term2 = (q0 / (mu**2 * k)) * ((x / Lx) * (np.exp(-mu * Lx) - 1) - (np.exp(-mu * x) - 1))
    return term1 + term2


def sigma_th(r, T_num, c=1.0):
    """Thermal stress profile (axisymmetric)."""
    T_mean = 2 / (Ro**2 - Ri**2) * np.trapezoid(T_num * r, r)
    return alpha * E / (1.0 - c * nu) * (T_mean - T_num)


def analytical_thermal_stress(r):
    """
    Analytical thermal stresses in a cylinder (no mechanical load)
    from integral formulation, valid for arbitrary T(r).

    Parameters
    ----------
    r : np.ndarray
        Radial coordinates (m)

    Returns
    -------
    sigma_r, sigma_t, sigma_z : np.ndarray
        Radial, hoop, and axial stresses (Pa)
    """

    r_star = np.linspace(Ri, Ro, 500)
    T_star = analytic_T(r_star)

    # Pre-compute global integral ∫_{Ri}^{Ro} T(r*) r* dr*
    I_global = np.trapezoid(T_star * r_star, r_star)

    sigma_r = np.zeros_like(r)
    sigma_t = np.zeros_like(r)
    sigma_z = np.zeros_like(r)

    for i, ri in enumerate(r):

        # Local integral ∫_{Ri}^{r} T(r*) r* dr*
        mask = r_star <= ri
        I_local = np.trapezoid(T_star[mask] * r_star[mask], r_star[mask])

        # Stress components
        sigma_r[i] = (
            (alpha * E / (1 - nu))
            * 1
            / ri**2
            * (((ri**2 - Ri**2) / (Ro**2 - Ri**2)) * I_global - I_local)
        )

        sigma_t[i] = (alpha * E / (1 - nu)) * (
            ((1 + (Ri / ri) ** 2) / (Ro**2 - Ri**2)) * I_global + I_local / ri**2 - analytic_T(ri)
        )

        sigma_z[i] = (alpha * E / (1 - nu)) * ((2 / (Ro**2 - Ri**2)) * I_global - analytic_T(ri))

    return sigma_r, sigma_t, sigma_z


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
x_T, y_T, z_T, T_all = extract_temperature(VTU_FILE)
r_T, T = average_section_radial(x_T, y_T, z_T, T_all, z_target=z_target, tol=z_tol, decimals=6)

r_s, sigma_rr, sigma_tt, sigma_zz = extract_cylindrical_stresses(
    filename=VTU_FILE,
    z_fixed=z_target,
    tol=z_tol,
    case_dir=CASE_DIR,
    stress_field_hint="Stress",
    data_source="point",
    average=True,
    decimals=6,
)

# Analytical results
sigma_th_ref = sigma_th(r_T, T, c=1.0)
T_ref = analytic_T(r_T)
sigma_rr_ana_th, sigma_tt_ana_th, sigma_zz_ana_th = analytical_thermal_stress(r_s)
max_sigma_T = np.max(sigma_tt)

# map (adiabatic, conservative reference)
print(f"Ro/Ri = {Ro/Ri:.2f}")
print(f"mu*Ri = {mu*Ri:.2f}")
sigma_T_map = 0.52
sigma_th_max_map = alpha * E * q0 / ((1 - nu) * k * mu**2) * sigma_T_map
print(f"From attenuation map, maximum thermal stress = {sigma_th_max_map/1e6:.2f} MPa")

# Plot
plotter_sigma_temperature_cylinder(
    r_s=r_s,
    sigma_tt=sigma_tt,
    sigma_zz=sigma_zz,
    r_T=r_T,
    T=T,
    T_ref=T_ref,
    label_T="Temperature (analytical, slab)",
    max_sigma_T=max_sigma_T,
    Ti=Ti,
    To=To,
    Ri=Ri,
    CASE_DIR=CASE_DIR,
    slenderness=slenderness,
    sigma_zz_ref=sigma_zz_ana_th,
    sigma_tt_ref=sigma_tt_ana_th,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_tt = np.sqrt(np.mean((sigma_tt - sigma_tt_ana_th) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_th**2))
err_zz = np.sqrt(np.mean((sigma_zz - sigma_zz_ana_th) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_th**2))

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
    "L2_error_sigma_tt": {
        "numerical": float(err_tt),
        "reference": 0.0,
        "abs_error": float(err_tt),
        "rel_error": float(err_tt),
    },
    "L2_error_sigma_zz": {
        "numerical": float(err_zz),
        "reference": 0.0,
        "abs_error": float(err_zz),
        "rel_error": float(err_zz),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
