#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 17_annular_cylinder

non-regression script
-----------------------
Analytical solution for radial temperature in an
annular cylinder with uniform volumetric heat generation.

"""

import os
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *
from z3st.utils.utils_plot import plotter_sigma_temperature_cylinder, plotter_sigma_cylinder

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Ri, Ro, Lz = 0.004, 0.04, 0.10        # m          inner and outer radius
Pi, Po = 0.0, 0.0                     # Pa         internal and external pressure
k, E, nu, alpha = 2.5, 1.7e11, 0.3, 1.45e-5    # W/m·K, Pa, -, 1/K (thermal conductivity, Young's modulus, Poisson's ratio, thermal expansion)
To = 500.0                            # K          outer surface temperature
Lx = Ro - Ri                          # m          wall thickness
slenderness = Ri / Lx                 # -          slenderness ratio
z_target, z_tol = 0.0, Lz/20          # m          z-plane for data extraction
LHR = 5.0e+2                          # linear heat generation (W/m)

TOLERANCE = 7.0e-1                    # -          tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
def analytic_T(r):
    """Analytical temperature profile in an annular cylinder with volumetric heat generation."""
    return To + LHR / (4 * k * np.pi) * ((Ro**2 - r**2) / (Ro**2 - Ri**2) + 2*Ri**2*np.log(r / Ro) / (Ro**2 - Ri**2))

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

    r_star = np.linspace(Ri, Ro, 20000)
    T_star = analytic_T(r_star)

    # Pre-compute global integral ∫_{Ri}^{Ro} T(r*) r* dr*
    I_global = np.trapezoid(T_star * r_star, r_star)

    sigma_r = np.zeros_like(r)
    sigma_t = np.zeros_like(r)
    sigma_z = np.zeros_like(r)

    for i, ri in enumerate(r):

        # Local integral ∫_{Ri}^{r} T(r*) r* dr*
        mask = r_star <= ri
        I_local  = np.trapezoid(T_star[mask] * r_star[mask], r_star[mask])

        # Stress components
        sigma_r[i] = (alpha * E / (1 - nu)) * 1/ri**2 * (
            ((ri**2 - Ri**2) / (Ro**2 - Ri**2)) * I_global
            - I_local
        )

        sigma_t[i] = (alpha * E / (1 - nu)) * (
            ((1 + (Ri/ri)**2) / (Ro**2 - Ri**2)) * I_global
            + I_local / ri**2
            - analytic_T(ri)
        )

        sigma_z[i] = (alpha * E / (1 - nu)) * (
            (2 / (Ro**2 - Ri**2)) * I_global
            - analytic_T(ri)
        )

    return sigma_r, sigma_t, sigma_z

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

x, y, z, T = extract_temperature(VTU_FILE)

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
T_ref = analytic_T(r_T)
sigma_rr_ana_th, sigma_tt_ana_th, sigma_zz_ana_th = analytical_thermal_stress(r_s)
max_sigma_T = np.max(sigma_tt)

# Plot
plotter_sigma_temperature_cylinder(
    r_s=r_s,
    sigma_rr=sigma_rr,
    sigma_tt=sigma_tt,
    sigma_zz=sigma_zz,
    r_T=r_T,
    T=T,
    T_ref=T_ref,
    max_sigma_T=max_sigma_T,
    Ti=None,
    To=To,
    Ri=Ri,
    CASE_DIR=CASE_DIR,
    slenderness=slenderness,
    sigma_rr_ref=sigma_rr_ana_th,
    sigma_zz_ref=sigma_zz_ana_th,
    sigma_tt_ref=sigma_tt_ana_th,
)

# Plot
# plotter_sigma_cylinder(
#     r_s=r_s,
#     sigma_rr=sigma_rr,
#     sigma_tt=sigma_tt,
#     sigma_zz=sigma_zz,
#     Ri=Ri, Ro=Ro, Pi=1e6, Po=Po, # Pi set to 1 MPa for visualization
#     CASE_DIR=CASE_DIR,
#     slenderness=slenderness,
# )

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_tt = np.sqrt(np.mean((sigma_tt - sigma_tt_ana_th)**2)) / np.sqrt(np.mean(sigma_tt_ana_th**2))
err_rr = np.sqrt(np.mean((sigma_rr - sigma_rr_ana_th)**2)) / np.sqrt(np.mean(sigma_rr_ana_th**2))
err_zz = np.sqrt(np.mean((sigma_zz - sigma_zz_ana_th)**2)) / np.sqrt(np.mean(sigma_zz_ana_th**2))

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
    "L2_error_sigma_rr": {
        "numerical": float(err_rr),
        "reference": 0.0,
        "abs_error": float(err_rr),
        "rel_error": float(err_rr)
    },
    "L2_error_sigma_tt": {
        "numerical": float(err_tt),
        "reference": 0.0,
        "abs_error": float(err_tt),
        "rel_error": float(err_tt)
    },
    "L2_error_sigma_zz": {
        "numerical": float(err_zz),
        "reference": 0.0,
        "abs_error": float(err_zz),
        "rel_error": float(err_zz)
    }
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
