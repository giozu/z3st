#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 10_cylindrical_shell_thick_plane_stress

non-regression script
---------------------
Analytical non-regression for a thick-walled cylindrical shell under
internal (Pi) and external (Po) pressure.
Reference is the Lamé solution under plane stress.

"""

import os

import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plotter_sigma_cylinder, plotter_strain_cylinder
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
Ri, Ro, Lz = 0.02, 0.03, 0.01  # m          inner and outer radius, height
Pi, Po = 1.0e6, 0.0  # Pa         internal and external pressure
E, nu = 2.0e11, 0.3  # Pa, -      Young modulus, Poisson ratio
t = Ro - Ri  # m          wall thickness
slenderness = Ri / t  # -          slenderness ratio
z_target, z_tol = Lz / 2, 0.01  # m          z-plane for data extraction

TOLERANCE = 3.0e-2  # -          tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
# Lamé solutions
A = (Pi * Ri**2 - Po * Ro**2) / (Ro**2 - Ri**2)
B = (Ri**2 * Ro**2 * (Pi - Po)) / (Ro**2 - Ri**2)

eps_zz_PS = -2 * nu * A / E
sigma_zz_ana_L = 0.0  # Plane stress


def strain_rr_ref(r):
    return (1 + nu) / E * (A * (1 - 2 * nu) - B / r**2) - nu * eps_zz_PS


def strain_tt_ref(r):
    return (1 + nu) / E * (A * (1 - 2 * nu) + B / r**2) - nu * eps_zz_PS


# Mariotte solutions (Po=0)
sigma_rr_ana_M = -Pi / 2
sigma_tt_ana_M = Pi * Ri / t
sigma_zz_ana_M = Pi * Ri / (2 * t)

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
# Numerical results
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

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

# In plane stress (σ_zz = 0), ε_zz is not kinematically defined from u_z
# but must be reconstructed from the constitutive law (Hooke 3D) via σ_rr and σ_tt
strain_rr = 1 / E * (sigma_rr - nu * sigma_tt)
strain_tt = 1 / E * (sigma_tt - nu * sigma_rr)
strain_zz = -nu * (sigma_rr + sigma_tt) / E

# Analytical results
# lamé
sigma_tt_ana_L = A + B / r_s**2
sigma_zz_ana_L = sigma_zz_ana_L * np.ones_like(r_s)
sigma_rr_ana_L = A - B / r_s**2

strain_rr_ana_L = strain_rr_ref(r_s)
strain_tt_ana_L = strain_tt_ref(r_s)
strain_zz_ana_L = np.ones_like(r_s) * eps_zz_PS

# average stress
sigma_zz_ana_M = sigma_zz_ana_M * np.ones_like(r_s)

sigma_tt_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_tt * r_s, r_s)
sigma_zz_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_zz * r_s, r_s)
sigma_rr_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_rr * r_s, r_s)

# Plot
plotter_sigma_cylinder(
    r_s=r_s,
    sigma_rr=sigma_rr,
    sigma_tt=sigma_tt,
    sigma_zz=sigma_zz,
    Ri=Ri,
    Ro=Ro,
    Pi=Pi,
    Po=Po,
    CASE_DIR=CASE_DIR,
    slenderness=slenderness,
    sigma_tt_ana_L=sigma_tt_ana_L,
    sigma_zz_ana_L=sigma_zz_ana_L,
    sigma_rr_ana_L=sigma_rr_ana_L,
    # sigma_tt_ana_M=sigma_tt_ana_M,
    # sigma_zz_ana_M=sigma_zz_ana_M,
    # sigma_rr_ana_M=sigma_rr_ana_M,
    # sigma_tt_avg=sigma_tt_avg,
    # sigma_zz_avg=sigma_zz_avg,
    # sigma_rr_avg=sigma_rr_avg,
)

plotter_strain_cylinder(
    r_s=r_s,
    strain_rr=strain_rr,
    strain_tt=strain_tt,
    strain_zz=strain_zz,
    Ri=Ri,
    Ro=Ro,
    Pi=Pi,
    Po=Po,
    CASE_DIR=CASE_DIR,
    slenderness=slenderness,
    strain_rr_ana_L=strain_rr_ana_L,
    strain_tt_ana_L=strain_tt_ana_L,
    strain_zz_ana_L=strain_zz_ana_L,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_rr = np.sqrt(np.mean((sigma_rr - sigma_rr_ana_L) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_L**2))
err_tt = np.sqrt(np.mean((sigma_tt - sigma_tt_ana_L) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_L**2))
err_zz = np.sqrt(np.mean((sigma_zz - sigma_zz_ana_L) ** 2)) / np.sqrt(np.mean(sigma_tt**2))

err_eps_rr = np.sqrt(np.mean((strain_rr - strain_rr_ana_L) ** 2)) / np.sqrt(
    np.mean(strain_rr_ana_L**2)
)
err_eps_tt = np.sqrt(np.mean((strain_tt - strain_tt_ana_L) ** 2)) / np.sqrt(
    np.mean(strain_tt_ana_L**2)
)
err_eps_zz = np.sqrt(np.mean((strain_zz - strain_zz_ana_L) ** 2))

errors = {
    "L2_error_sigma_rr": {
        "numerical": float(err_rr),
        "reference": 0.0,
        "abs_error": float(err_rr),
        "rel_error": float(err_rr),
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
    "L2_error_strain_rr": {
        "numerical": float(err_eps_rr),
        "reference": 0.0,
        "abs_error": float(err_eps_rr),
        "rel_error": float(err_eps_rr),
    },
    "L2_error_strain_tt": {
        "numerical": float(err_eps_tt),
        "reference": 0.0,
        "abs_error": float(err_eps_tt),
        "rel_error": float(err_eps_tt),
    },
    "L2_error_strain_zz": {
        "numerical": float(err_eps_zz),
        "reference": 0.0,
        "abs_error": float(err_eps_zz),
        "rel_error": float(err_eps_zz),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] Non-regression completed.\n")
