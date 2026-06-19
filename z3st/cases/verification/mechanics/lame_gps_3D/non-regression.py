#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/mechanics/lame_gps_3D

3D thick-walled cylindrical shell under internal (Pi) and external (Po)
pressure. Reference is the Lamé solution under generalized plane strain.
"""

import os

import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import plotter_sigma_cylinder, plotter_strain_cylinder
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Geometry and material
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ri, Ro, Lz = float(geom["Ri"]), float(geom["Ro"]), float(geom["Lz"])  # m  inner/outer radius, height

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
mat_path = os.path.join(CASE_DIR, next(iter(inp["materials"].values())))
with open(mat_path) as f:
    mat = yaml.safe_load(f)

Pi, Po = 1.0e6, 0.0  # Pa         internal and external pressure
E, nu = float(mat["E"]), float(mat["nu"])  # Pa, -      Young modulus, Poisson ratio
t = Ro - Ri  # m          wall thickness
slenderness = Ri / t  # -          slenderness ratio
z_target, z_tol = Lz / 2, 0.01  # m          z-plane for data extraction

TOLERANCE = 3.0e-1  # -          tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
eps_zz_GPS = -2.718310e-06

# Lamé solutions
A = (Pi * Ri**2 - Po * Ro**2) / (Ro**2 - Ri**2)
B = (Ri**2 * Ro**2 * (Pi - Po)) / (Ro**2 - Ri**2)
sigma_zz_ana_L = 2 * nu * A + E * eps_zz_GPS  # Generalized plane strain (epsilon_z = const)


def epsilon_rr_ref(r):
    return (1 + nu) / E * (A * (1 - 2 * nu) - B / r**2) - nu * eps_zz_GPS


def epsilon_tt_ref(r):
    return (1 + nu) / E * (A * (1 - 2 * nu) + B / r**2) - nu * eps_zz_GPS


# Mariotte solutions (Po=0)
sigma_rr_ana_M = -Pi / 2
sigma_tt_ana_M = Pi * Ri / t
sigma_zz_ana_M = Pi * Ri / (2 * t)

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
r_s, sigma_rr, sigma_tt, sigma_zz = extract_cylindrical_field(
    filename=VTU_FILE,
    z_fixed=z_target,
    tol=z_tol,
    case_dir=CASE_DIR,
    field_hint="Stress",
    data_source="cell",
    average=True,
    decimals=6,
)

r_s, epsilon_rr, epsilon_tt, epsilon_zz = extract_cylindrical_field(
    filename=VTU_FILE,
    z_fixed=z_target,
    tol=z_tol,
    case_dir=CASE_DIR,
    field_hint="Strain",
    data_source="cell",
    average=True,
    decimals=6,
)

# Analytical results
sigma_tt_ana_L = A + B / r_s**2
sigma_zz_ana_L = sigma_zz_ana_L * np.ones_like(r_s)
sigma_rr_ana_L = A - B / r_s**2

epsilon_rr_ana_L = epsilon_rr_ref(r_s)
epsilon_tt_ana_L = epsilon_tt_ref(r_s)
epsilon_zz_ana_L = np.ones_like(r_s) * eps_zz_GPS  # generalized plane strain

# Average stress
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
)

plotter_strain_cylinder(
    r_s=r_s,
    strain_rr=epsilon_rr,
    strain_tt=epsilon_tt,
    strain_zz=epsilon_zz,
    Ri=Ri,
    Ro=Ro,
    Pi=Pi,
    Po=Po,
    CASE_DIR=CASE_DIR,
    slenderness=slenderness,
    strain_rr_ana_L=epsilon_rr_ana_L,
    strain_tt_ana_L=epsilon_tt_ana_L,
    strain_zz_ana_L=epsilon_zz_ana_L,
)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_rr = np.sqrt(np.mean((sigma_rr - sigma_rr_ana_L) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_L**2))
err_tt = np.sqrt(np.mean((sigma_tt - sigma_tt_ana_L) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_L**2))
err_zz = np.sqrt(np.mean((sigma_zz - sigma_zz_ana_L) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_L**2))
err_eps_rr = np.sqrt(np.mean((epsilon_rr - epsilon_rr_ana_L) ** 2)) / np.sqrt(
    np.mean(epsilon_rr_ana_L**2)
)
err_eps_tt = np.sqrt(np.mean((epsilon_tt - epsilon_tt_ana_L) ** 2)) / np.sqrt(
    np.mean(epsilon_tt_ana_L**2)
)
err_eps_zz = np.sqrt(np.mean((epsilon_zz - epsilon_zz_ana_L) ** 2)) / np.sqrt(
    np.mean(epsilon_zz_ana_L**2)
)

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
