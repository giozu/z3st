#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: thin_cylindrical_shell_Mariotte

2D thin-walled cylindrical shell under internal (Pi) pressure. Reference is the
Mariotte solution.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import *
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
E, nu = float(mat["E"]), float(mat["nu"])  # Pa, -
t = Ro - Ri  # m          wall thickness
slenderness = Ri / t  # -          slenderness ratio
z_target, z_tol = Lz, 1.0  # m          z-plane for data extraction

TOLERANCE = 6.0e-1  # -          tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
# Mariotte solutions (Po=0)
sigma_rr_ana_M = -Pi / 2
sigma_tt_ana_M = Pi * Ri / t
sigma_zz_ana_M = Pi * Ri / (2 * t)

# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
# Numerical results
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

# Numerical results
# Stress
x_S, z_S, _, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")
mask = np.abs(z_S - z_target) < z_tol
sort_idx = np.argsort(x_S[mask])

r_s = x_S[mask][sort_idx]

sigma_rr = S_all[mask, 0][sort_idx]
sigma_tt = S_all[mask, 4][sort_idx]
sigma_zz = S_all[mask, 8][sort_idx]

# Strain
x_E, z_E, _, E_all = extract_field(VTU_FILE, field_name="Strain (cells)")
mask = np.abs(z_E - z_target) < z_tol
sort_idx = np.argsort(x_E[mask])

epsilon_rr = E_all[mask, 0][sort_idx]
epsilon_tt = E_all[mask, 4][sort_idx]
epsilon_zz = E_all[mask, 8][sort_idx]

# Analytical results
# Mariotte
sigma_rr_ana_M = sigma_rr_ana_M * np.ones_like(r_s)
sigma_tt_ana_M = sigma_tt_ana_M * np.ones_like(r_s)
sigma_zz_ana_M = sigma_zz_ana_M * np.ones_like(r_s)

epsilon_rr_ana_M = 1 / E * (sigma_rr_ana_M - nu * (sigma_tt_ana_M + sigma_zz_ana_M))
epsilon_tt_ana_M = 1 / E * (sigma_tt_ana_M - nu * (sigma_rr_ana_M + sigma_zz_ana_M))
epsilon_zz_ana_M = 1 / E * (sigma_zz_ana_M - nu * (sigma_rr_ana_M + sigma_tt_ana_M))

# average stress
sigma_rr_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_rr * r_s, r_s)
sigma_tt_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_tt * r_s, r_s)
sigma_zz_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_zz * r_s, r_s)


# Plot
Pa_to_MPa = 1e-6

plt.figure(figsize=(10, 7))

# Stress
ax1 = plt.gca()
ax1.plot(
    r_s, sigma_rr * Pa_to_MPa, "ro", label=r"Num. $\sigma_{rr}$ (Radial)", markersize=4, alpha=0.6
)
ax1.plot(r_s, sigma_rr_ana_M * Pa_to_MPa, "r-", label=r"Ana. $\sigma_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(
    r_s,
    sigma_tt * Pa_to_MPa,
    "go",
    label=r"Num. $\sigma_{\theta\theta}$ (Hoop)",
    markersize=4,
    alpha=0.6,
)
ax1.plot(
    r_s,
    sigma_tt_ana_M * Pa_to_MPa,
    "g-",
    label=r"Ana. $\sigma_{\theta\theta}$ (Hoop)",
    linewidth=1.5,
)
ax1.plot(
    r_s, sigma_zz * Pa_to_MPa, "bo", label=r"Num. $\sigma_{zz}$ (Axial)", markersize=4, alpha=0.6
)
ax1.plot(r_s, sigma_zz_ana_M * Pa_to_MPa, "b-", label=r"Ana. $\sigma_{zz}$ (Axial)", linewidth=1.5)

ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)

plt.legend()
plt.tight_layout()

plot_path = os.path.join(CASE_DIR, "output", "stress_comparison.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")


# Plot
plt.figure(figsize=(10, 7))

# Strain
ax1 = plt.gca()
ax1.plot(r_s, epsilon_rr, "ro", label=r"Num. $\varepsilon_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s, epsilon_rr_ana_M, "r-", label=r"Ana. $\varepsilon_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(
    r_s,
    epsilon_tt,
    "go",
    label=r"Num. $\varepsilon_{\theta\theta}$ (Hoop)",
    markersize=4,
    alpha=0.6,
)
ax1.plot(
    r_s, epsilon_tt_ana_M, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (Hoop)", linewidth=1.5
)
ax1.plot(r_s, epsilon_zz, "bo", label=r"Num. $\varepsilon_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s, epsilon_zz_ana_M, "b-", label=r"Ana. $\varepsilon_{zz}$ (Axial)", linewidth=1.5)

ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Strain (/)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)

plt.legend()
plt.tight_layout()

plot_path = os.path.join(CASE_DIR, "output", "strain_comparison.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")


# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_rr = np.sqrt(np.mean((sigma_rr - sigma_rr_ana_M) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_M**2))
err_tt = np.sqrt(np.mean((sigma_tt - sigma_tt_ana_M) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_M**2))
err_zz = np.sqrt(np.mean((sigma_zz - sigma_zz_ana_M) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_M**2))
err_eps_rr = np.sqrt(np.mean((epsilon_rr - epsilon_rr_ana_M) ** 2)) / np.sqrt(
    np.mean(epsilon_rr_ana_M**2)
)
err_eps_tt = np.sqrt(np.mean((epsilon_tt - epsilon_tt_ana_M) ** 2)) / np.sqrt(
    np.mean(epsilon_tt_ana_M**2)
)
err_eps_zz = np.sqrt(np.mean((epsilon_zz - epsilon_zz_ana_M) ** 2)) / np.sqrt(
    np.mean(epsilon_zz_ana_M**2)
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
    "L2_error_epsilon_rr": {
        "numerical": float(err_eps_rr),
        "reference": 0.0,
        "abs_error": float(err_eps_rr),
        "rel_error": float(err_eps_rr),
    },
    "L2_error_epsilon_tt": {
        "numerical": float(err_eps_tt),
        "reference": 0.0,
        "abs_error": float(err_eps_tt),
        "rel_error": float(err_eps_tt),
    },
    "L2_error_epsilon_zz": {
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
