#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: thick_cylindrical_shell_GPS_2D

non-regression script
---------------------
Analytical non-regression for two cylindrical shells in contact under
internal (Pi) and external (Po) pressure.
Reference is the Lamé solution under generalized plane strain.
"""

import os

import matplotlib
matplotlib.use('Agg')
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
geometry_path = os.path.join(CASE_DIR, "geometry.yaml")
with open(geometry_path, "r") as f:
    geom = yaml.safe_load(f)

Lx_ext = float(geom.get("Lx_ext", 0.20))
Ly_ext = float(geom.get("Ly_ext", 0.30))
Lx_int = float(geom.get("Lx_int", 0.19))
Ly_int = float(geom.get("Ly_int", 0.29))
Lx_mid = float(geom.get("Lx_mid", 0.095))
Ly_mid = float(geom.get("Ly_mid", 0.145))

Ri = Lx_mid / 2.0
Ro = Lx_int / 2.0
t = Lx_ext - Lx_int
slenderness = Ly_ext / t if t != 0 else np.nan

# This case is a 2D plane stress/strain problem in the XY plane.
# The VTU mesh has z = 0 for all points, so the mid-plane should be
# selected in the y direction, not z.
y_target = Ly_ext / 2.0
y_tol = 0.01
Lz = float(geom.get("Lz", Ly_ext))

Pi, Po = 1.0e6, 0.0  # Pa         internal and external pressure
E, nu = 2.0e11, 0.3  # Pa, -      Young modulus, Poisson ratio
P_gap= (Ro-Ri)/((1/E) * ( (Ro**2+Ri**2)/(Ro**2-Ri**2) + nu) + (1-nu)/E)  # (Pa)    
for i in range(1,5) :    
    print("\n")
print(P_gap)
TOLERANCE = 5.0e-3  # -          tolerance for non-regression

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
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

# Numerical results
# Stress
# Stress - outer (steel_o)
x_S_o, y_S_o, z_S_o, S_all_o = extract_field(VTU_FILE, field_name="Stress_steel_o (cells)")
mask_o = np.abs(y_S_o - y_target) < y_tol
sort_idx_o = np.argsort(x_S_o[mask_o])

r_s_o = x_S_o[mask_o][sort_idx_o]

sigma_rr_o = S_all_o[mask_o, 0][sort_idx_o]
sigma_tt_o = S_all_o[mask_o, 4][sort_idx_o]
sigma_zz_o = S_all_o[mask_o, 8][sort_idx_o]

# Stress - mid (steel_mid)
x_S_mid, y_S_mid, z_S_mid, S_all_mid = extract_field(VTU_FILE, field_name="Stress_steel_mid (cells)")
mask_mid = np.abs(y_S_mid - y_target) < y_tol
sort_idx_mid = np.argsort(x_S_mid[mask_mid])

r_s_mid = x_S_mid[mask_mid][sort_idx_mid]

sigma_rr_mid = S_all_mid[mask_mid, 0][sort_idx_mid]
sigma_tt_mid = S_all_mid[mask_mid, 4][sort_idx_mid]
sigma_zz_mid = S_all_mid[mask_mid, 8][sort_idx_mid]

# Strain
x_E, y_E, z_E, E_all = extract_field(VTU_FILE, field_name="Strain (cells)")
mask = np.abs(y_E - y_target) < y_tol
sort_idx = np.argsort(x_E[mask])

epsilon_rr = E_all[mask, 0][sort_idx]
epsilon_tt = E_all[mask, 4][sort_idx]
epsilon_zz = E_all[mask, 8][sort_idx]
# Analytical results for outer
sigma_rr_ana_L_o = A - B / r_s_o**2
sigma_tt_ana_L_o = A + B / r_s_o**2
sigma_zz_ana_L_o = sigma_zz_ana_L * np.ones_like(r_s_o)

epsilon_rr_ana_L_o = epsilon_rr_ref(r_s_o)
epsilon_tt_ana_L_o = epsilon_tt_ref(r_s_o)
epsilon_zz_ana_L_o = np.ones_like(r_s_o) * eps_zz_GPS

# Analytical results for mid
sigma_rr_ana_L_mid = A - B / r_s_mid**2
sigma_tt_ana_L_mid = A + B / r_s_mid**2
sigma_zz_ana_L_mid = sigma_zz_ana_L * np.ones_like(r_s_mid)

epsilon_rr_ana_L_mid = epsilon_rr_ref(r_s_mid)
epsilon_tt_ana_L_mid = epsilon_tt_ref(r_s_mid)
epsilon_zz_ana_L_mid = np.ones_like(r_s_mid) * eps_zz_GPS

# Average stres
sigma_zz_ana_M = sigma_zz_ana_M * np.ones_like(r_s_o)

sigma_rr_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_rr_o * r_s_o, r_s_o)
sigma_tt_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_tt_o * r_s_o, r_s_o)
sigma_zz_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_zz_o * r_s_o, r_s_o)

# Plot
Pa_to_MPa = 1e-6

plt.figure(figsize=(10, 7))
ax1 = plt.gca()
ax1.plot(r_s_o, sigma_rr_o * Pa_to_MPa, "ro", label=r"Num. $\sigma_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, sigma_rr_ana_L_o * Pa_to_MPa, "r-", label=r"Ana. $\sigma_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_s_o, sigma_tt_o * Pa_to_MPa, "go", label=r"Num. $\sigma_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, sigma_tt_ana_L_o * Pa_to_MPa, "g-", label=r"Ana. $\sigma_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_s_o, sigma_zz_o * Pa_to_MPa, "bo", label=r"Num. $\sigma_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, sigma_zz_ana_L_o * Pa_to_MPa, "b-", label=r"Ana. $\sigma_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "stress_comparison_o.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

plt.figure(figsize=(10, 7))
ax1 = plt.gca()
ax1.plot(r_s_o, epsilon_rr, "ro", label=r"Num. $\varepsilon_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, epsilon_rr_ana_L_o, "r-", label=r"Ana. $\varepsilon_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_s_o, epsilon_tt, "go", label=r"Num. $\varepsilon_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, epsilon_tt_ana_L_o, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_s_o, epsilon_zz, "bo", label=r"Num. $\varepsilon_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, epsilon_zz_ana_L_o, "b-", label=r"Ana. $\varepsilon_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Strain (/)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "strain_comparison_o.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- MID RECTANGLE ANALYSIS --.. ..- .-.. .-.. ---
plt.figure(figsize=(10, 7))
ax1 = plt.gca()
ax1.plot(r_s_mid, sigma_rr_mid * Pa_to_MPa, "ro", label=r"Num. $\sigma_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s_mid, sigma_rr_ana_L_mid * Pa_to_MPa, "r-", label=r"Ana. $\sigma_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_s_mid, sigma_tt_mid * Pa_to_MPa, "go", label=r"Num. $\sigma_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_s_mid, sigma_tt_ana_L_mid * Pa_to_MPa, "g-", label=r"Ana. $\sigma_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_s_mid, sigma_zz_mid * Pa_to_MPa, "bo", label=r"Num. $\sigma_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s_mid, sigma_zz_ana_L_mid * Pa_to_MPa, "b-", label=r"Ana. $\sigma_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "stress_comparison_mid.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

plt.figure(figsize=(10, 7))
ax1 = plt.gca()
ax1.plot(r_s_mid, epsilon_rr, "ro", label=r"Num. $\varepsilon_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s_mid, epsilon_rr_ana_L_mid, "r-", label=r"Ana. $\varepsilon_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_s_mid, epsilon_tt, "go", label=r"Num. $\varepsilon_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_s_mid, epsilon_tt_ana_L_mid, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_s_mid, epsilon_zz, "bo", label=r"Num. $\varepsilon_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s_mid, epsilon_zz_ana_L_mid, "b-", label=r"Ana. $\varepsilon_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Strain (/)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "strain_comparison_mid.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_rr = np.sqrt(np.mean((sigma_rr_o - sigma_rr_ana_L_o) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_L_o**2))
err_tt = np.sqrt(np.mean((sigma_tt_o - sigma_tt_ana_L_o) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_L_o**2))
err_zz = np.sqrt(np.mean((sigma_zz_o - sigma_zz_ana_L_o) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_L_o**2))
err_eps_rr = np.sqrt(np.mean((epsilon_rr - epsilon_rr_ana_L_o) ** 2)) / np.sqrt(np.mean(epsilon_rr_ana_L_o**2))
err_eps_tt = np.sqrt(np.mean((epsilon_tt - epsilon_tt_ana_L_o) ** 2)) / np.sqrt(np.mean(epsilon_tt_ana_L_o**2))
err_eps_zz = np.sqrt(np.mean((epsilon_zz - epsilon_zz_ana_L_o) ** 2)) / np.sqrt(np.mean(epsilon_zz_ana_L_o**2))

# --.. ..- .-.. .-.. --- metrics for MID rectangle --.. ..- .-.. .-.. ---
err_rr_mid = np.sqrt(np.mean((sigma_rr_mid - sigma_rr_ana_L_mid) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_L_mid**2))
err_tt_mid = np.sqrt(np.mean((sigma_tt_mid - sigma_tt_ana_L_mid) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_L_mid**2))
err_zz_mid = np.sqrt(np.mean((sigma_zz_mid - sigma_zz_ana_L_mid) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_L_mid**2))
err_eps_rr_mid = np.sqrt(np.mean((epsilon_rr - epsilon_rr_ana_L_mid) ** 2)) / np.sqrt(np.mean(epsilon_rr_ana_L_mid**2))
err_eps_tt_mid = np.sqrt(np.mean((epsilon_tt - epsilon_tt_ana_L_mid) ** 2)) / np.sqrt(np.mean(epsilon_tt_ana_L_mid**2))
err_eps_zz_mid = np.sqrt(np.mean((epsilon_zz - epsilon_zz_ana_L_mid) ** 2)) / np.sqrt(np.mean(epsilon_zz_ana_L_mid**2))

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
    # MID rectangle errors
    "L2_error_sigma_rr_mid": {
        "numerical": float(err_rr_mid),
        "reference": 0.0,
        "abs_error": float(err_rr_mid),
        "rel_error": float(err_rr_mid),
    },
    "L2_error_sigma_tt_mid": {
        "numerical": float(err_tt_mid),
        "reference": 0.0,
        "abs_error": float(err_tt_mid),
        "rel_error": float(err_tt_mid),
    },
    "L2_error_sigma_zz_mid": {
        "numerical": float(err_zz_mid),
        "reference": 0.0,
        "abs_error": float(err_zz_mid),
        "rel_error": float(err_zz_mid),
    },
    "L2_error_strain_rr_mid": {
        "numerical": float(err_eps_rr_mid),
        "reference": 0.0,
        "abs_error": float(err_eps_rr_mid),
        "rel_error": float(err_eps_rr_mid),
    },
    "L2_error_strain_tt_mid": {
        "numerical": float(err_eps_tt_mid),
        "reference": 0.0,
        "abs_error": float(err_eps_tt_mid),
        "rel_error": float(err_eps_tt_mid),
    },
    "L2_error_strain_zz_mid": {
        "numerical": float(err_eps_zz_mid),
        "reference": 0.0,
        "abs_error": float(err_eps_zz_mid),
        "rel_error": float(err_eps_zz_mid),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO]Non-regression completed.\n")