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

R_o = float(geom.get("R_o", 0.005))
Ly_rect = float(geom.get("Ly_rect", 0.01))
R_int = float(geom.get("R_int", 0.001))
R_mid = float(geom.get("R_mid", 0.003))
gap=0.001
delta=gap

# Les Lx sont maintenant des rayons directs (coordonnées x du maillage)
Ri = R_o - R_int  # Rayon interne de la gaine (0.004)
Ro = R_o          # Rayon externe de la gaine (0.005)
t = Ro - Ri
slenderness = Ly_rect / t if t != 0 else np.nan

# This case is a 2D plane stress/strain problem in the XY plane.
# The VTU mesh has z = 0 for all points, so the in-plane should be
# selected in the y direction, not z.
y_target = Ly_rect / 2.0
y_tol = 0.01
Lz = float(geom.get("Lz", Ly_rect))

Pi, Po = 0.0, 0.0  # Pa         internal and external pressure
E_in, nu_in = 2.0e11, 0.3  # Pa, -      Young modulus, Poisson ratio
E_o, nu_o = 2.0e11, 0.3  # Pa, -      Young modulus, Poisson ratio

P_gap=(gap+delta_exp) / (R_mid * ( (1/E_in) * (nu_in-(R_mid**2 + R_int**2) / (R_mid**2 - R_int**2)) + (1 / E_o) * ( nu_o+ (R_o**2 + R_mid**2) / (-R_mid**2 + R_o**2 ) )))   # Pa on the interface (inward)    


for i in range(1,5) :    
    print("\n")
print(P_gap)
TOLERANCE = 5.0e-3  # -          tolerance for non-regression

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
eps_zz_GPS = -2.718310e-06

# Lamé solutions
A_in = (Pi * R_int**2 - P_gap * R_mid**2) / (R_mid**2 - R_int**2)
B_in = -(R_int**2 * R_mid**2 * (Pi - P_gap)) / (R_mid**2 - R_int**2)

A_o = (P_gap * Ri**2 - Po * Ro**2) / (Ro**2 - Ri**2)
B_o = (Ri**2 * Ro**2 * (P_gap - Po)) / (Ro**2 - Ri**2)
sigma_zz_ana_L_in = 2 * nu_in * A_in + E_in * eps_zz_GPS  # Generalized plane strain (epsilon_z = const)
sigma_zz_ana_L_o = 2 * nu_o * A_o + E_o * eps_zz_GPS  # Generalized plane strain (epsilon_z = const)


def epsilon_rr_ref(r,E,nu,A,B):
    return (1 + nu) / E * (A * (1 - 2 * nu) - B / r**2) - nu * eps_zz_GPS


def epsilon_tt_ref(r,nu,E,A,B):
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
# Extraction of the Stress
x_S, y_S, z_S, S_all = extract_field(VTU_FILE, field_name="Stress_steel_o (cells)")

outer_mask_x = (x_S >= R_mid+gap/2) & (x_S <= R_o)
inner_mask_x = (x_S >= R_int) & (x_S <= R_mid-gap/2)

x_S_o, y_S_o, z_S_o, S_all_o = x_S[outer_mask_x], y_S[outer_mask_x], z_S[outer_mask_x], S_all[outer_mask_x]
x_S_in, y_S_in, z_S_in, S_all_in = x_S[inner_mask_x], y_S[inner_mask_x], z_S[inner_mask_x], S_all[inner_mask_x]

mask_o = np.abs(y_S_o - y_target) < y_tol
mask_in = np.abs(y_S_in - y_target) < y_tol

sort_idx_in = np.argsort(x_S_in[mask_in])
sort_idx_o = np.argsort(x_S_o[mask_o])

r_s_o = x_S_o[mask_o][sort_idx_o]
r_s_in = x_S_in[mask_in][sort_idx_in]


""" for i in range(1,5) :
    print("\n")

print(len(S_all_o), "taille de y_S_o")
print(len(x_S_o)," taille de x_S_o")
print(len(r_s_o), " taille de rs_o")
print(len(r_s_in), " taille de rs_in")
 """
sigma_rr_o = S_all_o[mask_o, 0][sort_idx_o]
sigma_tt_o = S_all_o[mask_o, 4][sort_idx_o]
sigma_zz_o = S_all_o[mask_o, 8][sort_idx_o]


sigma_rr_in = S_all_in[mask_in, 0][sort_idx_in]
sigma_tt_in = S_all_in[mask_in, 4][sort_idx_in]
sigma_zz_in = S_all_in[mask_in, 8][sort_idx_in]


""" print(len(sigma_rr_in), "taille de sigma_rr_in")
print(sigma_rr_in)
 """
# Extraction of the Strain
x_E, y_E, z_E, E_all = extract_field(VTU_FILE, field_name="Strain (cells)")
outer_mask_x_E = (x_E >= R_mid+gap/2) & (x_E <= R_o)
inner_mask_x_E = (x_E >= R_int) & (x_E <= R_mid-gap/2)

x_E_o, y_E_o, z_E_o, E_all_o = x_E[outer_mask_x_E], y_E[outer_mask_x_E], z_E[outer_mask_x_E], E_all[outer_mask_x_E]
x_E_in, y_E_in, z_E_in, E_all_in = x_E[inner_mask_x_E], y_E[inner_mask_x_E], z_E[inner_mask_x_E], E_all[inner_mask_x_E]

mask_e_o = np.abs(y_E_o - y_target) < y_tol
mask_e_in = np.abs(y_E_in - y_target) < y_tol

sort_idx_e_o = np.argsort(x_E_o[mask_e_o])
sort_idx_e_in = np.argsort(x_E_in[mask_e_in])

r_e_o = x_E_o[mask_e_o][sort_idx_e_o]
r_e_in = x_E_in[mask_e_in][sort_idx_e_in]

epsilon_rr_o = E_all_o[mask_e_o, 0][sort_idx_e_o]
epsilon_tt_o = E_all_o[mask_e_o, 4][sort_idx_e_o]
epsilon_zz_o = E_all_o[mask_e_o, 8][sort_idx_e_o]

epsilon_rr_in = E_all_in[mask_e_in, 0][sort_idx_e_in]
epsilon_tt_in = E_all_in[mask_e_in, 4][sort_idx_e_in]
epsilon_zz_in = E_all_in[mask_e_in, 8][sort_idx_e_in]

# Analytical results for the outer cylinder


sigma_rr_ana_L_o = A_o - B_o / r_s_o**2
sigma_tt_ana_L_o = A_o + B_o / r_s_o**2
sigma_zz_ana_L_o = sigma_zz_ana_L_o * np.ones_like(r_s_o)

print(len(sigma_rr_ana_L_o), " taille de sigma_rr_ana_L_o")


epsilon_rr_ana_L_o = epsilon_rr_ref(r_e_o, E_o, nu_o, A_o, B_o)
epsilon_tt_ana_L_o = epsilon_tt_ref(r_e_o, nu_o, E_o, A_o, B_o)
epsilon_zz_ana_L_o = np.ones_like(r_e_o) * eps_zz_GPS

# Analytical results for the inner cylinder

sigma_rr_ana_L_in = A_in - B_in / r_s_in**2 
sigma_tt_ana_L_in = A_in + B_in / r_s_in**2
sigma_zz_ana_L_in = sigma_zz_ana_L_in * np.ones_like(r_s_in)

epsilon_rr_ana_L_in = epsilon_rr_ref(r_e_in, E_in, nu_in, A_in, B_in)
epsilon_tt_ana_L_in = epsilon_tt_ref(r_e_in, nu_in, E_in, A_in, B_in)
epsilon_zz_ana_L_in = np.ones_like(r_e_in) * eps_zz_GPS

# Average stress
sigma_zz_ana_M = sigma_zz_ana_M * np.ones_like(r_s_o)

sigma_rr_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_rr_o * r_s_o, r_s_o)
sigma_tt_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_tt_o * r_s_o, r_s_o)
sigma_zz_avg = 2 / (Ro**2 - Ri**2) * np.trapezoid(sigma_zz_o * r_s_o, r_s_o)

# Plot
Pa_to_MPa = 1e-6

print(len(r_s_o), "taille de r_s_o")

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
ax1.plot(r_s_o, epsilon_rr_o, "ro", label=r"Num. $\varepsilon_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, epsilon_rr_ana_L_o, "r-", label=r"Ana. $\varepsilon_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_s_o, epsilon_tt_o, "go", label=r"Num. $\varepsilon_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, epsilon_tt_ana_L_o, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_s_o, epsilon_zz_o, "bo", label=r"Num. $\varepsilon_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s_o, epsilon_zz_ana_L_o, "b-", label=r"Ana. $\varepsilon_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Strain (/)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "strain_comparison_o.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- in RECTANGLE ANALYSIS --.. ..- .-.. .-.. ---
plt.figure(figsize=(10, 7))
ax1 = plt.gca()
ax1.plot(r_s_in, sigma_rr_in * Pa_to_MPa, "ro", label=r"Num. $\sigma_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_s_in, sigma_rr_ana_L_in * Pa_to_MPa, "r-", label=r"Ana. $\sigma_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_s_in, sigma_tt_in * Pa_to_MPa, "go", label=r"Num. $\sigma_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_s_in, sigma_tt_ana_L_in * Pa_to_MPa, "g-", label=r"Ana. $\sigma_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_s_in, sigma_zz_in * Pa_to_MPa, "bo", label=r"Num. $\sigma_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_s_in, sigma_zz_ana_L_in * Pa_to_MPa, "b-", label=r"Ana. $\sigma_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "stress_comparison_in.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

plt.figure(figsize=(10, 7))
ax1 = plt.gca()
ax1.plot(r_e_in, epsilon_rr_in, "ro", label=r"Num. $\varepsilon_{rr}$ (Radial)", markersize=4, alpha=0.6)
ax1.plot(r_e_in, epsilon_rr_ana_L_in, "r-", label=r"Ana. $\varepsilon_{rr}$ (Radial)", linewidth=1.5)
ax1.plot(r_e_in, epsilon_tt_in, "go", label=r"Num. $\varepsilon_{\theta\theta}$ (Hoop)", markersize=4, alpha=0.6)
ax1.plot(r_e_in, epsilon_tt_ana_L_in, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (Hoop)", linewidth=1.5)
ax1.plot(r_e_in, epsilon_zz_in, "bo", label=r"Num. $\varepsilon_{zz}$ (Axial)", markersize=4, alpha=0.6)
ax1.plot(r_e_in, epsilon_zz_ana_L_in, "b-", label=r"Ana. $\varepsilon_{zz}$ (Axial)", linewidth=1.5)
ax1.set_xlabel("Radius (m)", fontsize=12)
ax1.set_ylabel("Strain (/)", fontsize=12)
ax1.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "strain_comparison_in.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_rr = np.sqrt(np.mean((sigma_rr_o - sigma_rr_ana_L_o) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_L_o**2))
err_tt = np.sqrt(np.mean((sigma_tt_o - sigma_tt_ana_L_o) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_L_o**2))
err_zz = np.sqrt(np.mean((sigma_zz_o - sigma_zz_ana_L_o) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_L_o**2))
err_eps_rr = np.sqrt(np.mean((epsilon_rr_o - epsilon_rr_ana_L_o) ** 2)) / np.sqrt(np.mean(epsilon_rr_ana_L_o**2))
err_eps_tt = np.sqrt(np.mean((epsilon_tt_o - epsilon_tt_ana_L_o) ** 2)) / np.sqrt(np.mean(epsilon_tt_ana_L_o**2))
err_eps_zz = np.sqrt(np.mean((epsilon_zz_o - epsilon_zz_ana_L_o) ** 2)) / np.sqrt(np.mean(epsilon_zz_ana_L_o**2))

# --.. ..- .-.. .-.. --- metrics for in rectangle --.. ..- .-.. .-.. ---
err_rr_in = np.sqrt(np.mean((sigma_rr_in - sigma_rr_ana_L_in) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_L_in**2))
err_tt_in = np.sqrt(np.mean((sigma_tt_in - sigma_tt_ana_L_in) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_L_in**2))
err_zz_in = np.sqrt(np.mean((sigma_zz_in - sigma_zz_ana_L_in) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_L_in**2))
err_eps_rr_in = np.sqrt(np.mean((epsilon_rr_in - epsilon_rr_ana_L_in) ** 2)) / np.sqrt(np.mean(epsilon_rr_ana_L_in**2))
err_eps_tt_in = np.sqrt(np.mean((epsilon_tt_in - epsilon_tt_ana_L_in) ** 2)) / np.sqrt(np.mean(epsilon_tt_ana_L_in**2))
err_eps_zz_in = np.sqrt(np.mean((epsilon_zz_in - epsilon_zz_ana_L_in) ** 2)) / np.sqrt(np.mean(epsilon_zz_ana_L_in**2))

plt.plot(r_s_o,sigma_rr_o)
plt.show()


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
    # in rectangle errors
    "L2_error_sigma_rr_in": {
        "numerical": float(err_rr_in),
        "reference": 0.0,
        "abs_error": float(err_rr_in),
        "rel_error": float(err_rr_in),
    },
    "L2_error_sigma_tt_in": {
        "numerical": float(err_tt_in),
        "reference": 0.0,
        "abs_error": float(err_tt_in),
        "rel_error": float(err_tt_in),
    },
    "L2_error_sigma_zz_in": {
        "numerical": float(err_zz_in),
        "reference": 0.0,
        "abs_error": float(err_zz_in),
        "rel_error": float(err_zz_in),
    },
    "L2_error_strain_rr_in": {
        "numerical": float(err_eps_rr_in),
        "reference": 0.0,
        "abs_error": float(err_eps_rr_in),
        "rel_error": float(err_eps_rr_in),
    },
    "L2_error_strain_tt_in": {
        "numerical": float(err_eps_tt_in),
        "reference": 0.0,
        "abs_error": float(err_eps_tt_in),
        "rel_error": float(err_eps_tt_in),
    },
    "L2_error_strain_zz_in": {
        "numerical": float(err_eps_zz_in),
        "reference": 0.0,
        "abs_error": float(err_eps_zz_in),
        "rel_error": float(err_eps_zz_in),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO]Non-regression completed.\n")

delta_exp=epsilon_tt_o[0]*R_mid




print( abs(delta_exp**2-delta**2)/delta**2, "erreur en r=R_mid")
print((delta_exp), "delta_exp")
print(len(S_all_o),len(S_all_o[0]))
print(epsilon_tt_o[0], "epsilon_tt at r=R_mid")
print(x_S_o[-1], "x_S_0 last point")
print(x_S_in[-1], "x_S_in last point")

# Filtrer les coordonnées x (rayons) par la tolérance en y, puis les trier
mask_for_x_S_print = np.abs(y_S - y_target) < y_tol
x_S_filtered_sorted = np.sort(x_S[mask_for_x_S_print])
with np.printoptions(threshold=np.inf):
    print(x_S_filtered_sorted)
