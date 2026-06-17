#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: CONTACT_2D_GLOB/After_contact

non-regression script
---------------------
Analytical/verification script for the after_contact case.
This case models two concentric axisymmetric shells in contact,
with a pressure applied at the interface and a pressure on the outer face.
"""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.utils_extract_vtu import extract_field, list_fields
from z3st.utils.utils_verification import pass_fail_check, regression_check

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Load geometry
geometry_path = os.path.join(CASE_DIR, "geometry.yaml")
with open(geometry_path, "r") as f:
    geom = yaml.safe_load(f)

R_int = float(geom.get("R_int", 0.001))
R_mid = float(geom.get("R_mid", 0.004))
R_o = float(geom.get("R_o", 0.005))
Ly = float(geom.get("Ly", 0.30))
gap=1.0e-4 # gap for strain extraction (to avoid singularity at the interface)



# Material constants (use the same values as the case materials)
E_in = 2.0e11  # Pa
nu_in = 0.3
E_o = 358.0e9   # Pa
nu_o = 0.3

# Pressures
""" # Calcul de la pression de contact analytique (Pellet plein dans Gaine)
term_cladding = (1 / E_o) * ((R_o**2 + R_mid**2) / (R_o**2 - R_mid**2) + nu_o)
term_pellet = (1 - nu_in) / E_in
P_interface = delta / (R_mid * (term_cladding + term_pellet)) """

P_interface = gap / (R_mid * ( (1/E_in) * (R_mid**2 + R_int**2*(1+nu_in)) / (R_mid**2 - R_int**2) + (1 / (E_o * (R_o**2 - R_mid**2))) * (R_mid**2 + R_o**2 - nu_o*(2*R_mid**2 - R_o**2)) ))  # Pa on the interface (inward)    

P_outer = 1.0e6      # Pa on the outer radius (inward)
for i in range(1,5):
    print("\n")
print(P_interface)
TOLERANCE = 5.0e-3  # -

def epsilon_rr_ref(r, A, B, nu, E):
    return (1 + nu) / E * (A * (1 - 2 * nu) - B / r**2)


def epsilon_tt_ref(r, A, B, nu, E):
    return (1 + nu) / E * (A * (1 - 2 * nu) + B / r**2)


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# Mid-plane extraction of the axisymmetric model
y_target = Ly / 2.0
# NOTE: After_contact uses ny=500 divisions in y (vs 81 in case 9)

y_tol = 0.06  # Scaled to match case 9 behavior (1-2 points per radius)

# Helper that returns a stress field name candidate list
stress_field_names = [
    "Stress_steel_o (cells)",
    "Stress_steel_o",
    "Stress_steel_in (cells)",
    "Stress_steel_in",
]

# Extract fields from VTU (single stress field, then split by masks)
# Concatenate inner and outer stress fields to allow region-based masking (in, o, gap)
x_So, y_So, z_So, S_o = extract_field(VTU_FILE, field_name="Stress_steel_o (cells)")
x_Si, y_Si, z_Si, S_i = extract_field(VTU_FILE, field_name="Stress_steel_in (cells)")
x_Sgap, y_Sgap, z_Sgap, S_gap = extract_field(VTU_FILE, field_name="Stress_gap_region (cells)")

# Separate cells into outer and inner regions using radial ranges
outer_mask_x = (x_So >= R_mid) & (x_So <= R_o)
inner_mask_x = (x_Si >= R_int) & (x_Si <= R_mid-gap/2)
gap_mask_x = (x_Sgap > R_mid - gap/2) & (x_Sgap < R_mid + gap/2)

# Per-region coordinate and field arrays
x_S_o, y_S_o, z_S_o, S_all_o = x_So[outer_mask_x], y_So[outer_mask_x], z_So[outer_mask_x], S_o[outer_mask_x]
x_S_in, y_S_in, z_S_in, S_all_in = x_Si[inner_mask_x], y_Si[inner_mask_x], z_Si[inner_mask_x], S_i[inner_mask_x]
x_S_gap, y_S_gap, z_S_gap, S_all_gap = x_Sgap[gap_mask_x], y_Sgap[gap_mask_x], z_Sgap[gap_mask_x], S_gap[gap_mask_x]

# Select the mid-plane
mask_o = np.abs(y_S_o - y_target) < y_tol
mask_in = np.abs(y_S_in - y_target) < y_tol
mask_gap = np.abs(y_S_gap - y_target) < y_tol

# Debug: check for duplicate x values before sorting
""" x_S_o_filtered = x_S_o[mask_o]
x_unique_o, counts_o = np.unique(x_S_o_filtered, return_counts=True)
duplicates_o = x_unique_o[counts_o > 1] """

""" print(f"[DEBUG] Outer shell: total filtered points = {len(x_S_o_filtered)}")
print(f"[DEBUG] Outer shell: unique x values = {len(x_unique_o)}")
if len(duplicates_o) > 0:
    print(f"[DEBUG] WARNING: Found {len(duplicates_o)} duplicate x values in outer shell: {duplicates_o}")
    for x_val in duplicates_o[:3]:  # Affiche les 3 premiers doublons
        dup_mask = np.abs(x_S_o_filtered - x_val) < 1e-10
        print(f"  x = {x_val}: {np.sum(dup_mask)} points") """

sort_idx_o = np.argsort(x_S_o[mask_o])
sort_idx_in = np.argsort(x_S_in[mask_in])

# Utilise les rayons triés directement, sans agrégation
r_o = x_S_o[mask_o][sort_idx_o]
r_in = x_S_in[mask_in][sort_idx_in]


# Utilise les valeurs de contrainte triées directement, sans agrégation
S_all_o_sorted = S_all_o[mask_o][sort_idx_o]
S_all_in_sorted = S_all_in[mask_in][sort_idx_in]

sigma_rr_o = S_all_o_sorted[:, 0]
sigma_tt_o = S_all_o_sorted[:, 4]
sigma_zz_o = S_all_o_sorted[:, 8]

sigma_rr_in = S_all_in_sorted[:, 0]
sigma_tt_in = S_all_in_sorted[:, 4]
sigma_zz_in = S_all_in_sorted[:, 8]

# Extraction des déformations
x_E, y_E, z_E, E_all = extract_field(VTU_FILE, field_name="Strain (cells)")

outer_mask_x_E = (x_E >= R_mid) & (x_E <= R_o)
inner_mask_x_E = (x_E >= R_int) & (x_E <= R_mid-gap/2)
gap_mask_x_E = (x_E > R_mid - gap/2) & (x_E < R_mid + gap/2)

x_E_o, y_E_o, z_E_o, E_all_o = x_E[outer_mask_x_E], y_E[outer_mask_x_E], z_E[outer_mask_x_E], E_all[outer_mask_x_E]
x_E_in, y_E_in, z_E_in, E_all_in = x_E[inner_mask_x_E], y_E[inner_mask_x_E], z_E[inner_mask_x_E], E_all[inner_mask_x_E]
x_E_gap, y_E_gap, z_E_gap, E_all_gap = x_E[gap_mask_x_E], y_E[gap_mask_x_E], z_E[gap_mask_x_E], E_all[gap_mask_x_E]

mask_e_o = np.abs(y_E_o - y_target) < y_tol
mask_e_in = np.abs(y_E_in - y_target) < y_tol
mask_e_gap = np.abs(y_E_gap - y_target) < y_tol

sort_idx_e_o = np.argsort(x_E_o[mask_e_o])
sort_idx_e_in = np.argsort(x_E_in[mask_e_in])

# Utilise les rayons triés directement, sans agrégation
r_e_o = x_E_o[mask_e_o][sort_idx_e_o]
r_e_in = x_E_in[mask_e_in][sort_idx_e_in]



# Utilise les valeurs de déformation triées directement, sans agrégation
E_all_o_sorted = E_all_o[mask_e_o][sort_idx_e_o]
E_all_in_sorted = E_all_in[mask_e_in][sort_idx_e_in]

epsilon_rr_o = E_all_o_sorted[:, 0]
epsilon_tt_o = E_all_o_sorted[:, 4]
epsilon_zz_o = E_all_o_sorted[:, 8]

epsilon_rr_in = E_all_in_sorted[:, 0]
epsilon_tt_in = E_all_in_sorted[:, 4]
epsilon_zz_in = E_all_in_sorted[:, 8]

# Analytical values for the outer shell
A_o = (P_interface * R_mid**2 - P_outer * R_o**2) / (R_o**2 - R_mid**2)
B_o = (R_mid**2 * R_o**2 * (P_interface - P_outer)) / (R_o**2 - R_mid**2)

sigma_rr_ana_o = A_o - B_o / r_o**2
sigma_tt_ana_o = A_o + B_o / r_o**2
sigma_zz_ana_o = nu_o * (sigma_rr_ana_o + sigma_tt_ana_o)

epsilon_rr_ana_o = epsilon_rr_ref(r_e_o, A_o, B_o, nu_o, E_o)
epsilon_tt_ana_o = epsilon_tt_ref(r_e_o, A_o, B_o, nu_o, E_o)
epsilon_zz_ana_o = np.zeros_like(r_e_o)

# Analytical values for the inner shell
A_in = (0.0 * R_int**2 - P_interface * R_mid**2) / (R_mid**2 - R_int**2)
B_in = (R_int**2 * R_mid**2 * (0.0 - P_interface)) / (R_mid**2 - R_int**2)

sigma_rr_ana_in = A_in - B_in / r_in**2
sigma_tt_ana_in = A_in + B_in / r_in**2
sigma_zz_ana_in = nu_in * (sigma_rr_ana_in + sigma_tt_ana_in)

epsilon_rr_ana_in = epsilon_rr_ref(r_in, A_in, B_in, nu_in, E_in)
epsilon_tt_ana_in = epsilon_tt_ref(r_in, A_in, B_in, nu_in, E_in)
epsilon_zz_ana_in = np.zeros_like(r_in)

# Plot and save
os.makedirs(os.path.join(CASE_DIR, "output"), exist_ok=True)
Pa_to_MPa = 1e-6

plt.figure(figsize=(10, 7))
plt.plot(r_o, sigma_rr_o * Pa_to_MPa, "ro", label=r"Num. $\sigma_{rr}$ (outer)", markersize=4, alpha=0.6)
plt.plot(r_o, sigma_rr_ana_o * Pa_to_MPa, "r-", label=r"Ana. $\sigma_{rr}$ (outer)", linewidth=1.5)
plt.plot(r_o, sigma_tt_o * Pa_to_MPa, "go", label=r"Num. $\sigma_{\theta\theta}$ (outer)", markersize=4, alpha=0.6)
plt.plot(r_o, sigma_tt_ana_o * Pa_to_MPa, "g-", label=r"Ana. $\sigma_{\theta\theta}$ (outer)", linewidth=1.5)
plt.plot(r_o, sigma_zz_o * Pa_to_MPa, "bo", label=r"Num. $\sigma_{zz}$ (outer)", markersize=4, alpha=0.6)
plt.plot(r_o, sigma_zz_ana_o * Pa_to_MPa, "b-", label=r"Ana. $\sigma_{zz}$ (outer)", linewidth=1.5)
plt.xlabel("Radius (m)", fontsize=12)
plt.ylabel("Stress (MPa)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "stress_comparison_outer.png"), dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
plt.plot(r_e_o, epsilon_rr_o, "ro", label=r"Num. $\varepsilon_{rr}$ (outer)", markersize=4, alpha=0.6)
plt.plot(r_e_o, epsilon_rr_ana_o, "r-", label=r"Ana. $\varepsilon_{rr}$ (outer)", linewidth=1.5)
plt.plot(r_e_o, epsilon_tt_o, "go", label=r"Num. $\varepsilon_{\theta\theta}$ (outer)", markersize=4, alpha=0.6)
plt.plot(r_e_o, epsilon_tt_ana_o, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (outer)", linewidth=1.5)
plt.plot(r_e_o, epsilon_zz_o, "bo", label=r"Num. $\varepsilon_{zz}$ (outer)", markersize=4, alpha=0.6)
plt.plot(r_e_o, epsilon_zz_ana_o, "b-", label=r"Ana. $\varepsilon_{zz}$ (outer)", linewidth=1.5)
plt.xlabel("Radius (m)", fontsize=12)
plt.ylabel("Strain (/)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "strain_comparison_outer.png"), dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
plt.plot(r_in, sigma_rr_in * Pa_to_MPa, "ro", label=r"Num. $\sigma_{rr}$ (inner)", markersize=4, alpha=0.6)
plt.plot(r_in, sigma_rr_ana_in * Pa_to_MPa, "r-", label=r"Ana. $\sigma_{rr}$ (inner)", linewidth=1.5)
plt.plot(r_in, sigma_tt_in * Pa_to_MPa, "go", label=r"Num. $\sigma_{\theta\theta}$ (inner)", markersize=4, alpha=0.6)
plt.plot(r_in, sigma_tt_ana_in * Pa_to_MPa, "g-", label=r"Ana. $\sigma_{\theta\theta}$ (inner)", linewidth=1.5)
plt.plot(r_in, sigma_zz_in * Pa_to_MPa, "bo", label=r"Num. $\sigma_{zz}$ (inner)", markersize=4, alpha=0.6)
plt.plot(r_in, sigma_zz_ana_in * Pa_to_MPa, "b-", label=r"Ana. $\sigma_{zz}$ (inner)", linewidth=1.5)
plt.xlabel("Radius (m)", fontsize=12)
plt.ylabel("Stress (MPa)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "stress_comparison_inner.png"), dpi=300)
plt.close()

plt.figure(figsize=(10, 7))
plt.plot(r_in, epsilon_rr_in, "ro", label=r"Num. $\varepsilon_{rr}$ (inner)", markersize=4, alpha=0.6)
plt.plot(r_in, epsilon_rr_ana_in, "r-", label=r"Ana. $\varepsilon_{rr}$ (inner)", linewidth=1.5)
plt.plot(r_in, epsilon_tt_in, "go", label=r"Num. $\varepsilon_{\theta\theta}$ (inner)", markersize=4, alpha=0.6)
plt.plot(r_in, epsilon_tt_ana_in, "g-", label=r"Ana. $\varepsilon_{\theta\theta}$ (inner)", linewidth=1.5)
plt.plot(r_in, epsilon_zz_in, "bo", label=r"Num. $\varepsilon_{zz}$ (inner)", markersize=4, alpha=0.6)
plt.plot(r_in, epsilon_zz_ana_in, "b-", label=r"Ana. $\varepsilon_{zz}$ (inner)", linewidth=1.5)
plt.xlabel("Radius (m)", fontsize=12)
plt.ylabel("Strain (/)", fontsize=12)
plt.grid(True, linestyle="--", alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "strain_comparison_inner.png"), dpi=300)
plt.close()

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
err_rr_o = np.sqrt(np.mean((sigma_rr_o - sigma_rr_ana_o) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_o**2))
err_tt_o = np.sqrt(np.mean((sigma_tt_o - sigma_tt_ana_o) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_o**2))
err_zz_o = np.sqrt(np.mean((sigma_zz_o - sigma_zz_ana_o) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_o**2))
err_eps_rr_o = np.sqrt(np.mean((epsilon_rr_o - epsilon_rr_ana_o) ** 2)) / np.sqrt(np.mean(epsilon_rr_ana_o**2))
err_eps_tt_o = np.sqrt(np.mean((epsilon_tt_o - epsilon_tt_ana_o) ** 2)) / np.sqrt(np.mean(epsilon_tt_ana_o**2))
if np.allclose(epsilon_zz_ana_o, 0.0):
    err_eps_zz_o = np.sqrt(np.mean((epsilon_zz_o - epsilon_zz_ana_o) ** 2))
else:
    err_eps_zz_o = np.sqrt(np.mean((epsilon_zz_o - epsilon_zz_ana_o) ** 2)) / np.sqrt(np.mean(epsilon_zz_ana_o**2))

err_rr_in = np.sqrt(np.mean((sigma_rr_in - sigma_rr_ana_in) ** 2)) / np.sqrt(np.mean(sigma_rr_ana_in**2))
err_tt_in = np.sqrt(np.mean((sigma_tt_in - sigma_tt_ana_in) ** 2)) / np.sqrt(np.mean(sigma_tt_ana_in**2))
err_zz_in = np.sqrt(np.mean((sigma_zz_in - sigma_zz_ana_in) ** 2)) / np.sqrt(np.mean(sigma_zz_ana_in**2))
err_eps_rr_in = np.sqrt(np.mean((epsilon_rr_in - epsilon_rr_ana_in) ** 2)) / np.sqrt(np.mean(epsilon_rr_ana_in**2))
err_eps_tt_in = np.sqrt(np.mean((epsilon_tt_in - epsilon_tt_ana_in) ** 2)) / np.sqrt(np.mean(epsilon_tt_ana_in**2))
if np.allclose(epsilon_zz_ana_in, 0.0):
    err_eps_zz_in = np.sqrt(np.mean((epsilon_zz_in - epsilon_zz_ana_in) ** 2))
else:
    err_eps_zz_in = np.sqrt(np.mean((epsilon_zz_in - epsilon_zz_ana_in) ** 2)) / np.sqrt(np.mean(epsilon_zz_ana_in**2))

errors = {
    "L2_error_sigma_rr_outer": {
        "numerical": float(err_rr_o),
        "reference": 0.0,
        "abs_error": float(err_rr_o),
        "rel_error": float(err_rr_o),
    },
    "L2_error_sigma_tt_outer": {
        "numerical": float(err_tt_o),
        "reference": 0.0,
        "abs_error": float(err_tt_o),
        "rel_error": float(err_tt_o),
    },
    "L2_error_sigma_zz_outer": {
        "numerical": float(err_zz_o),
        "reference": 0.0,
        "abs_error": float(err_zz_o),
        "rel_error": float(err_zz_o),
    },
    "L2_error_strain_rr_outer": {
        "numerical": float(err_eps_rr_o),
        "reference": 0.0,
        "abs_error": float(err_eps_rr_o),
        "rel_error": float(err_eps_rr_o),
    },
    "L2_error_strain_tt_outer": {
        "numerical": float(err_eps_tt_o),
        "reference": 0.0,
        "abs_error": float(err_eps_tt_o),
        "rel_error": float(err_eps_tt_o),
    },
    "L2_error_strain_zz_outer": {
        "numerical": float(err_eps_zz_o),
        "reference": 0.0,
        "abs_error": float(err_eps_zz_o),
        "rel_error": float(err_eps_zz_o),
    },
    "L2_error_sigma_rr_inner": {
        "numerical": float(err_rr_in),
        "reference": 0.0,
        "abs_error": float(err_rr_in),
        "rel_error": float(err_rr_in),
    },
    "L2_error_sigma_tt_inner": {
        "numerical": float(err_tt_in),
        "reference": 0.0,
        "abs_error": float(err_tt_in),
        "rel_error": float(err_tt_in),
    },
    "L2_error_sigma_zz_inner": {
        "numerical": float(err_zz_in),
        "reference": 0.0,
        "abs_error": float(err_zz_in),
        "rel_error": float(err_zz_in),
    },
    "L2_error_strain_rr_inner": {
        "numerical": float(err_eps_rr_in),
        "reference": 0.0,
        "abs_error": float(err_eps_rr_in),
        "rel_error": float(err_eps_rr_in),
    },
    "L2_error_strain_tt_inner": {
        "numerical": float(err_eps_tt_in),
        "reference": 0.0,
        "abs_error": float(err_eps_tt_in),
        "rel_error": float(err_eps_tt_in),
    },
    "L2_error_strain_zz_inner": {
        "numerical": float(err_eps_zz_in),
        "reference": 0.0,
        "abs_error": float(err_eps_zz_in),
        "rel_error": float(err_eps_zz_in),
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)


delta_exp=epsilon_tt_o[0]*R_mid



""" print(abs(delta_exp-delta)/delta, "erreur en r=R_mid")
print((delta_exp), "delta_exp")
print(epsilon_tt_o[0], "epsilon_tt at r=R_mid")
print(sigma_rr_o[0])
print(r_o[0])
print(sigma_rr_in[-1])
print(r_in[-1]) """


target_r_interface = 0.04005
diffs = np.abs(r_o - target_r_interface)
idx_closest = np.argmin(diffs)

if diffs[idx_closest] < 1e-5:
    # Utiliser le tableau déjà filtré (mid-plane) et trié pour éviter l'erreur d'index
    sigma_rr_o_at_interface = sigma_rr_o[idx_closest]
    print(f"[EXTRACT] Valeur de sigma_rr_o à r={r_o[idx_closest]:.6f} (cible {target_r_interface}) : {sigma_rr_o_at_interface:.6e} Pa")

print(delta_exp)
print(gap)
P_interface = (gap + delta_exp) / (R_mid * ((1 / E_in) * (R_mid**2 + R_int**2 * (1 + nu_in)) / (R_mid**2 - R_int**2) + (1 / (E_o * (R_o**2 - R_mid**2))) * (R_mid**2 + R_o**2 - nu_o * (2 * R_mid**2 - R_o**2))))  # Pa on the interface (inward)
print(f"{P_interface:.6e} P_interface recalculated from gap and delta_exp")
print(f"{sigma_rr_o_at_interface:.6e} sigma_rr_o at r=R_mid")
print(abs((-sigma_rr_o_at_interface-P_interface)/P_interface), "relative error en r=R_mid")


print("\n[INFO] Non-regression completed.\n")
