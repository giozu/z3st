#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 28_stress_strain_curve_displacement

non-regression script
---------------------

This script processes multiple VTU files (fields_0000.vtu, ...),
extracts σ_xx and ε_xx for each step, and reconstructs the stress-strain curve.
"""

import os
from glob import glob

import matplotlib.pyplot as plt
import numpy as np

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
OUT_JSON = os.path.join(OUTPUT_DIR, "non-regression.json")

# geometry and material
Lx, Ly = 0.1, 1.0  # m (geometry dimensions)
E = 200e9  # (Pa) Young modulus

y_target, mask_tol = (
    Ly / 2,
    0.01,
)  # m, m, m (extraction plane selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
SIGMAS_REF = [0, 1e2, 1e4, 1e6, 1e8]  # (Pa)
STRAINS_REF = [sigma / E for sigma in SIGMAS_REF]

strain_ref_np = np.array(STRAINS_REF, dtype=float)
sigmas_ref_np = np.array(SIGMAS_REF, dtype=float)

TOLERANCE = 1e-2  # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

# stress–strain arrays
strains = []
sigmas = []
displacements = []

print(f"[INFO] Found {len(VTU_FILES)} VTU files")
print(f"[INFO] Extracting σ_xx and ε_xx on plane y={y_target:.3e}, z={z_target:.3e}")


for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] Processing {os.path.basename(vtufile)}")

    # --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
    list_fields(vtufile)

    # Stress extraction
    x_S, y_S, z_S, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
    mask = np.abs(y_S - y_target) < mask_tol
    sort_idx = np.argsort(x_S[mask])

    x_s = x_S[mask][sort_idx]

    sigma_xx = S_all[mask, 0][sort_idx]
    sigmas.append(float(np.mean(sigma_xx)))

    # Displacement extraction
    x_u_all, y_u_all, z_u_all, u_all = extract_field(VTU_FILE, field_name="Displacement")
    mask_u = np.abs(y_u_all - y_target) < mask_tol
    sort_idx_u = np.argsort(x_u_all[mask_u])

    x_u_line = x_u_all[mask_u][sort_idx_u]
    ux_line = u_all[mask_u, 0][sort_idx_u]

    u_max = np.max(ux_line)
    displacements.append(u_max)

    # Strain extraction
    x_s_all, y_s_all, z_s_all, S_all = extract_field(VTU_FILE, field_name="Strain (cells)")
    mask_s = np.abs(y_s_all - y_target) < mask_tol
    sort_idx_s = np.argsort(x_s_all[mask_s])

    x_s_line = x_s_all[mask_s][sort_idx_s]
    sigma_xx = S_all[mask_s, 0][sort_idx_s]


    eps_eng = float(np.mean(eps_xx_section))
    strains.append(eps_eng)

    print(f"  → σ_xx = {sigmas[-1]:.3e} Pa")
    print(f"  → ε_xx = {strains[-1]:.3e}")
    print(f"  → u_x =  {u_max:.3e} m")

# Stress–strain curve output
print("\n--. stress-strain-displacement values --..")
for e, s, u in zip(strains, sigmas, displacements):
    print(f"ε_xx = {e:.3e}\tσ_xx = {s:.3e}\tu_x = {u:.3e}")

plt.figure(figsize=(7, 5))
plt.plot(strains, sigmas, "--o", lw=2, label="Numerical")
plt.plot(strain_ref_np, sigmas_ref_np, "-", lw=2, label="Analytical")
plt.xlabel(r"strain $\epsilon_{xx}$ (/)")
plt.ylabel(r"stress $\sigma_{xx}$ (Pa)")
plt.grid(True)
plt.title("Stress-strain curve")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "stress_strain_curve.png"))
print("[INFO] stress_strain_curve.png saved\n")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
strains_np = np.array(strains, dtype=float)

mask_zero = np.isclose(strain_ref_np, 0.0)
rel_error = np.empty_like(strains_np)
rel_error[~mask_zero] = np.abs(strains_np[~mask_zero] - strain_ref_np[~mask_zero]) / np.abs(
    strain_ref_np[~mask_zero]
)
rel_error[mask_zero] = np.abs(strains_np[mask_zero] - strain_ref_np[mask_zero])
rel_error = rel_error.tolist()

errors = {
    "epsilon_xx": {
        "numerical": strains_np.tolist(),
        "reference": strain_ref_np.tolist(),
        "abs_error": abs(strains_np - strain_ref_np).tolist(),
        "rel_error": rel_error,
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
