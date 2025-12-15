#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 28_stress_strain_curve_damaged

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

if len(VTU_FILES) == 0:
    raise RuntimeError("[ERROR] No VTU files found in output/. Expected fields_*.vtu")

OUT_JSON = os.path.join(OUTPUT_DIR, "non-regression.json")

# geometry and material
Lx, Ly, Lz = 0.100, 0.100, 0.004  # m
E = 200e9  # (Pa) Young modulus

y_target, z_target, mask_tol = (
    Ly / 2,
    Lz / 2,
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
print(f"[INFO] Target z-plane for extraction: z = {z_target:.4e} m")

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
    x_s, y_s, z_s, s = extract_stress(vtufile, component="all", return_coords=True, prefer="cells")
    _, sigma_xx = average_section(
        x_s,
        y_s,
        z_s,
        s["xx"],
        y_target,
        z_target,
        mask_tol,
        decimals=5,
        label="sigma_xx",
    )
    sigmas.append(float(np.mean(sigma_xx)))

    # Displacement extraction
    x_d, y_d, z_d, u = extract_displacement(vtufile)
    ux = u[:, 0]

    _, ux_section = average_section(
        x_d,
        y_d,
        z_d,
        ux,
        y_target,
        z_target,
        mask_tol,
        decimals=5,
        label="ux",
    )
    u_max = np.max(ux_section)
    displacements.append(u_max)

    # Strain extraction
    eps_eng = u_max / Lx
    strains.append(eps_eng)

    print(f"  → σ_xx = {sigmas[-1]:.3e} Pa")
    print(f"  → ε_xx = {strains[-1]:.3e}")

# Stress–strain curve output
print("\n--. stress-strain values --..")
for e, s in zip(strains, sigmas):
    print(f"ε = {e:.3e}   σ = {s:.3e}")

plt.figure(figsize=(7, 5))
plt.plot(strains, sigmas, "--o", lw=2, label="Numerical")
# plt.plot(strain_ref_np, sigmas_ref_np, "-", lw=2, label='Analytical')
plt.xlabel(r"strain $\epsilon_{xx}$ (/)")
plt.ylabel(r"stress $\sigma_{xx}$ (Pa)")
plt.grid(True)
plt.title("Stress-strain curve")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "stress_strain_curve.png"))
print("[INFO] stress_strain_curve.png saved\n")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
strains_np = np.array(strains, dtype=float)
stresses_np = np.array(sigmas, dtype=float)

mask_zero = np.isclose(strain_ref_np, 0.0)

rel_error = np.empty_like(strains_np)
rel_error_s = np.empty_like(stresses_np)

rel_error[~mask_zero] = np.abs(strains_np[~mask_zero] - strain_ref_np[~mask_zero]) / np.abs(
    strain_ref_np[~mask_zero]
)
rel_error[mask_zero] = np.abs(strains_np[mask_zero] - strain_ref_np[mask_zero])
rel_error = rel_error.tolist()

rel_error_s[~mask_zero] = np.abs(strains_np[~mask_zero] - strain_ref_np[~mask_zero]) / np.abs(
    strain_ref_np[~mask_zero]
)
rel_error_s[mask_zero] = np.abs(strains_np[mask_zero] - strain_ref_np[mask_zero])
rel_error_s = rel_error_s.tolist()

errors = {
    "sigma_xx": {
        "numerical": strains_np.tolist(),
        "reference": strain_ref_np.tolist(),
        "abs_error": abs(strains_np - strain_ref_np).tolist(),
        "rel_error": rel_error,
    },
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
