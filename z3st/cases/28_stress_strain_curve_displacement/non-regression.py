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
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
OUT_JSON = os.path.join(OUTPUT_DIR, "non-regression.json")

# geometry and material
Lx, Ly = 1.0, 0.010         # m (geometry dimensions)
E = 200e9                   # (Pa) Young modulus

n_fe = 40                   # number of vertical finite elements (from mesh.geo)

y_target, mask_tol = (
    Ly / 2,
    1 / (2 * n_fe),
)                           # m, m, m (extraction plane selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
U_X_REF = [0.0, 5e-11, 5e-9, 5e-7, 5e-5]                # (m)
STRESSES_REF = [ ux / Lx * E for ux in U_X_REF ]        # (Pa)
STRAINS_REF = [sigma / E for sigma in STRESSES_REF]     # (/)

strain_ref_np = np.array(STRAINS_REF, dtype=float)
stresses_ref_np = np.array(STRESSES_REF, dtype=float)

TOLERANCE = 1e-1            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

# stress–strain arrays
strains = []
sigmas = []
displacements = []

for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] Processing {os.path.basename(vtufile)}")

    # Stress extraction - usa vtufile (minuscolo)
    x_S, y_S, z_S, S_all = extract_field(vtufile, field_name="Stress_steel (cells)")
    mask = np.abs(y_S - y_target) < mask_tol
    sigmas.append(float(np.mean(S_all[mask, 0])))

    # Displacement extraction
    x_u_all, y_u_all, z_u_all, u_all = extract_field(vtufile, field_name="Displacement")
    mask_u = np.abs(y_u_all - y_target) < mask_tol
    u_max = np.max(u_all[mask_u, 0])
    displacements.append(u_max)

    # Strain extraction
    x_eps_all, y_eps_all, z_eps_all, E_all = extract_field(vtufile, field_name="Strain (cells)")
    mask_e = np.abs(y_eps_all - y_target) < mask_tol
    eps_eng = float(np.mean(E_all[mask_e, 0]))
    strains.append(eps_eng)

# Stress–strain curve output
print("\n--. stress-strain-displacement values --..")
for e, s, u in zip(strains, sigmas, displacements):
    print(f"ε_xx = {e:.3e}\tσ_xx = {s:.3e}\tu_x = {u:.3e}")

plt.figure(figsize=(7, 5))
plt.plot(strains, sigmas, "--o", lw=2, label="Numerical")
plt.plot(strain_ref_np, stresses_ref_np, "-", lw=2, label="Analytical")
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
stresses_np = np.array(sigmas, dtype=float)

# Strain
mask_zero = np.isclose(strain_ref_np, 0.0)
rel_error_E = np.zeros_like(strains_np)
rel_error_E[~mask_zero] = np.abs(strains_np[~mask_zero] - strain_ref_np[~mask_zero]) / np.abs(strain_ref_np[~mask_zero])
rel_error_E[mask_zero] = np.abs(strains_np[mask_zero] - strain_ref_np[mask_zero])

# Stress
rel_error_S = np.zeros_like(stresses_np)
mask_zero_S = np.isclose(stresses_ref_np, 0.0)
rel_error_S[~mask_zero_S] = np.abs(stresses_np[~mask_zero_S] - stresses_ref_np[~mask_zero_S]) / np.abs(stresses_ref_np[~mask_zero_S])
rel_error_S[mask_zero_S] = np.abs(stresses_np[mask_zero_S] - stresses_ref_np[mask_zero_S])

errors = {
    "epsilon_xx": {
        "numerical": strains_np.tolist(),
        "reference": strain_ref_np.tolist(),
        "abs_error": np.abs(strains_np - strain_ref_np).tolist(),
        "rel_error": rel_error_E.tolist(),
    },
    "sigma_xx": {
        "numerical": stresses_np.tolist(),
        "reference": stresses_ref_np.tolist(),
        "abs_error": np.abs(stresses_np - stresses_ref_np).tolist(),
        "rel_error": rel_error_S.tolist(),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
