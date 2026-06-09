#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: stress_strain_curve_displacement

non-regression script
---------------------

This script processes multiple VTU files (fields_0000.vtu, ...),
extracts σ_xx and ε_xx for each step, and reconstructs the stress-strain curve.
"""

import os
import yaml
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

# geometry and material (loaded from YAML)
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly = float(geom["Lx"]), float(geom["Ly"])   # m (geometry dimensions)

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
mat_path = os.path.join(CASE_DIR, next(iter(inp["materials"].values())))
with open(mat_path) as f:
    mat = yaml.safe_load(f)
E = float(mat["E"])          # (Pa) Young modulus
nu = float(mat["nu"])        # (/) Poisson modulus

n_fe = 40                   # number of vertical finite elements (from mesh.geo)

y_target, mask_tol = (
    Ly / 2,
    1 / (2 * n_fe),
)                           # m, m, m (extraction plane selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
U_X_REF = [0.0, 5e-11, 5e-9, 5e-7, 5e-5]    # (m) imposed displacement
STRAINS_REF = [ux / Lx for ux in U_X_REF]   # (/)
# IMPOSED DISPLACEMENT: (1 - nu**2) accounts for the 2D PLANE STRAIN condition
STRESSES_REF = [ eps_xx * E / (1 - nu**2) for eps_xx in STRAINS_REF ] # (Pa)

TOLERANCE = 1e-3            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

# stress–strain arrays
strains = []
stresses = []
displacements = []

for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] Processing {os.path.basename(vtufile)}")

    # Stress extraction - usa vtufile (minuscolo)
    x_S, y_S, z_S, S_all = extract_field(vtufile, field_name="Stress (cells)")
    mask = np.abs(y_S - y_target) < mask_tol
    stresses.append(float(np.mean(S_all[mask, 0])))

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
strain_ref_np = np.array(STRAINS_REF, dtype=float)
stresses_ref_np = np.array(STRESSES_REF, dtype=float)
displ_x_ref_np = np.array(U_X_REF, dtype=float)

strains_np = np.array(strains, dtype=float)
stresses_np = np.array(stresses, dtype=float)
displ_x_np = np.array(displacements, dtype=float)

print("\n--. stress-strain-displacement values --..")
for e, s, u in zip(strains, stresses, displacements):
    print(f"ε_xx = {e:.3e}\tσ_xx = {s:.3e}\tu_x = {u:.3e}")

plt.figure(figsize=(7, 5))
plt.plot(strains, stresses, "--o", lw=2, label="Numerical")
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

# Displacement
rel_error_U = np.zeros_like(displ_x_np)
mask_zero_U = np.isclose(displ_x_ref_np, 0.0)
rel_error_U[~mask_zero_U] = np.abs(displ_x_np[~mask_zero_U] - displ_x_ref_np[~mask_zero_U]) / np.abs(displ_x_ref_np[~mask_zero_U])
rel_error_U[mask_zero_U] = np.abs(displ_x_np[mask_zero_U] - displ_x_ref_np[mask_zero_U])


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
    "u_xx": {
        "numerical": displ_x_np.tolist(),
        "reference": displ_x_ref_np.tolist(),
        "abs_error": np.abs(displ_x_np - displ_x_ref_np).tolist(),
        "rel_error": rel_error_U.tolist(),
    },
}

# --.. ..- .-.. .-.. --- pass/fail + regression --.. ..- .-.. .-.. ---
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
