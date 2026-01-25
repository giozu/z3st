#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: stress_strain_curve_damage

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
Lx, Ly = 1.0, 0.500         # m (geometry dimensions)
E = 200e9                   # (Pa) Young modulus
nu = 0.3                    # Poisson modulus

n_fe = 40                   # number of vertical finite elements (from mesh.geo)

y_target, mask_tol = (
    Ly / 2,
    1 / (2 * n_fe),
)                           # m, m, m (extraction plane selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
STRESSES_REF = [ 0, 1e2, 1e4, 1e6, 1e8 ]                            # (Pa)
# IMPOSED STRESS: (1 - nu**2) accounts for the 2D PLANE STRAIN conditions (epsilon_zz = 0)
STRAINS_REF = [(1 - nu**2) / E * sigma for sigma in STRESSES_REF]   # (/)
U_X_REF = [eps * Lx for eps in STRAINS_REF]                         # (m)

TOLERANCE = 1e-2            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

strains = []
stresses = []
displacements = []

for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] Processing {os.path.basename(vtufile)}")

    # Stress extraction - usa vtufile (minuscolo)
    x_S, y_S, z_S, S_all = extract_field(vtufile, field_name="Stress_steel (cells)")
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
plt.yscale("linear")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "stress_strain_curve.png"))
print("[INFO] stress_strain_curve.png saved\n")