#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: stress_strain_curve_damage

non-regression script
---------------------

This script processes multiple VTU files (fields_0000.vtu, ...),
extracts σ_xx and ε_xx for each step, and reconstructs the stress-strain curve.
"""

import os, yaml, re
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))

MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

# Geometry, material, boundary conditions:
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E  = float(mat_data.get('E'))
nu = float(mat_data.get('nu'))
sigma_c = float(mat_data.get('sigma_c'))

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()
W_notch = float(re.search(rf'W_notch\s*=\s*([\d\.]+);', content).group(1))
D_notch = float(re.search(rf'D_notch\s*=\s*([\d\.]+);', content).group(1))

print(f"[INFO] Geometry loaded: Lx = {Lx} m, Ly = {Ly} m")
print(f"[INFO] Material loaded: E = {E:.2e} Pa, nu = {nu}")
print(f"[INFO]                : sigma_c = {sigma_c:.2e} Pa")
print(f"[INFO] W_notch loaded: W_notch = {W_notch}")

# geometry and material
y_target, mask_tol = (Ly - D_notch, 1 / (2 * 40)) # m, m, m (extraction line selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
STRESSES_REF = [0.0, 1.0e6, 1.0e7, 1.0e8, 2.0e8, 3.0e8]             # (Pa)
# IMPOSED STRESS: (1 - nu**2) accounts for the 2D PLANE STRAIN conditions (epsilon_zz = 0)
STRAINS_REF = [(1 - nu**2) / E * sigma for sigma in STRESSES_REF]   # (/)
U_X_REF = [eps * Lx for eps in STRAINS_REF]                         # (m)

TOLERANCE = 1e-2            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
print(f"[INFO] Target y-plane for extraction: y = {y_target:.4e} m")

strains = []
stresses = []
displacements = []
d_max_list = []
h_max_list = []

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

    # Damage
    x_d, y_d, _, D_all = extract_field(vtufile, field_name="Damage")
    d_max_list.append(np.max(D_all))

    # Crack driving force
    _, _, _, H_all = extract_field(vtufile, field_name="CrackDrivingForce")
    h_max_list.append(np.max(H_all))

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

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
# Sigma-epsilon
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

# Damage
steps = np.arange(len(VTU_FILES))
fig, ax1 = plt.subplots(figsize=(9, 6))

ax1.set_xlabel('Step')
ax1.set_ylabel('Damage $D$ / History $H$', color='tab:red')
lns1 = ax1.plot(steps, d_max_list, 'r-o', lw=2, label='Max Damage $D$')
lns2 = ax1.plot(steps, h_max_list / np.max(h_max_list) if np.max(h_max_list)>0 else h_max_list, 
                'g--', lw=1.5, label='Normalized $H$')
ax1.set_ylim(-0.05, 1.1)
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid(True, ls=':', alpha=0.6)

ax2 = ax1.twinx()
ax2.set_ylabel(r'Stress $\sigma_{xx}$ [MPa]', color='tab:blue')
lns3 = ax2.plot(steps, np.array(stresses)*1e-6, 'b-s', lw=2, label=r'$\sigma_{xx}$ Mean at Tip')
ax2.axhline(sigma_c * 1e-6, color='black', ls=':', alpha=0.4, label=r'Critical $\sigma_c$')
ax2.tick_params(axis='y', labelcolor='tab:blue')

lns = lns1 + lns2 + lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='center left', frameon=True, shadow=True)

plt.title(f"Z3ST: Crack initiation & softening\nNotch tip evolution (y={y_target:.3f}m)")
fig.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "damage_evolution.png"), dpi=300)
