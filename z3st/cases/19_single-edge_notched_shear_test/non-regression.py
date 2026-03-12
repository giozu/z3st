#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: single-edge notched shear test

non-regression script
---------------------

This script processes multiple VTU files (fields_0000.vtu, ...),
extracts σ and ε for each step, and reconstructs the stress-strain curve.
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

INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/high_carbon_steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

# Geometry, material, boundary conditions:
with open(INPUT_FILE, 'r') as f:
    input_data = yaml.safe_load(f)
dmg_cfg = input_data.get("damage", {})
lc = float(dmg_cfg["lc"])

with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E  = float(mat_data.get('E'))
nu = float(mat_data.get('nu'))
Gc = float(mat_data.get('Gc'))
sigma_c = np.sqrt(E * Gc / lc)**0.5

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()
Dn = float(re.search(rf'Dn\s*=\s*([\d\.]+);', content).group(1))

print(f"[INFO] Geometry loaded: Lx = {Lx} m, Ly = {Ly} m")
print(f"[INFO] Material loaded: E = {E:.2e} Pa, nu = {nu}")
print(f"[INFO]                : sigma_c = {sigma_c:.2e} Pa")
print(f"[INFO]                : Gc = {Gc:.2e} J/m2")
print(f"[INFO] Dn loaded      : Dn = {Dn}")

# geometry and material
x_target, mask_tol = (Dn, 1 / (2 * 40)) # m, m, m (extraction line selection and tolerance)

# --.. ..- .-.. .-.. --- analytic functions  --.. ..- .-.. .-.. ---
# Fracture energy, initial value
crack_length = Dn  # m
E_frac_i = Gc * crack_length
print(f'Theoretical initial fracture energy = {E_frac_i} J')

TOLERANCE = 1e-2            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
strains = []
stresses = []
displacements = []
d_max_list = []
h_max_list = []

for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] Processing {os.path.basename(vtufile)}")

    # Stress extraction - usa vtufile (minuscolo)
    x_S, y_S, z_S, S_all = extract_field(vtufile, field_name="Stress_steel (cells)")
    mask = np.abs(x_S - x_target) < mask_tol
    stresses.append(float(np.mean(S_all[mask, 4])))

    # Displacement extraction
    x_u_all, y_u_all, z_u_all, u_all = extract_field(vtufile, field_name="Displacement")
    mask_u = np.abs(x_u_all - x_target) < mask_tol
    u_max = np.max(u_all[mask_u, 1])
    displacements.append(u_max)

    # Strain extraction
    x_eps_all, y_eps_all, z_eps_all, E_all = extract_field(vtufile, field_name="Strain (cells)")
    mask_e = np.abs(x_eps_all - x_target) < mask_tol
    eps_eng = float(np.mean(E_all[mask_e, 4]))
    strains.append(eps_eng)

    # Damage
    x_d, y_d, _, D_all = extract_field(vtufile, field_name="Damage")
    d_max_list.append(np.max(D_all))

    # Crack driving force
    _, _, _, H_all = extract_field(vtufile, field_name="CrackDrivingForce")
    h_max_list.append(np.max(H_all))

# Stress–strain curve output
strains_np = np.array(strains, dtype=float)
stresses_np = np.array(stresses, dtype=float)
displ_x_np = np.array(displacements, dtype=float)

print("\n--. stress-strain-displacement values --..")
for e, s, u in zip(strains, stresses, displacements):
    print(f"ε_yy = {e:.3e}\tσ_yy = {s:.3e}\tu_x = {u:.3e}")

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
# crack profile
plt.figure(figsize=(8, 6))

x_line = Dn / 2.0
tol_x = Lx / 100.0
yc = Ly / 2.0

colors = plt.cm.jet(np.linspace(0, 1, len(VTU_FILES)))

for i, vtufile in enumerate(VTU_FILES):
    print(vtufile)
    x_d, y_d, _, D_all = extract_field(vtufile, field_name="Damage")
    mask_line = np.abs(x_d - x_line) < tol_x
    
    if np.any(mask_line):
        y_coords = y_d[mask_line]
        d_values = D_all[mask_line]
        
        sort_idx = np.argsort(y_coords)
        y_plot = y_coords[sort_idx]
        d_plot = d_values[sort_idx]
        
        label = f"Step {i}"
        plt.plot(y_plot, d_plot, color=colors[i], lw=1.5, label=label, alpha=0.8)

d_analytical = np.exp(-np.abs(y_plot - yc) / lc)
plt.plot(y_plot, d_analytical, "k--", lw=2, label=r"Theoretical $\exp(-|y-y_c|/l_c)$", zorder=10)

plt.xlabel("Coordinate (m)")
plt.ylabel("Damage $D$ (/)")
plt.title(f"Evolution of the crack profile at $x={x_line:.3f}$ m, $l_c = {lc}$ m")
plt.grid(True, ls=':', alpha=0.6)
plt.legend(loc='upper right', fontsize='small', ncol=2)
plt.tight_layout()
# plt.show()
plt.savefig(os.path.join(OUTPUT_DIR, "crack_profile_evolution.png"))
print(f"[INFO] crack_profile_evolution.png saved showing {len(VTU_FILES)} steps.")

# energy balance
energy_file = os.path.join("energies.txt")
if os.path.exists(energy_file):
    data = np.genfromtxt(energy_file, names=True, skip_header=0)
    
    plt.figure(figsize=(8, 5))
    plt.plot(data['Step'], data['E_el'], 'b-o', label='Elastic Energy ($E_{el}$)')
    plt.plot(data['Step'], data['E_frac'], 'r-s', label='Fracture Energy ($E_{frac}$)')
    plt.plot(data['Step'], data['E_tot'], 'k--', label='Total Energy ($E_{tot}$)')
    
    plt.xlabel('Step')
    plt.ylabel('Energy (J)')
    plt.title('Z3ST: Global energy balance')
    plt.grid(True, ls=':', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "energy_balance.png"))
    print("[INFO] energy_balance.png saved")

# Sigma-epsilon
plt.figure(figsize=(7, 5))
plt.plot(strains, stresses, "--o", lw=2, label="Numerical")
plt.xlabel(r"strain $\epsilon_{yy}$ (/)")
plt.ylabel(r"stress $\sigma_{yy}$ (Pa)")
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
ax2.set_ylabel(r'Stress $\sigma_{yy}$ (MPa)', color='tab:blue')
lns3 = ax2.plot(steps, np.array(stresses)*1e-6, 'b-s', lw=2, label=r'$\sigma_{yy}$ Mean at Tip')
ax2.axhline(sigma_c * 1e-6, color='black', ls=':', alpha=0.4, label=r'Critical $\sigma_c$')
ax2.tick_params(axis='y', labelcolor='tab:blue')

lns = lns1 + lns2 + lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='center left', frameon=True, shadow=True)

plt.title(f"Z3ST: Crack initiation & softening\nNotch tip evolution (x={x_target:.3f}m)")
fig.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "damage_evolution.png"), dpi=300)
