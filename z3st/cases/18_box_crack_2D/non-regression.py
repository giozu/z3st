#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: box_crack_2D

non-regression script
---------------------

"""

import os, re
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/high_carbon_steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")

# Phase-field / damage
with open(INPUT_FILE, 'r') as f:
    input_data = yaml.safe_load(f)
dmg_cfg = input_data.get("damage", {})
lc = float(dmg_cfg["lc"])

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()

Dn = float(re.search(rf'Dn\s*=\s*([\d\.]+);', content).group(1))

X_tip = Dn
y_target = Ly/2

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)

E = float(mat_data.get('E'))
nu = float(mat_data.get('nu'))
Gc = float(mat_data.get('Gc'))
# sigma_c = float(mat_data.get('sigma_c'))
sigma_c = ((27 * E * Gc) / (256 * lc))**0.5

TOLERANCE = 7e-2

# --.. ..- .-.. .-.. --- Data --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# Damage
x_d, y_d, _, D_all = extract_field(VTU_FILE, field_name="Damage")
d_max = np.max(D_all)

# Stress
x_s, y_s, _, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
sigma_yy_max = np.max(S_all[:, 4]) 

mask_horiz = np.abs(y_d - y_target) < (Ly/1500)
idx_h = np.argsort(x_d[mask_horiz])
x_prof = x_d[mask_horiz][idx_h]
D_prof = D_all[mask_horiz][idx_h]

mask_s = np.abs(y_s - y_target) < (Ly/600)
idx_s = np.argsort(x_s[mask_s])
x_s_prof = x_s[mask_s][idx_s]
sigma_yy_prof = S_all[mask_s, 4][idx_s]

x_slice = X_tip - 2*lc
mask_vert = np.abs(x_d - x_slice) < (Lx/500) 
idx_v = np.argsort(y_d[mask_vert]) 
y_slice_prof = y_d[mask_vert][idx_v]
d_slice_prof = D_all[mask_vert][idx_v]

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---

# PLOT 1:
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(y_slice_prof, d_slice_prof, 'r-', markersize=4, label=f"Damage at x={x_slice:.2f}")
plt.axhline(1.0, color='k', linestyle='--', alpha=0.3)
plt.xlabel("x (m)")
plt.ylabel("Damage $D$")
plt.grid(True, ls=':')
plt.legend()

plt.subplot(1, 2, 2)
sc = plt.scatter(x_d, y_d, c=D_all, cmap='jet', s=2)
plt.colorbar(sc, label="Damage $D$")
plt.axvline(X_tip, color='white', linestyle=':', alpha=0.5, label="Notch tip line")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title(f"Max damage: {d_max:.3f}")
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "damage_check.png"), dpi=300)

# PLOT 2:
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Damage
ax1.plot(x_prof, D_prof, "r-", lw=2, label="Damage $D$")
ax1.fill_between(x_prof, D_prof, color='red', alpha=0.1)
ax1.axvline(X_tip, color='k', ls=':', label="Notch tip")
ax1.set_ylabel("Damage $D$", color='red')
ax1.set_ylim(-0.05, 1.05)
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend()
ax1.set_title(f"Z3ST analysis: centerline profile (y = {y_target:.2f} m)\n"
              rf"Steel: $\sigma_c$={sigma_c*1e-6:.0f} MPa, $G_c$={Gc:.1f} J/m²")

# Stress
ax2.plot(x_s_prof, sigma_yy_prof * 1e-6, "b-", lw=2, label=r"$\sigma_{yy}$")
ax2.axhline(sigma_c * 1e-6, color='black', ls='--', lw=1, label=rf"$\sigma_c$")
ax2.axvline(X_tip, color='red', ls=':', label="Notch tip")
ax2.set_xlabel("Vertical position $y$ (m)")
ax2.set_ylabel(fr"Vertical stress $\sigma$ (MPa)", color='blue')
ax2.grid(True, ls=":", alpha=0.6)
ax2.legend(loc='best')

plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "damage_stress_profile.png"), dpi=300)
print(f"[INFO] Detailed profiles saved in: output/damage_stress_profile.png")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "max_damage": {
        "numerical": float(d_max),
        "reference": 1.0, 
        "rel_error": float(abs(d_max - 1.0))
    },
    "max_stress_yy": {
        "numerical": float(sigma_yy_max),
        "reference": sigma_c,
        "rel_error": float(abs(sigma_yy_max - sigma_c)/sigma_c)
    }
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print(f"\n[INFO] Max damage: {d_max:.4f}")
print(f"[INFO] Max sigma_yy: {sigma_yy_max*1e-6:.2f} MPa")

print("[INFO] non-regression completed.\n")