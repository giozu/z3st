#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: d_box_knotch_2D

non-regression script
---------------------

"""

import os
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

# Knotch
W_notch = 0.02
D_notch = 0.10
Y_tip   = Ly - D_notch
x_target = Lx / 2

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))
sigma_c = float(mat_data.get('sigma_c'))

# Phase-field / damage
lc = 0.002  
Gc_ref = (3.0 * sigma_c**2 * lc) / (2.0 * E)

TOLERANCE = 1e-2 

# --.. ..- .-.. .-.. --- Data --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# Damage
x_d, y_d, _, D_all = extract_field(VTU_FILE, field_name="Damage")
d_max = np.max(D_all)

# Stress
x_s, y_s, _, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
sigma_xx_max = np.max(S_all[:, 0]) 

# Stress profile
mask_d = np.abs(x_d - x_target) < (Lx/200)
idx_d = np.argsort(y_d[mask_d])
y_prof = y_d[mask_d][idx_d]
D_prof = D_all[mask_d][idx_d]

mask_s = np.abs(x_s - x_target) < (Lx/200)
idx_s = np.argsort(y_s[mask_s])
y_s_prof = y_s[mask_s][idx_s]
sigma_yy_prof = S_all[mask_s, 4][idx_s]

y_slice = Y_tip - 2*lc
mask_x = np.abs(y_d - y_slice) < (Ly/100)
idx_x = np.argsort(x_d[mask_x])
x_slice_prof = x_d[mask_x][idx_x]
d_slice_prof = D_all[mask_x][idx_x]

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---

# PLOT 1:
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(x_slice_prof, d_slice_prof, 'r-o', markersize=4, label=f"Damage at y={y_slice:.2f}")
plt.axhline(1.0, color='k', linestyle='--', alpha=0.3)
plt.xlabel("x (m)")
plt.ylabel("Damage $d$")
plt.grid(True, ls=':')
plt.legend()

plt.subplot(1, 2, 2)
sc = plt.scatter(x_d, y_d, c=D_all, cmap='jet', s=2)
plt.colorbar(sc, label="Damage $d$")
plt.axhline(Y_tip, color='white', linestyle=':', alpha=0.5, label="Notch Tip Line")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title(f"Max Damage: {d_max:.3f}")
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "damage_check.png"), dpi=300)

# PLOT 2:
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Damge
ax1.plot(y_prof, D_prof, "r-", lw=2, label="Damage $d$")
ax1.fill_between(y_prof, D_prof, color='red', alpha=0.1)
ax1.axvline(Y_tip, color='k', ls=':', label="Notch tip")
ax1.set_ylabel("Damage $d$ [-]", color='red')
ax1.set_ylim(-0.05, 1.05)
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend()
ax1.set_title(rf"Z3ST analysis: centerline profile (x = {x_target:.2f} m)\n"
              rf"Steel: $\sigma_c$={sigma_c*1e-6:.0f} MPa, $G_c$={Gc_ref:.1f} J/m²")

# Stress
ax2.plot(y_s_prof, sigma_yy_prof * 1e-6, "b-", lw=2, label=r"$\sigma_{yy}$")
ax2.axhline(sigma_c * 1e-6, color='black', ls='--', lw=1, label=rf"$\sigma_c$ Limit")
ax2.axvline(Y_tip, color='k', ls=':', label="Notch Tip")
ax2.set_xlabel("Vertical position $y$ [m]")
ax2.set_ylabel(fr"Vertical stress $\sigma_yy$ [MPa]", color='blue')
ax2.grid(True, ls=":", alpha=0.6)
ax2.legend(loc='best')

plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "damage_stress_profile.png"), dpi=300)
print(f"[INFO] Detailed profiles saved in: output/damage_stress_profile.png")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "max_damage": {
        "numerical": float(d_max),
        "reference": 0.5, # Soglia indicativa per innesco cricca
        "rel_error": float(abs(d_max - 0.5)) if d_max < 0.5 else 0.0
    },
    "max_stress_xx": {
        "numerical": float(sigma_xx_max),
        "reference": sigma_c,
        "rel_error": float(abs(sigma_xx_max - sigma_c)/sigma_c)
    }
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print(f"\n[INFO] Max damage: {d_max:.4f}")
print(f"[INFO] Max sigma_xx: {sigma_xx_max*1e-6:.2f} MPa")

print("[INFO] non-regression completed.\n")