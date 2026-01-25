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

# Parameters
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))
sigma_c = float(mat_data.get('sigma_c'))

TOLERANCE = 1e-2 

# --.. ..- .-.. .-.. --- data --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# Damage
x_d, y_d, z_d, D_all = extract_field(VTU_FILE, field_name="Damage")
d_max = np.max(D_all)

# Stress
x_s, y_s, z_s, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
sigma_xx_max = np.max(S_all[:, 0]) 

# Crack profile
y_target = Ly * 0.8 
mask = np.abs(y_d - y_target) < (Ly/50)
sort_idx = np.argsort(x_d[mask])

x_profile = x_d[mask][sort_idx]
d_profile = D_all[mask][sort_idx]

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(x_profile, d_profile, 'r-o', markersize=4, label="Damage profile")
plt.axhline(1.0, color='k', linestyle='--', alpha=0.5)
plt.xlabel("x (m)")
plt.ylabel("Damage $d$")
plt.title(f"Damage profile at y = {y_target:.2f}m")
plt.grid(True, ls=':')
plt.legend()

# Distribution, scatter plot
plt.subplot(1, 2, 2)
sc = plt.scatter(x_d, y_d, c=D_all, cmap='jet', s=5)
plt.colorbar(sc, label="Damage $d$")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title(f"Max Damage: {d_max:.3f}")

plt.tight_layout()
plot_path = os.path.join(CASE_DIR, "output", "damage_check.png")
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "max_damage": {
        "numerical": float(d_max),
        "reference": 1.0,
        "rel_error": float(abs(d_max - 1.0))
    },
    "max_stress_xx": {
        "numerical": float(sigma_xx_max),
        "reference": sigma_c,
        "rel_error": float(abs(sigma_xx_max - sigma_c)/sigma_c)
    }
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print(f"\n[INFO] Danno Massimo rilevato: {d_max:.4f}")
print("[INFO] non-regression completed.\n")