#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: single_elliptical_cavity_2D

non-regression script
---------------------
"""

import os, re
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import *

# Units of measure for micromechanics
# Length: micrometer (μm)
# Time: second (s)
# Mass: kilogram (kg)
# Force: micronewton (μN)
# Pressure: megapascal (MPa)
# Energy: picojoule (pJ)

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)

OUTPUT_DIR = os.path.join(CASE_DIR, "output")
vtu_files = [f for f in os.listdir(OUTPUT_DIR) if f.startswith("fields_") and f.endswith(".vtu")]
if vtu_files:
    latest_vtu = sorted(vtu_files)[-1]
    VTU_FILE = os.path.join(OUTPUT_DIR, latest_vtu)
    print(f"[INFO] Using VTU file: {VTU_FILE}")
else:
    VTU_FILE = os.path.join(OUTPUT_DIR, "fields.vtu")

OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/oxide.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
ax = float(geom_data.get('ax'))
ay = float(geom_data.get('ay'))
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

theta = 30 * np.pi / 180 # angle in radians, semi-dihedral angle
ay_ax = (1 - np.cos(theta)) / np.sin(theta)
intensification_factor = 2 / ay_ax - 1 # theoretical pressure intensification
print(f"[INFO] Theoretical intensification factor: {intensification_factor:.2f}")

p_applied = 1.0 # MPa
p_target = p_applied * intensification_factor

print(f"\n[GEOMETRY ANALYSIS]")
print(f"  → Domain (Lx, Ly): {Lx:.3f} x {Ly:.3f} μm^2")
print(f"  → Bubble (ax, ay): {ax:.3f} x {ay:.3f} μm^2")

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))

# --.. ..- .-.. .-.. --- extract fields --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

X_tip = ax
y_target = 0.0

x, y, _, sigma = extract_field(VTU_FILE, field_name="Stress_solid (cells)")
sigma_yy_max = np.max(sigma[:, 4]) 

# mask = (np.abs(y - y_target) < (Ly/500)) & (x >= X_tip)
mask = (np.abs(y - y_target) < (Ly/100))
idx_line = np.argsort(x[mask])
x_line = x[mask][idx_line]
sigma_yy_line = sigma[mask, 4][idx_line]

# --.. ..- .-.. .-.. --- plot --.. ..- .-.. .-.. ---
plt.figure(figsize=(10, 6))

plt.plot(x_line, sigma_yy_line, 'b-o', markersize=4, label=r"$\sigma_{yy}$")
# plt.axvline(X_tip, color='r', linestyle='--', label="Bubble tip")

plt.xlabel(r"Distance $x$ ($\mu$m)")
plt.ylabel(r"Stress (MPa)")
plt.grid(True, ls=':', alpha=0.6)
plt.legend()

plot_path = os.path.join(CASE_DIR, "output", "stress_profile_tip.png")
plt.tight_layout()
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

errors = {
    "max_stress_yy": {
        "numerical": float(sigma_yy_max),
        "reference": float(p_target),
        "rel_error": float(abs(sigma_yy_max - p_target)/p_target)
    }
}

TOLERANCE = 1.0e-2
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)

print(f"\n[RESULTS]")
print(f"  → Max Stress: {sigma_yy_max:.2f} MPa (Target: {p_target:.2f} MPa)")
