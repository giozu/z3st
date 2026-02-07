#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: box_elliptical_cavity_2D

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
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/oxide.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")

# Input
with open(INPUT_FILE, 'r') as f:
    input_data = yaml.safe_load(f)
lc = float(input_data.get("damage", {}).get("lc", 0.001))

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
ax = float(geom_data.get('ax'))
ay = float(geom_data.get('ay'))
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

theta = 50 * np.pi / 180 # angle in radians, semi-dihedral angle
ay_ax = (1 - np.cos(theta)) / np.sin(theta)
intensification_factor = 2 / ay_ax - 1 # theoretical pressure intensification
print(f"[INFO] Theoretical intensification factor: {intensification_factor:.2f}")

p_applied = 1.0 # MPa
p_target = p_applied * intensification_factor

n_bubbles = 2
Fc_linear = n_bubbles * 2.0 * ax / Lx
Fc_area = n_bubbles**2 * (np.pi * ax**2) / (Lx**2)

print(f"\n[GEOMETRY ANALYSIS]")
print(f"  → Fc (2D): {Fc_linear:.4f}")
print(f"  → Fc (3D):  {Fc_area:.4f}")
print(f"  → Domain (Lx, Ly): {Lx:.3f} x {Ly:.3f} μm^2")
print(f"  → Bubble (ax, ay): {ax:.3f} x {ay:.3f} μm^2")

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))

def Gc(y_coords):
    Gc_gb = 2.0     # J/m2
    Gc_bulk = 100.0 # J/m2
    half_width = 1e-3 # 1 nm in micron
    transition = np.tanh(np.abs(y_coords) / half_width)
    return (Gc_gb + (Gc_bulk - Gc_gb) * transition) * 1e-6

# --.. ..- .-.. .-.. --- extract fields --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

X_tip = ax
y_target = 0.0

x_d, y_d, _, D_all = extract_field(VTU_FILE, field_name="Damage")
d_max = np.max(D_all)

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
plt.ylabel(r"Stress$ (MPa)")
plt.grid(True, ls=':', alpha=0.6)
plt.legend()

plot_path = os.path.join(CASE_DIR, "output", "stress_profile_tip.png")
plt.tight_layout()
plt.savefig(plot_path, dpi=300)
print(f"[INFO] Plot saved in: {plot_path}")

# Gc profile
# y_line = np.linspace(-Ly/2, Ly/2, 500) # (micron)
# gc_line = Gc(y_line)

# sigma_c = ((27 * E * gc_line) / (256 * lc))**0.5

# plt.figure(figsize=(8, 5))
# plt.plot(y_line * 1000, gc_line * 1e6, 'g-', label="$G_c$ profile (GB zone)")
# plt.xlabel(r"Distance $y$ (nm)")
# plt.ylabel("$G_c$ ($J/m^2$)")
# plt.grid(True, ls=':')
# plt.legend()
# plt.savefig(os.path.join(CASE_DIR, "output", "gc_profile_check.png"))

TOLERANCE = 0.1
errors = {
    "max_damage": {
        "numerical": float(d_max),
        "reference": 1.0, 
        "rel_error": float(abs(d_max - 1.0))
    },
    "max_stress_yy": {
        "numerical": float(sigma_yy_max),
        "reference": float(p_target),
        "rel_error": float(abs(sigma_yy_max - p_target)/p_target)
    }
}

print(f"\n[RESULTS]")
# print(f"  → Max Damage: {d_max:.4f}")
print(f"  → Max Stress: {sigma_yy_max:.2f} MPa (Target: {p_target:.2f} MPa)")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
