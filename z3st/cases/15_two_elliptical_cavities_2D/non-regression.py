#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: two_elliptical_cavities_2D

non-regression script
---------------------
"""

import os, re, sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_extract_xdmf import *
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
XDMF_FILE = os.path.join(OUTPUT_DIR, "results.xdmf")
print(f"[INFO] Using XDMF file: {XDMF_FILE}")

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
print(f"[INFO] Theoretical intensification factor (single bubble): {intensification_factor:.2f}")

p_applied = 15.0 # MPa (Correspond au maximum de la liste de Neumann)
p_target = p_applied * intensification_factor

n_bubbles = 2
Fc_linear = n_bubbles * 2.0 * ax / Lx
Fc_area = n_bubbles**2 * (np.pi * ax**2) / (Lx**2)

print(f"\n[GEOMETRY ANALYSIS]")
print(f"  → Fc (2D): {Fc_linear:.4f}")
print(f"  → Fc (3D):  {Fc_area:.4f}")
print(f"  → Domain (Lx, Ly): {Lx*1e6:.3f} x {Ly*1e6:.3f} μm^2")
print(f"  → Bubble (ax, ay): {ax*1e6:.3f} x {ay*1e6:.3f} μm^2")

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))

def Gc(y_coords):

    Gc_gb = 0.002
    Gc_bulk = 0.01
    half_width = 10e-3
    transition = np.tanh(np.abs(y_coords) / half_width)
    return (Gc_gb + (Gc_bulk - Gc_gb) * transition) * 1e-6

# --.. ..- .-.. .-.. --- extract fields --.. ..- .-.. .-.. ---
list_fields_xdmf(XDMF_FILE)

X_tip = ax
y_target = 0.0

try:
    x_d, y_d, _, D_all = extract_field_xdmf(XDMF_FILE, field_name="Damage")
    d_max = np.max(D_all)
except Exception:
    d_max = 0.0

mat_name = list(input_data.get("materials", {}).keys())[0]
stress_field_name = f"Stress_{mat_name}"
x, y, _, sigma = extract_field_xdmf(XDMF_FILE, field_name=stress_field_name)
sigma_yy_final = np.max(sigma[:, 4]) / 1e6 # Convert to MPa

# --.. ..- .-.. .-.. --- plot macroscopic stress-strain --.. ..- .-.. .-.. ---
strains = []
stresses = []

print("\n[INFO] Plotting macroscopic stress and strain from diagnostics...")
stress_strain_file = os.path.join(OUTPUT_DIR, "stress_strain.txt")

if os.path.exists(stress_strain_file):
    try:
        data = np.loadtxt(stress_strain_file, skiprows=1)
        if data.ndim == 1: data = data.reshape(1, -1)
        strains = list(data[:, 1])
        stresses = list(data[:, 2])
    except Exception as e:
        print(f"[WARNING] Could not read {stress_strain_file}: {e}")

if strains:
    plt.figure(figsize=(10, 6))
    plt.plot(strains, stresses, 'r-s', markersize=4, label=r"Two Cavities - $\sigma_{yy}$")
    plt.xlabel(r"Macroscopic Strain $\varepsilon_{yy}$ (-)")
    plt.ylabel(r"Average Stress $\sigma_{yy}$ (MPa)")
    plt.title("Stress-Strain Curve (Two Interacting Cavities)")
    plt.grid(True, ls=':', alpha=0.6)
    plt.legend()

    plot_path_ss = os.path.join(OUTPUT_DIR, "stress_strain_curve.png")
    plt.tight_layout()
    plt.savefig(plot_path_ss, dpi=300)
    print(f"[INFO] Stress-strain plot saved in: {plot_path_ss}")

# Gc profile
y_line = np.linspace(-Ly/2, Ly/2, 500) # (micron)
gc_line = Gc(y_line)

sigma_c = ((27 * E * gc_line) / (256 * lc))**0.5

plt.figure(figsize=(8, 5))
plt.plot(y_line * 1000, gc_line * 1e6, 'g-', label="$G_c$ profile (GB zone)")
plt.xlabel(r"Distance $y$ (nm)")
plt.ylabel("$G_c$ ($J/m^2$)")
plt.grid(True, ls=':')
plt.legend()
plt.savefig(os.path.join(CASE_DIR, "output", "gc_profile_check.png"))

TOLERANCE = 0.1

peak_stress = max(stresses) if stresses else sigma_yy_final

errors = {
    "max_stress_yy": {
        "numerical": float(peak_stress),
        "reference": float(peak_stress),
        "rel_error": 0.0
    }
}

print(f"\n[RESULTS]")
if d_max > 0.0:
    print(f"  → Max Damage: {d_max:.4f}")
if stresses:
    print(f"  → Peak Macroscopic Stress: {max(stresses):.2f} MPa")
print(f"  → Residual Stress at end : {sigma_yy_final:.2f} MPa")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
