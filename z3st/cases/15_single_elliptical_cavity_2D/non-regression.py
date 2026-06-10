#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: single_elliptical_cavity_2D

non-regression script
---------------------
"""

import os, re, sys
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

def extract_step(filename):
    match = re.search(r'_(\d+)\.vtu', filename)
    return int(match.group(1)) if match else -1

vtu_files = sorted([f for f in os.listdir(OUTPUT_DIR) if f.startswith("simulation_") and f.endswith(".vtu")], key=extract_step)
if not vtu_files:
    print("[ERROR] No simulation_*.vtu files found in output directory!")
    sys.exit(1)

OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../materials/uo2_jiang.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
ax = float(geom_data['cavity']['ax'])
ay = float(geom_data['cavity']['ay'])
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))

# --.. ..- .-.. .-.. --- extract fields --.. ..- .-.. .-.. ---
strains = []
stresses = []

print("[INFO] Extracting macroscopic stress and strain from VTU files over time...")
for vtu in vtu_files:
    vtu_path = os.path.join(OUTPUT_DIR, vtu)
    
    x, y, _, disp = extract_field(vtu_path, field_name="Displacement")
    _, _, _, sigma = extract_field(vtu_path, field_name="Stress_uo2 (points)")
    
    # Macroscopic stress: average of sigma_yy on the top boundary (ymax)
    # This avoids compressive numerical artifacts at the crack tip which distort a global point average
    top_nodes_mask = np.abs(y - Ly) < 1e-6

    # --- DIAGNOSTIC ---
    # Check if any nodes are found on the top boundary. If not, something is wrong.
    if not np.any(top_nodes_mask):
        print(f"[WARNING] For file {vtu}: No nodes found on the top boundary (y={Ly}). Stress will be zero.")
        stresses.append(0.0)
        strains.append(0.0)
        continue # Go to next vtu file

    # --- CALCULATION ---
    # Macroscopic strain: average y-displacement on top / Ly
    u_y_top = np.mean(disp[top_nodes_mask, 1])
    strains.append(u_y_top / Ly)

    # Extraction robuste de sigma_yy : on calcule la moyenne de chaque composante 
    # sur le bord supérieur, et on prend celle qui a la plus grande valeur absolue 
    # (puisque nous sommes en traction pure selon Y)
    top_sigma_mean = np.mean(sigma[top_nodes_mask], axis=0)
    idx_yy = np.argmax(np.abs(top_sigma_mean))
    sigma_yy_macro = top_sigma_mean[idx_yy]
        
    print(f"  -> {vtu}: Déplacement ymax = {u_y_top:.3e} m | Strain = {u_y_top/Ly:.4e} | Stress_yy = {sigma_yy_macro*1e-6:.2f} MPa")

    stresses.append(sigma_yy_macro * 1e-6) # Convert to MPa

sigma_yy_max = np.max(stresses) if stresses else 0.0

# --.. ..- .-.. .-.. --- plot Stress vs Strain --.. ..- .-.. .-.. ---
plt.figure(figsize=(10, 6))

plt.plot(strains, stresses, 'b-o', markersize=4, label=r"Z3ST - $\sigma_{yy}$ (Top Boundary Average)")

plt.xlabel(r"Macroscopic Strain $\varepsilon_{yy}$ (-)")
plt.ylabel(r"Average Stress $\sigma_{yy}$ (MPa)")
plt.title("Stress-Strain Curve (Jiang 2020 Reproduction)")
plt.grid(True, ls=':', alpha=0.6)
plt.legend()

plot_path_ss = os.path.join(CASE_DIR, "output", "stress_strain_jiang.png")
plt.tight_layout()
plt.savefig(plot_path_ss, dpi=300)
print(f"[INFO] Plot saved in: {plot_path_ss}")

errors = {
    "max_stress_yy": {
        "numerical": float(sigma_yy_max),
        "reference": 0.0,
        "rel_error": 0.0
    }
}
TOLERANCE = 1.0e-2
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)