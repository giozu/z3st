#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: elliptical_cavity_tension_2D

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
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../../materials/uo2_jiang.yaml")
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
max_damage = []          # per step, for the initiation metric

DAMAGE_CRITICAL_THRESHOLD = 0.5

print("[INFO] Extracting macroscopic stress and strain from VTU files over time...")
for vtu in vtu_files:
    vtu_path = os.path.join(OUTPUT_DIR, vtu)
    
    x, y, _, disp = extract_field(vtu_path, field_name="Displacement")
    _, _, _, sigma = extract_field(vtu_path, field_name="Stress (points)")
    _, _, _, D = extract_field(vtu_path, field_name="Damage")
    
    # Macroscopic stress: average of sigma_yy on the top boundary (ymax)
    # This avoids compressive numerical artifacts at the crack tip which distort a global point average
    top_nodes_mask = np.abs(y - Ly) < 1e-6

    # --- DIAGNOSTIC ---
    # Check if any nodes are found on the top boundary. If not, something is wrong.
    if not np.any(top_nodes_mask):
        print(f"[WARNING] For file {vtu}: No nodes found on the top boundary (y={Ly}). Stress will be zero.")
        stresses.append(0.0)
        strains.append(0.0)
        max_damage.append(float(np.max(D)))
        continue # Go to next vtu file

    # --- CALCULATION ---
    # Macroscopic strain: average y-displacement on top / Ly
    u_y_top = np.mean(disp[top_nodes_mask, 1])
    strains.append(u_y_top / Ly)

    # Extraction robuste de sigma_yy : on calcule la moyenne de chaque composante 
    # sur le bord supérieur, et on prend celle qui a la plus grande valeur absolue 
    # (puisque nous sommes en traction pure selon Y)
    # Stress is written as a flattened 3x3 tensor (9 components, row-major),
    # so sigma_yy is component 4. The previous largest-magnitude heuristic
    # would silently return sigma_xx wherever it happened to dominate.
    top_sigma_mean = np.mean(sigma[top_nodes_mask], axis=0)
    sigma_yy_macro = top_sigma_mean[4]
        
    print(f"  -> {vtu}: Déplacement ymax = {u_y_top:.3e} m | Strain = {u_y_top/Ly:.4e} | Stress_yy = {sigma_yy_macro*1e-6:.2f} MPa")

    stresses.append(sigma_yy_macro * 1e-6) # Convert to MPa
    max_damage.append(float(np.max(D)))

sigma_yy_max = np.max(stresses) if stresses else 0.0

# --.. ..- .-.. .-.. --- Display Max Rupture Stress --.. ..- .-.. .-.. ---
print("\n" + "="*70)
print(f"{'█'*70}")
print(f"  📊 MAXIMUM RUPTURE STRESS: {sigma_yy_max:.2f} MPa")
print(f"{'█'*70}")
print("="*70 + "\n")

# --.. ..- .-.. .-.. --- degenerate-run guards --.. ..- .-.. .-.. ---
# The reference below is a placeholder (0.0), so pass_fail_check cannot fail:
# without these guards a run that produced no load at all still reports PASS.
if len(stresses) < 3:
    sys.exit(f"[ERROR] only {len(stresses)} steps recovered; nothing to characterise")
if sigma_yy_max <= 1.0:
    sys.exit(f"[ERROR] peak macroscopic sigma_yy is {sigma_yy_max:.3e} MPa -- the specimen "
             "never loaded up. Check that a remote load is actually applied on ymax; a "
             "cavity-pressure-driven case leaves the top surface traction-free and must "
             "use the critical-cracking-pressure post-processing instead.")
if int(np.argmax(stresses)) == len(stresses) - 1:
    sys.exit("[ERROR] stress peaks at the last step: the ramp ends before fracture, so no "
             "softening branch was captured -- extend the ramp before blessing")

# --.. ..- .-.. .-.. --- headline metric: initiation, not peak --.. ..- .-.. .-.. ---
# The quantity this case guards is the REMOTE REACTION STRESS AT CRACK
# INITIATION -- the macroscopic sigma_yy at the first step where max(Damage)
# crosses the threshold. That is what pairs with the pressurised sibling's
# critical cracking pressure. The peak macroscopic stress over the whole ramp
# is a different (post-initiation) quantity and is reported for information
# only; see CALIBRATION.md.
max_damage = np.asarray(max_damage, dtype=float)
above = np.where(max_damage >= DAMAGE_CRITICAL_THRESHOLD)[0]

if above.size == 0:
    sys.exit(
        f"[ERROR] max(Damage) never reached {DAMAGE_CRITICAL_THRESHOLD} over the ramp "
        f"(peaked at {max_damage.max():.4f}): the crack never initiated, so there is no "
        f"initiation stress to report. Extend the ramp before blessing."
    )

i_init = int(above[0])
initiation_stress = float(stresses[i_init])
initiation_strain = float(strains[i_init])

print("\n" + "=" * 70)
print(f"  CRACK INITIATION at step {i_init}: remote sigma_yy = {initiation_stress:.2f} MPa "
      f"(strain {initiation_strain:.4e}, max D = {max_damage[i_init]:.4f})")
print(f"  peak macroscopic sigma_yy over the ramp: {sigma_yy_max:.2f} MPa  [informational]")
print("=" * 70 + "\n")

errors = {
    "initiation_stress_MPa": {
        "numerical": initiation_stress,
        "reference": 0.0,
        "rel_error": 0.0
    },
    "initiation_strain": {
        "numerical": initiation_strain,
        "reference": 0.0,
        "rel_error": 0.0
    }
}
TOLERANCE = 1.0e-2
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

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
