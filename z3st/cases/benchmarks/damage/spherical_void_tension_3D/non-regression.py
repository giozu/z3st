#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: spherical_void_tension_3D

RVE tensile test on a 120 um UO2 cube containing a single unpressurised
spherical void (R = 7.9 um), after Jiang et al. (2020). A Dirichlet_z ramp on
zmax drives the cube; the crack initiates at the void equator, so that plane
acts as the grain boundary.

Reported metrics -- there is no closed-form reference for this configuration,
so `reference` is left at 0.0 and the GOLD regression check is what actually
guards the case. The assertions below catch a degenerate run (no output, no
load-up, no softening) before it can be blessed as a gold.
"""

import os
import glob
import sys

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_verification import pass_fail_check, regression_check

CASE_DIR = os.path.dirname(__file__)
OUT_DIR = os.path.join(CASE_DIR, "output")

Lz = 120e-6  # RVE height

vtu_files = sorted(glob.glob(os.path.join(OUT_DIR, "fields_*.vtu")))
if not vtu_files:
    print("No VTU files found in output/")
    exit()

macroscopic_strain = []
macroscopic_stress = []

print("Extracting macroscopic stress-strain curve...")
for vtu in vtu_files:
    # Extract average cell stress (Sigma_zz)
    _, _, _, S_all = extract_field(vtu, field_name="Stress (cells)")
    s_zz = S_all[:, 8] # Index 8 corresponds to the ZZ component
    Sigma_zz = np.mean(s_zz)
    
    # Extract applied displacement (E_zz)
    mesh = pv.read(vtu)
    U = mesh.point_data.get("Displacement", np.zeros((mesh.n_points, 3)))
    U_z_max = np.max(U[:, 2]) # Maximum displacement at the top
    E_zz = U_z_max / Lz
    
    macroscopic_strain.append(E_zz)
    macroscopic_stress.append(Sigma_zz)

plt.figure(figsize=(6, 4))
plt.plot(macroscopic_strain, np.array(macroscopic_stress) / 1e6, "r-o", lw=1.5, markersize=4)
plt.xlabel("Macroscopic Strain $E_{zz}$ (-)")
plt.ylabel(r"Macroscopic Stress $\Sigma_{zz}$ (MPa)")
plt.title("RVE Macroscopic Stress-Strain Curve")
plt.grid(True, which="both", ls=":", lw=0.5)
plt.tight_layout(rect=[0, 0, 1, 1])
plot_path = os.path.join(OUT_DIR, "macroscopic_stress_strain.png")
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
print(f"Curve saved: {plot_path}")

# Save data to a CSV file
csv_path = os.path.join(OUT_DIR, "macroscopic_stress_strain.csv")
data_to_save = np.column_stack((macroscopic_strain, np.array(macroscopic_stress) / 1e6))
np.savetxt(csv_path, data_to_save, delimiter=",", header="Strain,Stress_MPa", comments="")
print(f"Raw CSV data saved: {csv_path}")

# --.. ..- .-.. .-.. --- metrics --.. ..- .-.. .-.. ---
strain = np.asarray(macroscopic_strain, dtype=float)
stress_mpa = np.asarray(macroscopic_stress, dtype=float) / 1e6

i_peak = int(np.argmax(stress_mpa))
sigma_peak = float(stress_mpa[i_peak])
strain_peak = float(strain[i_peak])
sigma_final = float(stress_mpa[-1])
softening_ratio = sigma_final / sigma_peak if sigma_peak > 0 else np.nan

print("\n" + "=" * 70)
print(f"  RVE fracture strength : {sigma_peak:8.2f} MPa at E_zz = {strain_peak:.4e}")
print(f"  final / peak stress   : {softening_ratio:8.4f}")
print("=" * 70)

# --.. ..- .-.. .-.. --- degenerate-run guards --.. ..- .-.. .-.. ---
if len(stress_mpa) < 3:
    sys.exit(f"[ERROR] only {len(stress_mpa)} steps recovered; nothing to characterise")
if sigma_peak <= 0.0:
    sys.exit(f"[ERROR] peak macroscopic stress is {sigma_peak:.3e} MPa; the RVE never loaded up")
if i_peak == len(stress_mpa) - 1:
    sys.exit("[ERROR] stress peaks at the last step: the ramp ends before fracture, "
             "so no softening branch was captured -- extend the ramp before blessing")

errors = {
    "sigma_zz_peak_MPa":  {"numerical": sigma_peak,      "reference": 0.0, "rel_error": 0.0},
    "strain_at_peak":     {"numerical": strain_peak,     "reference": 0.0, "rel_error": 0.0},
    "softening_ratio":    {"numerical": softening_ratio, "reference": 0.0, "rel_error": 0.0},
}

OUT_JSON = os.path.join(OUT_DIR, "non-regression.json")
TOLERANCE = 1.0e-2
pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)
