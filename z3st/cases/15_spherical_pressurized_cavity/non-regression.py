#!/usr/bin/env python3
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from z3st.utils.utils_extract_vtu import *

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
    _, _, _, S_all = extract_field(vtu, field_name="Stress_solid (cells)")
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
