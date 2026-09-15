#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: notched_plate_2D

non-regression script
---------------------

V-notch at mid-span of the top edge, loaded by a prescribed displacement ramp
u_x on the right edge (left Clamp_x, bottom Clamp_y). The ramp passes the peak
load, so the reaction force F_x(u_x) on the right edge traces the softening
branch while the crack runs from the notch tip towards the bottom edge.

Checks (per-step VTUs and energies.txt):
  - crack_path_offset  the crack grows on x = Lx/2: max |x - Lx/2| of the D > 0.9
                       nodes below the notch tip, over Lx
  - peak_force, u_at_peak, force_ratio_final, crack_depth_final, E_frac_final
                       no closed form, tracked by the gold
"""

import os, yaml, re
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from z3st.utils.non_regression import edge_force, error_metric, finish, tracked
from z3st.utils.utils_extract_vtu import extract_field

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../../materials/high_carbon_steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")

with open(INPUT_FILE, 'r') as f:
    lc = float(yaml.safe_load(f)["damage"]["lc"])

with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))
Gc = float(mat_data.get('Gc'))
sigma_c = ((27 * E * Gc) / (256 * lc))**0.5

with open(MESH_GEO_FILE, 'r') as f:
    D_notch = float(re.search(r'D_notch\s*=\s*([\d\.]+);', f.read()).group(1))
Y_tip = Ly - D_notch

with open(BC_FILE, 'r') as f:
    bc_data = yaml.safe_load(f)
u_ramp = next(bc["displacement"] for bc in bc_data["mechanical"]["steel"] if bc["type"] == "Dirichlet_x")

print(f"[INFO] Geometry: Lx = {Lx} m, Ly = {Ly} m, notch tip at y = {Y_tip} m")
print(f"[INFO] Material: E = {E:.2e} Pa, Gc = {Gc} J/m2, sigma_c (AT2) = {sigma_c:.3e} Pa")

TOLERANCE = 1e-2            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
u_steps = np.array(u_ramp[:len(VTU_FILES)], dtype=float)
F_steps = np.array([edge_force(v, "x", Lx, Lx / 2000) for v in VTU_FILES])
i_peak = int(np.argmax(F_steps))
F_peak, u_peak = F_steps[i_peak], u_steps[i_peak]
print(f"[INFO] Peak reaction {F_peak*1e-6:.3f} MN/m at u_x = {u_peak*1e6:.1f} um (step {i_peak}), "
      f"final {F_steps[-1]*1e-6:.3f} MN/m ({F_steps[-1]/F_peak:.3f} of peak)")

energies = np.genfromtxt(os.path.join(CASE_DIR, "energies.txt"), names=True)
E_frac = np.atleast_1d(energies["E_frac"])

x_d, y_d, _, D_final = extract_field(VTU_FILES[-1], field_name="Damage")
below = (D_final > 0.9) & (y_d < Y_tip - 2 * lc)
crack_offset = float(np.max(np.abs(x_d[below] - Lx / 2))) if np.any(below) else 0.0
crack_depth = float(Ly - np.min(y_d[D_final > 0.9])) if np.any(D_final > 0.9) else 0.0
print(f"[INFO] Crack depth from the top edge: {crack_depth*1e3:.1f} mm (notch {D_notch*1e3:.0f} mm), "
      f"max offset from x = Lx/2: {crack_offset*1e3:.2f} mm")

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
plt.figure(figsize=(7, 5))
plt.plot(u_steps * 1e6, F_steps * 1e-6, "-o", color="#0072B2", lw=2, markersize=4)
plt.plot(u_peak * 1e6, F_peak * 1e-6, "s", color="#D55E00", label="Peak")
plt.xlabel(r"Prescribed displacement $u_x$ ($\mu$m)")
plt.ylabel(r"Reaction force $F_x$ (MN/m)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "force_displacement.png"))

plt.figure(figsize=(8, 5))
plt.plot(energies['Step'], energies['E_el'], "-o", color="#0072B2", label='Elastic Energy ($E_{el}$)')
plt.plot(energies['Step'], energies['E_frac'], "-s", color="#D55E00", label='Fracture Energy ($E_{frac}$)')
plt.plot(energies['Step'], energies['E_tot'], 'k--', label='Total Energy ($E_{tot}$)')
plt.xlabel('Step')
plt.ylabel('Energy (J)')
plt.grid(True, ls=':', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "energy_balance.png"))

plt.figure(figsize=(10, 5))
sc = plt.scatter(x_d, y_d, c=D_final, cmap='jet', s=1)
plt.colorbar(sc, label="Damage $D$")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title(f"Final damage, u_x = {u_steps[-1]*1e6:.0f} um")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "damage_final.png"), dpi=200)
print("[INFO] plots saved in output/")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "crack_path_offset": error_metric(crack_offset, crack_offset / Lx),
    "peak_force": tracked(F_peak),
    "u_at_peak": tracked(u_peak),
    "force_ratio_final": tracked(F_steps[-1] / F_peak),
    "crack_depth_final": tracked(crack_depth),
    "E_frac_final": tracked(E_frac[-1]),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
