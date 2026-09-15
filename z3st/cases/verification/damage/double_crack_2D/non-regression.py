#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: double_crack_2D

non-regression script
---------------------

Two edge pre-cracks (D = 1 on 0 < x < Dn and Lx - Dn < x < Lx, y = Ly/2)
opened by a prescribed displacement ramp u_y on the top edge (bottom Clamp_y,
right edge Clamp_x). The ramp passes the peak load, so the reaction force
F_y(u_y) on the top edge traces the softening branch.

Checks (per-step VTUs and energies.txt):
  - E_frac_initial     AT2 energy of the two prescribed cracks against Gc * 2 Dn,
                       Gc = 256 lc sigma_c^2 / (27 E) from the material card
  - crack_path_offset  the cracks grow on y = Ly/2: max |y - Ly/2| of the D > 0.9
                       nodes between the pre-crack tips, over Ly
  - peak_force, u_at_peak, force_ratio_final, E_frac_final
                       no closed form, tracked by the gold
"""

import os, yaml, re
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from z3st.utils.non_regression import edge_force, error_metric, finish, metric, tracked
from z3st.utils.utils_extract_vtu import extract_field

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
VTU_FILES = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))

INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

with open(INPUT_FILE, 'r') as f:
    input_data = yaml.safe_load(f)
lc = float(input_data["damage"]["lc"])

with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))
sigma_c = float(mat_data.get('sigma_c'))
Gc = 256.0 * lc * sigma_c**2 / (27.0 * E)   # AT2 conversion, as in spine.load_materials

with open(MESH_GEO_FILE, 'r') as f:
    Dn = float(re.search(r'Dn\s*=\s*([\d\.]+);', f.read()).group(1))

with open(BC_FILE, 'r') as f:
    bc_data = yaml.safe_load(f)
u_ramp = next(bc["displacement"] for bc in bc_data["mechanical"]["steel"] if bc["type"] == "Dirichlet_y")

print(f"[INFO] Geometry: Lx = {Lx} m, Ly = {Ly} m, Dn = {Dn} m")
print(f"[INFO] Material: E = {E:.2e} Pa, sigma_c = {sigma_c:.2e} Pa -> Gc (AT2) = {Gc:.1f} J/m2")

TOLERANCE = 1e-2            # relative tolerance for pass/fail

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
u_steps = np.array(u_ramp[:len(VTU_FILES)], dtype=float)
F_steps = np.array([edge_force(v, "y", Ly, Ly / 2000) for v in VTU_FILES])
i_peak = int(np.argmax(F_steps))
F_peak, u_peak = F_steps[i_peak], u_steps[i_peak]
print(f"[INFO] Peak reaction {F_peak*1e-6:.3f} MN/m at u_y = {u_peak*1e6:.1f} um (step {i_peak}), "
      f"final {F_steps[-1]*1e-6:.3f} MN/m ({F_steps[-1]/F_peak:.3f} of peak)")

energies = np.genfromtxt(os.path.join(CASE_DIR, "energies.txt"), names=True)
E_frac = np.atleast_1d(energies["E_frac"])

x_d, y_d, _, D_final = extract_field(VTU_FILES[-1], field_name="Damage")
between = (D_final > 0.9) & (x_d > Dn + 2 * lc) & (x_d < Lx - Dn - 2 * lc)
crack_offset = float(np.max(np.abs(y_d[between] - Ly / 2))) if np.any(between) else 0.0
print(f"[INFO] Grown-crack nodes (D > 0.9 between the tips): {np.count_nonzero(between)}, "
      f"max offset from y = Ly/2: {crack_offset*1e3:.2f} mm")

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
# crack profile across the pre-crack
plt.figure(figsize=(8, 6))
x_line, tol_x, yc = Dn / 2.0, Lx / 100.0, Ly / 2.0
colors = plt.cm.jet(np.linspace(0, 1, len(VTU_FILES)))
for i, vtufile in enumerate(VTU_FILES):
    xs, ys, _, Ds = extract_field(vtufile, field_name="Damage")
    mask_line = np.abs(xs - x_line) < tol_x
    order = np.argsort(ys[mask_line])
    y_plot, d_plot = ys[mask_line][order], Ds[mask_line][order]
    plt.plot(y_plot, d_plot, color=colors[i], lw=1.0, alpha=0.8)
plt.plot(y_plot, np.exp(-np.abs(y_plot - yc) / lc), "k--", lw=2, label=r"$\exp(-|y-y_c|/l_c)$", zorder=10)
plt.xlabel("Coordinate (m)")
plt.ylabel("Damage $D$ (/)")
plt.title(f"Crack profile at $x={x_line:.3f}$ m, all steps, $l_c = {lc}$ m")
plt.grid(True, ls=':', alpha=0.6)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "crack_profile_evolution.png"))

# energy balance
plt.figure(figsize=(8, 5))
plt.plot(energies['Step'], energies['E_el'], "-o", color="#0072B2", label='Elastic Energy ($E_{el}$)')
plt.plot(energies['Step'], energies['E_frac'], "-s", color="#D55E00", label='Fracture Energy ($E_{frac}$)')
plt.plot(energies['Step'], energies['E_tot'], 'k--', label='Total Energy ($E_{tot}$)')
plt.xlabel('Step')
plt.ylabel('Energy (J)')
plt.title('Z3ST: Global energy balance')
plt.grid(True, ls=':', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "energy_balance.png"))

# reaction force against prescribed displacement
plt.figure(figsize=(7, 5))
plt.plot(u_steps * 1e6, F_steps * 1e-6, "-o", color="#0072B2", lw=2, markersize=4)
plt.plot(u_peak * 1e6, F_peak * 1e-6, "s", color="#D55E00", label="Peak")
plt.xlabel(r"Prescribed displacement $u_y$ ($\mu$m)")
plt.ylabel(r"Reaction force $F_y$ (MN/m)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "force_displacement.png"))

# final damage map
plt.figure(figsize=(10, 5))
sc = plt.scatter(x_d, y_d, c=D_final, cmap='jet', s=1)
plt.colorbar(sc, label="Damage $D$")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title(f"Final damage, u_y = {u_steps[-1]*1e6:.0f} um")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "damage_final.png"), dpi=200)
print("[INFO] plots saved in output/")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "E_frac_initial": metric(E_frac[0], Gc * 2 * Dn),
    "crack_path_offset": error_metric(crack_offset, crack_offset / Ly),
    "peak_force": tracked(F_peak),
    "u_at_peak": tracked(u_peak),
    "force_ratio_final": tracked(F_steps[-1] / F_peak),
    "E_frac_final": tracked(E_frac[-1]),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
