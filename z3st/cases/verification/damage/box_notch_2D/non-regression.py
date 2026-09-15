#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: box_knotch_2D

non-regression script
---------------------

V-notch at mid-span of the top edge, loaded by a coarse prescribed
displacement ramp u_x on the right edge (left Clamp_x, bottom Clamp_y), past
the peak load. notched_plate_2D runs the same specimen on a finer ramp and
checks the force-displacement curve; this case checks the final fields:
  - max_damage          a fully formed crack, D = 1
  - max_stress_xx       the largest sigma_xx over the whole history (tension
                        across the notch) against the AT2 strength
                        sigma_c = sqrt(27 E Gc / (256 lc))
  - crack_path_offset   the crack grows on x = Lx/2 (max |x - Lx/2| of the
                        D > 0.9 nodes below the notch tip, over Lx)
  - peak_force          tracked by the gold
"""

import os, re
from glob import glob
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.non_regression import case_paths, edge_force, error_metric, finish, tracked
from z3st.utils.utils_extract_vtu import extract_field, list_fields

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR, _, OUT_JSON = case_paths(__file__)
VTU_FILES = sorted(glob(os.path.join(CASE_DIR, "output", "fields_*.vtu")))
VTU_FILE = VTU_FILES[-1]   # final state of the displacement ramp
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../../materials/high_carbon_steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")

# Phase-field / damage
with open(INPUT_FILE, 'r') as f: 
    input_data = yaml.safe_load(f)
dmg_cfg = input_data.get("damage", {})
lc = float(dmg_cfg["lc"])

# Geometry
with open(GEOMETRY_FILE, 'r') as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get('Lx'))
Ly = float(geom_data.get('Ly'))

with open(MESH_GEO_FILE, 'r') as f:
    content = f.read()

W_notch = float(re.search(rf'W_notch\s*=\s*([\d\.]+);', content).group(1))
D_notch = float(re.search(rf'D_notch\s*=\s*([\d\.]+);', content).group(1))

Y_tip   = Ly - D_notch
x_target = Lx / 2

# Material
with open(MATERIAL_FILE, 'r') as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get('E'))
nu = float(mat_data.get('nu'))
Gc = float(mat_data.get('Gc'))
# sigma_c = float(mat_data.get('sigma_c'))
sigma_c = ((27 * E * Gc) / (256 * lc))**0.5

TOLERANCE = 5e-1 

# --.. ..- .-.. .-.. --- Data --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# Damage
x_d, y_d, _, D_all = extract_field(VTU_FILE, field_name="Damage")
d_max = np.max(D_all)

# Stress
x_s, y_s, _, S_all = extract_field(VTU_FILE, field_name="Stress (cells)")
sigma_xx_max = max(np.max(extract_field(v, field_name="Stress (cells)")[3][:, 0]) for v in VTU_FILES)

# Reaction on the loaded edge and crack path
F_steps = np.array([edge_force(v, "x", Lx, Lx / 2000) for v in VTU_FILES])
F_peak = float(np.max(F_steps))
print(f"[INFO] Peak reaction {F_peak*1e-6:.3f} MN/m (step {int(np.argmax(F_steps))}), "
      f"final {F_steps[-1]*1e-6:.3f} MN/m")
below = (D_all > 0.9) & (y_d < Y_tip - 2 * lc)
crack_offset = float(np.max(np.abs(x_d[below] - x_target))) if np.any(below) else 0.0

# Profiles
mask_d = np.abs(x_d - x_target) < (Lx/200)
idx_d = np.argsort(y_d[mask_d])
y_prof = y_d[mask_d][idx_d]
D_prof = D_all[mask_d][idx_d]

mask_s = np.abs(x_s - x_target) < (Lx/500)
idx_s = np.argsort(y_s[mask_s])
y_s_prof = y_s[mask_s][idx_s]
sigma_xx_prof = S_all[mask_s, 0][idx_s]

y_slice = Y_tip - 2*lc
mask_x = np.abs(y_d - y_slice) < (Ly/100)
idx_x = np.argsort(x_d[mask_x])
x_slice_prof = x_d[mask_x][idx_x]
d_slice_prof = D_all[mask_x][idx_x]

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---

# PLOT 1:
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(x_slice_prof, d_slice_prof, "-o", color="#D55E00", markersize=4, label=f"Damage at y={y_slice:.2f}")
plt.axhline(1.0, color='k', linestyle='--', alpha=0.3)
plt.xlabel("x (m)")
plt.ylabel("Damage $d$")
plt.grid(True, ls=':')
plt.legend()

plt.subplot(1, 2, 2)
sc = plt.scatter(x_d, y_d, c=D_all, cmap='jet', s=2)
plt.colorbar(sc, label="Damage $d$")
plt.axhline(Y_tip, color='white', linestyle=':', alpha=0.5, label="Notch tip")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title(f"Max Damage: {d_max:.3f}")
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "damage_check.png"), dpi=300)

# PLOT 2:
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Damage
ax1.plot(y_prof, D_prof, "-", color="#D55E00", lw=2, label="Damage $d$")
ax1.fill_between(y_prof, D_prof, color="#D55E00", alpha=0.1)
ax1.axvline(Y_tip, color='k', ls=':', label="Notch tip")
ax1.set_ylabel("Damage $d$ ", color="#D55E00")
ax1.set_ylim(-0.05, 1.05)
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend()
ax1.set_title(rf"Z3ST analysis: centerline profile (x = {x_target:.2f} m)\n"
              rf"Steel: $\sigma_c$={sigma_c*1e-6:.0f} MPa, $G_c$={Gc:.1f} J/m²")

# Stress
ax2.plot(y_s_prof, sigma_xx_prof * 1e-6, "-", color="#0072B2", lw=2, label=r"$\sigma_{xx}$")
ax2.axhline(sigma_c * 1e-6, color='black', ls='--', lw=1, label=rf"$\sigma_c$ Limit")
ax2.axvline(Y_tip, color='k', ls=':', label="Notch tip")
ax2.set_xlabel("Vertical position $y$ (m)")
ax2.set_ylabel(r"Normal stress $\sigma_{xx}$ (MPa)", color="#0072B2")
ax2.grid(True, ls=":", alpha=0.6)
ax2.legend(loc='best')

plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "damage_stress_profile.png"), dpi=300)
print(f"[INFO] Detailed profiles saved in: output/damage_stress_profile.png")

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
        "rel_error": float(abs(sigma_xx_max - sigma_c) / sigma_c)
    },
    "crack_path_offset": error_metric(crack_offset, crack_offset / Lx),
    "peak_force": tracked(F_peak),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)

print(f"\n[INFO] Max damage: {d_max:.4f}")
print(f"[INFO] Max sigma_xx over the history: {sigma_xx_max*1e-6:.2f} MPa")

print("[INFO] non-regression completed.\n")