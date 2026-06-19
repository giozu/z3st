#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: two_elliptical_cavities_2D

non-regression script
---------------------
"""

import os, re
import yaml
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_extract_xdmf import *
from z3st.utils.utils_verification import *

from z3st.materials.oxide import Gc_numpy as Gc

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
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../materials/oxide.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
INPUT_FILE = os.path.join(CASE_DIR, "input.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")

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

# Applied bubble pressure = peak of the Neumann ramp on the cavity, read from
# boundary_conditions.yaml (the ramp ends at 15 MPa).
with open(BC_FILE, 'r') as f:
    bc_data = yaml.safe_load(f)
p_applied = max(
    abs(float(v))
    for e in bc_data["mechanical"]["solid"]
    if e.get("type") == "Neumann" and e.get("region") == "cavity"
    for v in e["traction"]
)
print(f"[INFO] Applied bubble pressure (peak of ramp): {p_applied:.1f} MPa")

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

# --.. ..- .-.. .-.. --- extract fields --.. ..- .-.. .-.. ---
list_fields_xdmf(XDMF_FILE)

X_tip = ax
y_target = 0.0

x_d, y_d, _, D_all = extract_field_xdmf(XDMF_FILE, field_name="Damage")
d_max = np.max(D_all)

x, y, _, sigma = extract_field_xdmf(XDMF_FILE, field_name="Stress")
sigma_yy_max = np.max(sigma[:, 4]) 

# mask = (np.abs(y - y_target) < (Ly/500)) & (x >= X_tip)
mask = (np.abs(y - y_target) < (Ly/100))
idx_line = np.argsort(x[mask])
x_line = x[mask][idx_line]
sigma_yy_line = sigma[mask, 4][idx_line]

# --.. ..- .-.. .-.. --- sigma_yy along the GB ligament (y≈0) --.. ..- .-..
edges = np.linspace(x_line.min(), x_line.max(), 160)
xc = 0.5 * (edges[:-1] + edges[1:])
syy_mean = np.array([
    sigma_yy_line[(x_line >= a) & (x_line < b)].mean()
    if np.any((x_line >= a) & (x_line < b)) else np.nan
    for a, b in zip(edges[:-1], edges[1:])
])

plt.figure(figsize=(10, 5))
for xb in (-Lx / 4, +Lx / 4):                        # cavity footprints + tips
    plt.axvspan(xb - ax, xb + ax, color="0.92", zorder=0)
    for s in (-1, +1):
        plt.axvline(xb + s * ax, color="r", ls="--", lw=0.8, alpha=0.6)
plt.scatter(x_line, sigma_yy_line, s=4, c="0.7", alpha=0.4, label="band samples")
plt.plot(xc, syy_mean, "b-", lw=2, label=r"$\sigma_{yy}$ (mean along $y\approx0$)")
plt.xlabel(r"x along grain boundary ($\mu$m)")
plt.ylabel(r"$\sigma_{yy}$ (MPa)")
plt.title("GB stress profile — bubble tips (red dashed) concentrate the stress")
plt.grid(True, ls=":", alpha=0.5)
plt.legend(loc="upper right")
plot_path = os.path.join(CASE_DIR, "output", "stress_profile_tip.png")
plt.tight_layout()
plt.savefig(plot_path, dpi=200)
print(f"[INFO] Plot saved in: {plot_path}")

# Gc profile
y_line = np.linspace(-Ly/2, Ly/2, 500) # (micron)
gc_line = Gc(y_line)

sigma_c = ((27 * E * gc_line) / (256 * lc))**0.5

plt.figure(figsize=(8, 5))
# gc_line is in pJ/µm², which equals J/m² numerically in this unit system.
plt.plot(y_line * 1000, gc_line, 'g-', label="$G_c$ profile (GB zone)")
plt.xlabel(r"Distance $y$ (nm)")
plt.ylabel("$G_c$ ($J/m^2$)")
plt.grid(True, ls=':')
plt.legend()
plt.savefig(os.path.join(CASE_DIR, "output", "gc_profile_check.png"))

# --.. ..- .-.. .-.. --- field overview: mesh + damage + stress + Gc --.. ..- .-..
try:
    with h5py.File(XDMF_FILE.replace(".xdmf", ".h5"), "r") as f:
        topo = np.array(f["Mesh/mesh/topology"])            # cell connectivity
    triang = Triangulation(x_d, y_d, topo)                  # nodal coords + topology
    gc_nodal = Gc(y_d)                                      # shared toughness on nodes

    fig, ax = plt.subplots(2, 2, figsize=(12, 11))
    ax[0, 0].triplot(triang, lw=0.15, color="0.4")
    ax[0, 0].set_title("Mesh (2 elliptical cavities)")

    cf = ax[0, 1].tricontourf(triang, D_all, levels=20, cmap="inferno")
    fig.colorbar(cf, ax=ax[0, 1]); ax[0, 1].set_title(f"Damage  (max = {d_max:.3g})")

    pc = ax[1, 0].tripcolor(triang, facecolors=sigma[:, 4], cmap="viridis")
    fig.colorbar(pc, ax=ax[1, 0]); ax[1, 0].set_title(r"$\sigma_{yy}$ (MPa)")

    cg = ax[1, 1].tricontourf(triang, gc_nodal, levels=20, cmap="cividis")
    fig.colorbar(cg, ax=ax[1, 1]); ax[1, 1].set_title(r"$G_c$ (pJ/$\mu$m$^2$ = J/m$^2$)")

    for a in ax.flat:
        a.set_aspect("equal"); a.set_xlabel(r"x ($\mu$m)"); a.set_ylabel(r"y ($\mu$m)")
    fig.tight_layout()
    overview_path = os.path.join(OUTPUT_DIR, "fields_overview.png")
    fig.savefig(overview_path, dpi=200)
    print(f"[INFO] Field overview saved in: {overview_path}")
except Exception as e:
    print(f"[WARNING] field overview plot skipped: {e}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
# The config is the equilibrium baseline (15 MPa, no fracture), so the analytic
# references describe that regime:
#   - max_damage: baseline must not crack -> d_max ~ 0 (percolation is 0.9).
#   - max_stress_yy: bubble-tip stress concentration. The isolated pressurised-
#     ellipse factor K_t underpredicts ~2x because the two bubbles share the load
#     through the intact ligament; the ligament-corrected estimate
#     K_t/(1-Fc)*p_applied is accurate to ~20% (the 2-bubble interaction and
#     finite domain have no closed form). Drift is guarded by regression_check
#     vs gold (rtol 1e-3), not this estimate.
sigma_yy_ref = intensification_factor / (1.0 - Fc_area) * p_applied

# Loose analytic tolerance: the tip-concentration estimate is approximate. A
# pass means the right regime (no fracture, expected concentration); drift vs
# code/mesh is guarded by regression_check vs gold.
TOLERANCE = 0.25

errors = {
    "max_damage": {
        "numerical": float(d_max),
        "reference": 0.0,                  # equilibrium baseline: no fracture
        "rel_error": float(abs(d_max))     # absolute deviation from 0
    },
    "max_stress_yy": {
        "numerical": float(sigma_yy_max),
        "reference": float(sigma_yy_ref),
        "rel_error": float(abs(sigma_yy_max - sigma_yy_ref)/sigma_yy_ref)
    }
}

print(f"\n[RESULTS]")
print(f"  → Applied bubble pressure: {p_applied:.1f} MPa")
print(f"  → Max Damage: {d_max:.3e}  (reference 0 — no fracture at equilibrium)")
print(f"  → Max Stress_yy: {sigma_yy_max:.2f} MPa "
      f"(analytic tip estimate K_t/(1-Fc)*p = {sigma_yy_ref:.2f} MPa)")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)
