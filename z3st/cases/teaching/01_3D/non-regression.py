#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: teaching/01_3D  --  1D bar in uniaxial tension on a 3D hex mesh.

Reference: Bower, "Applied Mechanics of Solids" (CRC, 2010),
ch. 7.2 + ch. 8.1.5.

Companion to 01_1D. The 3D mesh leaves the transverse Poisson contraction
in y and z completely free, so the FE recovers the *engineering* bar:

    sigma_xx = P,   sigma_yy = sigma_zz = 0
    sigma_xy = sigma_yz = sigma_xz = 0
    u_x(x, y, z) = (P / E) x
    u_y(x, y, z) = -nu (P / E) y
    u_z(x, y, z) = -nu (P / E) z
    u_x(L)       = P L / E

The displacement field is linear in (x, y, z), so CG1 Lagrange elements
recover it exactly (to machine precision). Mesh refinement changes nothing.

Contrast with 01_1D (same physical problem, 1D mesh):
    u_x_3D / u_x_1D = (lambda + 2G) / E
                    = (1 - nu) / [(1 + nu)(1 - 2 nu)]
                    ≈ 1.346 for steel (nu = 0.3)
The 1D-mesh result underestimates the engineering elongation by ~26 %
because it implicitly enforces a transverse-rigid sleeve around the bar.
"""

import os
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt

from z3st.utils.utils_extract_vtu import (
    list_fields,
    extract_field,
    extract_displacement,
)
from z3st.utils.utils_verification import pass_fail_check, regression_check

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")
MATERIAL_FILE = os.path.join(CASE_DIR, "../../../materials/steel.yaml")
GEOMETRY_FILE = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE = os.path.join(CASE_DIR, "mesh.geo")

with open(GEOMETRY_FILE, "r") as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data.get("Lx"))
Ly = float(geom_data.get("Ly"))
Lz = float(geom_data.get("Lz"))

with open(MATERIAL_FILE, "r") as f:
    mat_data = yaml.safe_load(f)
E = float(mat_data.get("E"))
nu = float(mat_data.get("nu"))

with open(BC_FILE, "r") as f:
    bc_data = yaml.safe_load(f)
mech_list = bc_data.get("mechanical", {}).get("steel", [])
P = next(
    (float(bc["traction"]) for bc in mech_list
     if bc.get("type") == "Neumann" and bc.get("region") == "xmax"),
    None,
)

with open(MESH_GEO_FILE, "r") as f:
    content = f.read()
nx = int(re.search(r"nx\s*=\s*(\d+);", content).group(1)) - 1
ny = int(re.search(r"ny\s*=\s*(\d+);", content).group(1)) - 1
nz = int(re.search(r"nz\s*=\s*(\d+);", content).group(1)) - 1

print(f"[INFO] Geometry  : Lx = {Lx} m, Ly = {Ly} m, Lz = {Lz} m")
print(f"[INFO] Material  : E = {E:.2e} Pa, nu = {nu}")
print(f"[INFO] Load      : P = {P:.2e} Pa ({P * 1e-6:.1f} MPa)")
print(f"[INFO] Mesh      : nx x ny x nz = {nx} x {ny} x {nz} hex elements")

TOLERANCE = 1e-4

# --.. ..- .-.. .-.. --- analytical reference (engineering bar) ----------------
sigma_xx_ref = P
sigma_yy_ref = 0.0
sigma_zz_ref = 0.0
eps_xx_ref = P / E
eps_yy_ref = -nu * P / E
u_xL_ref = eps_xx_ref * Lx

print(f"[INFO] Engineering-bar ref: eps_xx = {eps_xx_ref:.4e}, "
      f"u_x(L) = {u_xL_ref * 1e3:.4f} mm, "
      f"eps_yy = eps_zz = {eps_yy_ref:.4e}")

# --.. ..- .-.. .-.. --- inspect VTU --.. ..- .-.. .-.. ---
list_fields(VTU_FILE)

# --.. ..- .-.. .-.. --- extract FE fields --.. ..- .-.. .-.. ---
# Stress is per-cell. Sample at mid-cross-section (y ≈ Ly/2, z ≈ Lz/2) and
# sort by x to get a clean axial profile.
x_S, y_S, z_S, S_all = extract_field(VTU_FILE, field_name="Stress_steel (cells)")
y_mid = Ly / 2.0
z_mid = Lz / 2.0
mask_tol_y = Ly / (2 * max(ny, 1))
mask_tol_z = Lz / (2 * max(nz, 1))
mask = (np.abs(y_S - y_mid) < mask_tol_y) & (np.abs(z_S - z_mid) < mask_tol_z)
sort_idx = np.argsort(x_S[mask])
x_s = x_S[mask][sort_idx]
sigma_xx = S_all[mask, 0][sort_idx]
sigma_yy = S_all[mask, 4][sort_idx]
sigma_zz = S_all[mask, 8][sort_idx]
sigma_xy = S_all[mask, 1][sort_idx]

# Displacement is per-node.
x_n, y_n, z_n, u_vec = extract_displacement(VTU_FILE)
u_x = u_vec[:, 0]
u_y = u_vec[:, 1]
u_z = u_vec[:, 2]

# Tip displacement: nodes at x = Lx.
tip_mask = np.abs(x_n - Lx) < 1e-9
u_xL_num = float(np.mean(u_x[tip_mask]))

# u_x along the central axis (y ≈ Ly/2, z ≈ Lz/2).
axis_mask = (np.abs(y_n - y_mid) < mask_tol_y) & (np.abs(z_n - z_mid) < mask_tol_z)
sort_n = np.argsort(x_n[axis_mask])
x_n_axis = x_n[axis_mask][sort_n]
u_x_axis = u_x[axis_mask][sort_n]

print(f"[INFO] FE tip disp: u_x(L) = {u_xL_num * 1e3:.6f} mm "
      f"(ref {u_xL_ref * 1e3:.6f} mm)")

# Lateral contraction at the tip (any face): u_y(y=Ly) should be -nu*(P/E)*Ly
side_mask = (np.abs(x_n - Lx) < 1e-9) & (np.abs(y_n - Ly) < 1e-9)
if side_mask.any():
    u_y_top_corner = float(np.mean(u_y[side_mask]))
    u_y_ref = eps_yy_ref * Ly
    print(f"[INFO] Poisson contraction at (Lx, Ly, *): u_y = {u_y_top_corner*1e6:.3f} um "
          f"(ref {u_y_ref*1e6:.3f} um)")

# --.. ..- .-.. .-.. --- plots --.. ..- .-.. .-.. ---
Pa_to_MPa = 1e-6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(x_s, sigma_xx * Pa_to_MPa, "bo-", label=r"FE $\sigma_{xx}$",
         markersize=6, alpha=0.8)
ax1.plot(x_s, sigma_yy * Pa_to_MPa, "rs-", label=r"FE $\sigma_{yy}$",
         markersize=6, alpha=0.8)
ax1.plot(x_s, sigma_zz * Pa_to_MPa, "m^-", label=r"FE $\sigma_{zz}$",
         markersize=6, alpha=0.8)
ax1.plot(x_s, sigma_xy * Pa_to_MPa, "gv-", label=r"FE $\sigma_{xy}$",
         markersize=6, alpha=0.8)
ax1.axhline(P * Pa_to_MPa, color="k", linestyle="--", linewidth=1.2,
            label=rf"Analytic $\sigma_{{xx}} = P = {P*Pa_to_MPa:.0f}$ MPa")
ax1.axhline(0.0, color="grey", linestyle=":", linewidth=0.8,
            label=r"Analytic $\sigma_{yy}=\sigma_{zz}=\sigma_{xy}=0$")
ax1.set_xlabel("x (m)", fontsize=12)
ax1.set_ylabel("Stress (MPa)", fontsize=12)
ax1.set_title("Stress along the bar axis (y=Ly/2, z=Lz/2)", fontsize=12)
ax1.legend(loc="best", fontsize=9)
ax1.grid(True, linestyle="--", alpha=0.6)

ax2.plot(x_n_axis, u_x_axis * 1e3, "bo", label=r"FE $u_x$", markersize=7)
x_dense = np.linspace(0.0, Lx, 200)
ax2.plot(x_dense, eps_xx_ref * x_dense * 1e3, "k--", linewidth=1.5,
         label=r"Analytic $u_x = P\,x/E$")
ax2.set_xlabel("x (m)", fontsize=12)
ax2.set_ylabel(r"$u_x$ (mm)", fontsize=12)
ax2.set_title("Axial displacement along the bar axis", fontsize=12)
ax2.legend(loc="best")
ax2.grid(True, linestyle="--", alpha=0.6)

fig.suptitle(
    rf"1D bar in tension (3D hex mesh)  |  $P$ = {P*Pa_to_MPa:.0f} MPa,  "
    rf"$L$ = {Lx} m,  $E$ = {E*1e-9:.0f} GPa,  $\nu$ = {nu}",
    fontsize=13,
)
plt.tight_layout()

plot_path = os.path.join(CASE_DIR, "output", "bar_tension.png")
plt.savefig(plot_path, dpi=200)
print(f"[INFO] Plot saved to: {plot_path}")

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
sigma_xx_err = float(np.max(np.abs(sigma_xx - sigma_xx_ref)) / P)
sigma_yy_err = float(np.max(np.abs(sigma_yy)) / P)
sigma_zz_err = float(np.max(np.abs(sigma_zz)) / P)
sigma_xy_err = float(np.max(np.abs(sigma_xy)) / P)
u_xL_err = float(abs(u_xL_num - u_xL_ref) / u_xL_ref)

print("\n[RESULT] Relative errors vs engineering-bar analytical:")
print(f"   sigma_xx ~ P     : {sigma_xx_err:.3e}")
print(f"   sigma_yy ~ 0     : {sigma_yy_err:.3e}")
print(f"   sigma_zz ~ 0     : {sigma_zz_err:.3e}")
print(f"   sigma_xy ~ 0     : {sigma_xy_err:.3e}")
print(f"   u_x(L)  ~ P L/E  : {u_xL_err:.3e}")

errors = {
    "sigma_xx_max_err": {
        "numerical": float(np.max(sigma_xx)),
        "reference": float(sigma_xx_ref),
        "abs_error": float(np.max(np.abs(sigma_xx - sigma_xx_ref))),
        "rel_error": sigma_xx_err,
    },
    "sigma_yy_max_abs": {
        "numerical": float(np.max(np.abs(sigma_yy))),
        "reference": 0.0,
        "abs_error": float(np.max(np.abs(sigma_yy))),
        "rel_error": sigma_yy_err,
    },
    "sigma_zz_max_abs": {
        "numerical": float(np.max(np.abs(sigma_zz))),
        "reference": 0.0,
        "abs_error": float(np.max(np.abs(sigma_zz))),
        "rel_error": sigma_zz_err,
    },
    "u_xL": {
        "numerical": u_xL_num,
        "reference": u_xL_ref,
        "abs_error": abs(u_xL_num - u_xL_ref),
        "rel_error": u_xL_err,
    },
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
