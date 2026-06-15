#!/usr/bin/env python3
"""Quick diagnostic for the SENT tension test.

Processes only the **last** available VTU file (no per-step loop), so it
runs in a few seconds and is safe to invoke while the simulation is still
writing new steps. Generates the ParaView-style field plots of the final
state (damage, sigma_yy, sigma_xx, sigma_vm, crack-driving force) plus the
global energy balance read directly from energies.txt.

For per-step diagnostics (F-u curve, damage evolution panels), use a
separate post-processing script after the simulation completes.
"""

import os
import re
from glob import glob

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
import yaml

# ----- configuration --------------------------------------------------------
CASE_DIR   = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
VTU_FILES  = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))
if not VTU_FILES:
    raise FileNotFoundError(f"No VTU files found in {OUTPUT_DIR}")
LAST_VTU = VTU_FILES[-1]

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    cfg = yaml.safe_load(f)
mat_rel = list(cfg["materials"].values())[0]
MATERIAL_FILE = os.path.normpath(os.path.join(CASE_DIR, mat_rel))

dmg_cfg  = cfg.get("damage", {})
lc       = float(dmg_cfg["lc"])
dmg_type = str(dmg_cfg.get("type", "AT2")).upper()

with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lx, Ly = float(geom["Lx"]), float(geom["Ly"])

with open(MATERIAL_FILE) as f:
    mat = yaml.safe_load(f)
E, nu, Gc = float(mat["E"]), float(mat["nu"]), float(mat["Gc"])

# AT2/AT1 analytical sigma_c (matches spine.py reconciliation at material load).
if dmg_type == "AT2":
    sigma_c = float(np.sqrt(27.0 / 256.0 * Gc * E / lc))
elif dmg_type == "AT1":
    sigma_c = float(np.sqrt(3.0 / 8.0 * Gc * E / lc))
else:
    raise ValueError(f"Unknown damage type: {dmg_type}")

with open(os.path.join(CASE_DIR, "mesh.geo")) as f:
    Dn = float(re.search(r"Dn\s*=\s*([\d\.eE+-]+)", f.read()).group(1))

print(f"[INFO] Case      : {os.path.basename(CASE_DIR)}")
print(f"[INFO] Last VTU  : {os.path.basename(LAST_VTU)}  (of {len(VTU_FILES)} files)")
print(f"[INFO] Geometry  : Lx = {Lx*1e3:.3f} mm, Ly = {Ly*1e3:.3f} mm, Dn = {Dn*1e3:.3f} mm")
print(f"[INFO] Material  : {os.path.basename(MATERIAL_FILE)} (E={E:.2e} Pa, nu={nu}, Gc={Gc} J/m^2)")
print(f"[INFO] Phase-fld : {dmg_type}, lc = {lc*1e6:.2f} um -> sigma_c = {sigma_c/1e9:.3f} GPa")


# ----- ParaView-style 2D field plots (last VTU only) ------------------------
def _build_triangulation(pv_mesh):
    pts = pv_mesh.points
    x, y = pts[:, 0], pts[:, 1]
    tri = None
    try:
        cells = pv_mesh.cells_dict
        if 5 in cells:                        # VTK_TRIANGLE
            tri = np.asarray(cells[5])
        elif 9 in cells:                      # VTK_QUAD -> split into two tris
            q = np.asarray(cells[9])
            tri = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    except Exception:
        tri = None
    if tri is None:
        # dolfinx 0.11 writes VTK_LAGRANGE_TRIANGLE; cells_dict raises. Rebuild
        # from the raw connectivity (P1 -> 3 corner nodes per cell). Required here:
        # the notched domain is non-convex, so a Delaunay fallback would bridge
        # the slit.
        try:
            conn = np.asarray(pv_mesh.cell_connectivity)
            sizes = np.diff(np.asarray(pv_mesh.offset))
            if sizes.size and np.all(sizes == 3):
                tri = conn.reshape(-1, 3)
        except Exception:
            tri = None
    return (mtri.Triangulation(x, y, tri) if tri is not None and len(tri) > 0
            else mtri.Triangulation(x, y)), x, y


def _draw_notch(ax):
    ax.plot([0.0, Dn], [Ly / 2.0, Ly / 2.0], color="cyan", linewidth=2.5,
            label=f"Pre-crack slit (0 -- {Dn*1e3:.2f} mm at y = Ly/2)")


def plot_field(triang, field, *, output_path, title, cbar_label, cmap, vmin, vmax):
    levels = np.linspace(vmin, vmax, 41)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    cf = ax.tricontourf(triang, field, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, extend="both")
    ax.triplot(triang, color="black", linewidth=0.05, alpha=0.3)
    _draw_notch(ax)
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(cf, ax=ax, label=cbar_label)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlim(-0.05 * Lx, 1.05 * Lx)
    ax.set_ylim(-0.05 * Ly, 1.05 * Ly)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close(fig)


def von_mises(s):
    sxx, sxy, sxz = s[:, 0], s[:, 1], s[:, 2]
    syy, syz      = s[:, 4], s[:, 5]
    szz           = s[:, 8]
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2
                          + 6.0 * (sxy ** 2 + syz ** 2 + sxz ** 2)))


def _get(pv_mesh, names):
    for n in names:
        if n in pv_mesh.point_data:
            return np.asarray(pv_mesh.point_data[n])
    return None


m = pv.read(LAST_VTU)
triang, _, _ = _build_triangulation(m)
src = os.path.basename(LAST_VTU)

D_field = _get(m, ["Damage", "damage", "D"])
if D_field is not None:
    plot_field(triang, D_field.reshape(-1),
               output_path=os.path.join(OUTPUT_DIR, "damage_field.png"),
               title=f"Damage field D (Ambati Fig. 8c reproducer)  -- {src}",
               cbar_label="Damage D", cmap="hot_r", vmin=0.0, vmax=1.0)
    print("[INFO] damage_field.png saved")
    print(f"          max D = {float(np.max(D_field)):.4f}")

S = _get(m, ["Stress (points)", "Stress (points)"])
if S is not None and S.ndim == 2 and S.shape[1] >= 9:
    vm = von_mises(S) / 1e6
    vm_hi = max(float(np.nanpercentile(vm, 99.0)), 1.0)
    plot_field(triang, vm,
               output_path=os.path.join(OUTPUT_DIR, "stress_vm_field.png"),
               title="Von Mises equivalent stress (99th-pct clip)",
               cbar_label="sigma_vm (MPa)", cmap="viridis", vmin=0.0, vmax=vm_hi)
    print("[INFO] stress_vm_field.png saved")

    syy = S[:, 4] / 1e6
    syy_abs = max(float(np.nanpercentile(np.abs(syy), 99.0)), 1.0)
    plot_field(triang, syy,
               output_path=os.path.join(OUTPUT_DIR, "stress_yy_field.png"),
               title="Tensile stress sigma_yy (Mode-I driver from notch tip)",
               cbar_label="sigma_yy (MPa)", cmap="seismic", vmin=-syy_abs, vmax=+syy_abs)
    print("[INFO] stress_yy_field.png saved")

    sxx = S[:, 0] / 1e6
    sxx_abs = max(float(np.nanpercentile(np.abs(sxx), 99.0)), 1.0)
    plot_field(triang, sxx,
               output_path=os.path.join(OUTPUT_DIR, "stress_xx_field.png"),
               title="Transverse stress sigma_xx (Poisson-blocked at corners)",
               cbar_label="sigma_xx (MPa)", cmap="seismic", vmin=-sxx_abs, vmax=+sxx_abs)
    print("[INFO] stress_xx_field.png saved")

H_arr = _get(m, ["CrackDrivingForce", "crack_driving_force", "H"])
if H_arr is not None:
    H_arr = H_arr.reshape(-1)
    H_hi = max(float(np.nanpercentile(H_arr, 99.0)), 1e-30)
    plot_field(triang, H_arr,
               output_path=os.path.join(OUTPUT_DIR, "crack_driving_force_field.png"),
               title="Crack-driving force H (non-dim. AT2: H = (2 lc/Gc) psi+)",
               cbar_label="H", cmap="magma", vmin=0.0, vmax=H_hi)
    print("[INFO] crack_driving_force_field.png saved")


# ----- global energy balance from energies.txt ------------------------------
# For SENT (Mode-I straight crack), the expected final E_frac is Gc * Lx,
# i.e. full ligament traversal of the prescribed notch (0.5 mm) plus the
# propagated arc (0.5 mm) -- the crack reaches the right edge in Ambati's
# Fig. 8c. The step-0 baseline is Gc * Dn (notch only, before propagation).
energy_file = os.path.join(CASE_DIR, "energies.txt")
if os.path.exists(energy_file):
    data = np.genfromtxt(energy_file, names=True, skip_header=0)
    E_frac_notch    = Gc * Dn
    E_frac_full_lig = Gc * Lx
    plt.figure(figsize=(8, 5))
    plt.plot(data["Step"], data["E_el"],   "b-o", markersize=3, label=r"Elastic $E_{el}$")
    plt.plot(data["Step"], data["E_frac"], "r-s", markersize=3, label=r"Fracture $E_{frac}$")
    plt.plot(data["Step"], data["E_tot"],  "k--", lw=1.5, label=r"Total $E_{tot}$")
    plt.axhline(E_frac_notch, color="gray", ls="--", alpha=0.5,
                label=rf"Gc * Dn = {E_frac_notch:.2f} J (notch baseline, step 0)")
    plt.axhline(E_frac_full_lig, color="gray", ls=":", alpha=0.7,
                label=rf"Gc * Lx = {E_frac_full_lig:.2f} J (full ligament traversal, Ambati Fig. 8c)")
    plt.xlabel("Step")
    plt.ylabel("Energy (J)")
    plt.title("Global energy balance")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend(loc="upper left", fontsize="small")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "energy_balance.png"), dpi=200)
    plt.close()
    print("[INFO] energy_balance.png saved")
    print(f"          final  E_el = {data['E_el'][-1]:.3f} J, "
          f"E_frac = {data['E_frac'][-1]:.3f} J  (step {int(data['Step'][-1])})")
else:
    print(f"[WARN] {energy_file} not found; skipping energy plot.")
