#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 19_SENT_2D (single-edge notched tension, 2D plane strain).

Reproduces the Ambati et al. 2015 (Comput Mech 55:383-405) §4.1 SENT
tension benchmark using the AT2 hybrid phase-field formulation (their
Eq. 27). Processes the per-step VTU snapshots and reconstructs:
  - sigma_yy vs eps_yy at the notch tip (the canonical Ambati Fig. 6
    force-displacement curve, expressed as macroscopic stress-strain),
  - crack-profile evolution D(y) at x = Dn/2 (cf. Ambati Fig. 5),
  - global energy balance E_el(t), E_frac(t),
  - notch-tip damage and crack-driving-force history.
"""

import os
import re
from glob import glob

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
import yaml

from z3st.utils.utils_extract_vtu import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR    = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(CASE_DIR, "output")
VTU_FILES   = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))

INPUT_FILE     = os.path.join(CASE_DIR, "input.yaml")
GEOMETRY_FILE  = os.path.join(CASE_DIR, "geometry.yaml")
BC_FILE        = os.path.join(CASE_DIR, "boundary_conditions.yaml")
MESH_GEO_FILE  = os.path.join(CASE_DIR, "mesh.geo")

# Resolve the actual material card declared in input.yaml (not a hardcoded
# default). The previous version of this script read steel.yaml while the
# simulation ran high_carbon_steel.yaml -- the resulting analytical curve
# and sigma_c reference line were both wrong.
with open(INPUT_FILE, "r") as f:
    input_data = yaml.safe_load(f)
mat_rel_path = list(input_data["materials"].values())[0]
MATERIAL_FILE = os.path.normpath(os.path.join(CASE_DIR, mat_rel_path))

dmg_cfg = input_data.get("damage", {})
lc = float(dmg_cfg["lc"])
dmg_type = str(dmg_cfg.get("type", "AT2")).upper()

with open(GEOMETRY_FILE, "r") as f:
    geom_data = yaml.safe_load(f)
Lx = float(geom_data["Lx"])
Ly = float(geom_data["Ly"])

with open(MATERIAL_FILE, "r") as f:
    mat_data = yaml.safe_load(f)
E  = float(mat_data["E"])
nu = float(mat_data["nu"])
Gc = float(mat_data["Gc"])

# AT2 / AT1 sigma_c reconciliation -- the framework's spine.py applies the
# same identities at material-load time. We replicate them here for the
# diagnostic threshold line in damage_evolution.png.
if dmg_type == "AT2":
    sigma_c = float(np.sqrt(27.0 / 256.0 * Gc * E / lc))
elif dmg_type == "AT1":
    sigma_c = float(np.sqrt(3.0 / 8.0 * Gc * E / lc))
else:
    raise ValueError(f"Unknown damage type: {dmg_type}")

with open(MESH_GEO_FILE, "r") as f:
    geo_content = f.read()
Dn = float(re.search(r"Dn\s*=\s*([\d\.eE+-]+)", geo_content).group(1))

print(f"[INFO] Geometry  : Lx = {Lx*1e3:.3f} mm, Ly = {Ly*1e3:.3f} mm, Dn = {Dn*1e3:.3f} mm")
print(f"[INFO] Material  : {os.path.basename(MATERIAL_FILE)} (E = {E:.2e} Pa, nu = {nu}, Gc = {Gc} J/m^2)")
print(f"[INFO] Phase-fld : {dmg_type}, lc = {lc*1e6:.2f} um -> sigma_c (analytic) = {sigma_c/1e9:.3f} GPa")

# --.. ..- .-.. .-.. --- analytical reference --.. ..- .-.. .-.. ---
# Plane-strain uniaxial extension at the centreline: sigma_yy = E/(1-nu^2) * eps_yy.
STRESSES_REF = np.array([0.0, 2.0e8, 4.0e8, 5.5e8, 6.0e8])
STRAINS_REF  = (1.0 - nu**2) / E * STRESSES_REF
U_X_REF      = STRAINS_REF * Lx

# Initial fracture energy (per unit out-of-plane thickness): Gc * notch_length.
E_frac_initial = Gc * Dn * 1.0
print(f"[INFO] Initial fracture energy (Gc * Dn * 1 m) = {E_frac_initial:.2f} J")

TOLERANCE = 1.0e-2

# --.. ..- .-.. .-.. --- per-step extraction --.. ..- .-.. .-.. ---
# Sample on a thin vertical strip through the notch tip x = Dn.
x_target = Dn
mask_tol = Lx / 80.0

strains, stresses, displacements, d_max_list, h_max_list = [], [], [], [], []
for step, vtufile in enumerate(VTU_FILES):
    print(f"\n[STEP {step}] {os.path.basename(vtufile)}")

    x_S, _, _, S_all = extract_field(vtufile, field_name="Stress_steel (cells)")
    mask_S = np.abs(x_S - x_target) < mask_tol
    stresses.append(float(np.mean(S_all[mask_S, 4])))   # sigma_yy

    x_u, _, _, u_all = extract_field(vtufile, field_name="Displacement")
    mask_u = np.abs(x_u - x_target) < mask_tol
    displacements.append(float(np.max(u_all[mask_u, 1])))

    x_E, _, _, E_all = extract_field(vtufile, field_name="Strain (cells)")
    mask_E = np.abs(x_E - x_target) < mask_tol
    strains.append(float(np.mean(E_all[mask_E, 4])))    # eps_yy

    _, _, _, D_all = extract_field(vtufile, field_name="Damage")
    d_max_list.append(float(np.max(D_all)))

    _, _, _, H_all = extract_field(vtufile, field_name="CrackDrivingForce")
    h_max_list.append(float(np.max(H_all)))

strains_np  = np.array(strains)
stresses_np = np.array(stresses)
displ_x_np  = np.array(displacements)

print("\n--. step  eps_yy        sigma_yy       u_y_max")
for i, (e, s, u) in enumerate(zip(strains, stresses, displacements)):
    print(f"  {i:3d}   {e:.3e}   {s:.3e}   {u:.3e}")

# --.. ..- .-.. .-.. --- crack profile evolution D(y) at x = Dn/2 --.. ..- .-.. .-.. ---
plt.figure(figsize=(8, 6))
x_line = Dn / 2.0
tol_x  = Lx / 100.0
yc     = Ly / 2.0
colors = plt.cm.jet(np.linspace(0, 1, len(VTU_FILES)))

y_plot_last = None
for i, vtufile in enumerate(VTU_FILES):
    x_d, y_d, _, D_all = extract_field(vtufile, field_name="Damage")
    mask_line = np.abs(x_d - x_line) < tol_x
    if not np.any(mask_line):
        continue
    y_coords = y_d[mask_line]
    d_values = D_all[mask_line]
    sort_idx = np.argsort(y_coords)
    y_plot = y_coords[sort_idx]
    d_plot = d_values[sort_idx]
    plt.plot(y_plot, d_plot, color=colors[i], lw=1.0, alpha=0.7, label=f"Step {i}" if i % 20 == 0 else None)
    y_plot_last = y_plot

if y_plot_last is not None:
    # AT2 analytical 1D profile around a fully-formed crack: D(y) = exp(-|y-yc|/lc).
    # (AT1 would be D(y) = (1 - |y-yc|/(2 lc))^2 within |y-yc| < 2 lc, else 0.)
    if dmg_type == "AT2":
        d_analytical = np.exp(-np.abs(y_plot_last - yc) / lc)
        label_a = r"AT2 analytical $\exp(-|y-y_c|/\ell_c)$"
    else:
        d_analytical = np.clip(1.0 - np.abs(y_plot_last - yc) / (2.0 * lc), 0.0, 1.0) ** 2
        label_a = r"AT1 analytical $(1 - |y-y_c|/(2\ell_c))_+^2$"
    plt.plot(y_plot_last, d_analytical, "k--", lw=2, label=label_a, zorder=10)

plt.xlabel("y coordinate (m)")
plt.ylabel(r"Damage $D$ (-)")
plt.title(f"Crack profile evolution at x = Dn/2 = {x_line*1e3:.3f} mm  ($\\ell_c = {lc*1e6:.1f}$ um)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(loc="upper right", fontsize="small", ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "crack_profile_evolution.png"), dpi=200)
plt.close()
print(f"\n[INFO] crack_profile_evolution.png saved ({len(VTU_FILES)} steps).")

# --.. ..- .-.. .-.. --- energy balance --.. ..- .-.. .-.. ---
energy_file = os.path.join(CASE_DIR, "energies.txt")
if os.path.exists(energy_file):
    data = np.genfromtxt(energy_file, names=True, skip_header=0)
    plt.figure(figsize=(8, 5))
    plt.plot(data["Step"], data["E_el"],   "b-o", markersize=3, label=r"Elastic $E_{el}$")
    plt.plot(data["Step"], data["E_frac"], "r-s", markersize=3, label=r"Fracture $E_{frac}$")
    plt.plot(data["Step"], data["E_tot"],  "k--", lw=1.5, label=r"Total $E_{tot}$")
    plt.axhline(Gc * Lx, color="gray", ls=":", alpha=0.6,
                label=rf"Gc * Lx = {Gc*Lx:.2f} J (full ligament traversal)")
    plt.xlabel("Step")
    plt.ylabel("Energy (J)")
    plt.title("Global energy balance")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend(loc="upper left", fontsize="small")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "energy_balance.png"), dpi=200)
    plt.close()
    print("[INFO] energy_balance.png saved")
else:
    print(f"[WARN] {energy_file} not found; skipping energy plot.")

# --.. ..- .-.. .-.. --- stress-strain curve --.. ..- .-.. .-.. ---
plt.figure(figsize=(7, 5))
plt.plot(strains, np.array(stresses) / 1e6, "C0--o", markersize=4, label="Numerical (notch tip)")
plt.plot(STRAINS_REF, STRESSES_REF / 1e6, "k-", lw=1.5, label=r"Plane-strain uniaxial $\sigma = E/(1-\nu^2)\,\epsilon$")
plt.xlabel(r"$\epsilon_{yy}$ (-)")
plt.ylabel(r"$\sigma_{yy}$ (MPa)")
plt.title("Stress-strain at the notch tip (x = Dn)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "stress_strain_curve.png"), dpi=200)
plt.close()
print("[INFO] stress_strain_curve.png saved")

# --.. ..- .-.. .-.. --- damage / driving-force history --.. ..- .-.. .-.. ---
steps = np.arange(len(VTU_FILES))
fig, ax1 = plt.subplots(figsize=(9, 6))

ax1.set_xlabel("Step")
ax1.set_ylabel(r"Damage $D$ / Normalised $H$", color="tab:red")
lns1 = ax1.plot(steps, d_max_list, "r-o", markersize=4, label=r"Max $D$")
h_norm = np.array(h_max_list) / max(np.max(h_max_list), 1e-30)
lns2 = ax1.plot(steps, h_norm, "g--", lw=1.5, label=r"Normalised $H$")
ax1.set_ylim(-0.05, 1.1)
ax1.tick_params(axis="y", labelcolor="tab:red")
ax1.grid(True, ls=":", alpha=0.6)

ax2 = ax1.twinx()
ax2.set_ylabel(r"$\sigma_{yy}$ at notch tip (MPa)", color="tab:blue")
lns3 = ax2.plot(steps, np.array(stresses) / 1e6, "b-s", markersize=4, label=r"$\sigma_{yy}$ (notch tip)")
ax2.axhline(sigma_c / 1e6, color="black", ls=":", alpha=0.5,
            label=rf"Analytic $\sigma_c$ ({dmg_type}) = {sigma_c/1e9:.2f} GPa")
ax2.tick_params(axis="y", labelcolor="tab:blue")

lns = lns1 + lns2 + lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc="center left", frameon=True)

plt.title(f"Crack initiation and softening at the notch tip (x = {x_target*1e3:.3f} mm)")
fig.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "damage_evolution.png"), dpi=200)
plt.close()
print("[INFO] damage_evolution.png saved")


# --.. ..- .-.. .-.. --- ParaView-like 2D field plots --.. ..- .-.. .-.. ---
# Same plotting helpers as case 14_full_cylinder_cracking_2D_xy. Renders
# filled-contour fields on the triangulation with the notch slit marked
# in cyan; saves at the final time and at a 4-step evolution panel.

def _build_triangulation(pv_mesh):
    """matplotlib Triangulation from a PyVista 2D mesh, handling tri/quad cells."""
    pts = pv_mesh.points
    x_pts, y_pts = pts[:, 0], pts[:, 1]
    cells_dict = pv_mesh.cells_dict
    triangles = None
    if 5 in cells_dict:                       # VTK_TRIANGLE
        triangles = np.asarray(cells_dict[5])
    elif 9 in cells_dict:                     # VTK_QUAD -> split into two triangles
        q = np.asarray(cells_dict[9])
        triangles = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    if triangles is not None and len(triangles) > 0:
        return mtri.Triangulation(x_pts, y_pts, triangles), x_pts, y_pts
    return mtri.Triangulation(x_pts, y_pts), x_pts, y_pts


def _draw_notch(ax, color="cyan", lw=2.5):
    """Overlay the pre-crack notch slit as a horizontal segment."""
    ax.plot([0.0, Dn], [Ly / 2.0, Ly / 2.0],
            color=color, linewidth=lw,
            label=f"Pre-crack slit (0 -- {Dn*1e3:.2f} mm at y = Ly/2)")


def plot_field_2d(triang, field, *, output_path, title, cbar_label, cmap,
                  vmin, vmax, show_mesh=True, levels=40, mark_notch=True):
    """ParaView-like 2D field plot: filled contours + optional mesh overlay + notch."""
    levels_arr = np.linspace(vmin, vmax, levels + 1)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    cf = ax.tricontourf(triang, field, levels=levels_arr, cmap=cmap,
                        vmin=vmin, vmax=vmax, extend="both")
    if show_mesh:
        ax.triplot(triang, color="black", linewidth=0.05, alpha=0.3)
    if mark_notch:
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


def von_mises_from_flat(sigma_flat):
    """Von Mises equivalent stress from an (N, 9) row-major flattened 3x3 sigma."""
    s_xx, s_xy, s_xz = sigma_flat[:, 0], sigma_flat[:, 1], sigma_flat[:, 2]
    s_yy, s_yz       = sigma_flat[:, 4], sigma_flat[:, 5]
    s_zz             = sigma_flat[:, 8]
    return np.sqrt(0.5 * (
        (s_xx - s_yy) ** 2 + (s_yy - s_zz) ** 2 + (s_zz - s_xx) ** 2
        + 6.0 * (s_xy ** 2 + s_yz ** 2 + s_xz ** 2)
    ))


def _get_point_array(pv_mesh, candidates):
    for name in candidates:
        if name in pv_mesh.point_data:
            return name, np.asarray(pv_mesh.point_data[name])
    return None, None


try:
    last_mesh_pv = pv.read(VTU_FILES[-1])
    triang, x_pts, y_pts = _build_triangulation(last_mesh_pv)
    src_name = os.path.basename(VTU_FILES[-1])

    # ------ Damage field (Ambati Fig. 5 reproducer) ------
    dname, D_field = _get_point_array(last_mesh_pv, ["Damage", "damage", "D"])
    if D_field is not None:
        D_field = D_field.reshape(-1)
        plot_field_2d(
            triang, D_field,
            output_path=os.path.join(OUTPUT_DIR, "damage_field.png"),
            title=(f"Damage field D at the final loading step\n"
                   f"(Ambati Fig. 5 reproducer)  -- Source: {src_name}"),
            cbar_label="Damage D", cmap="hot_r", vmin=0.0, vmax=1.0,
        )
        print(f"[INFO] damage_field.png saved")

        # ------ 4-step damage-evolution panel ------
        n_files = len(VTU_FILES)
        if n_files >= 4:
            panel_idx = sorted(set(np.linspace(0, n_files - 1, 4, dtype=int).tolist()))
            fig, axes = plt.subplots(1, 4, figsize=(20, 5.5))
            for ax, idx in zip(axes, panel_idx):
                m = pv.read(VTU_FILES[idx])
                tri, _, _ = _build_triangulation(m)
                D = np.asarray(m.point_data[dname]).reshape(-1)
                cf = ax.tricontourf(tri, D, levels=np.linspace(0, 1, 41),
                                    cmap="hot_r", vmin=0, vmax=1, extend="both")
                ax.plot([0.0, Dn], [Ly / 2.0, Ly / 2.0], color="cyan", lw=2)
                ax.set_aspect("equal")
                ax.set_xlim(-0.05 * Lx, 1.05 * Lx)
                ax.set_ylim(-0.05 * Ly, 1.05 * Ly)
                ax.set_xlabel("x (m)")
                if ax is axes[0]:
                    ax.set_ylabel("y (m)")
                ax.set_title(f"Step {idx}")
            fig.colorbar(cf, ax=axes.tolist(), label="Damage D",
                         orientation="horizontal", fraction=0.06, pad=0.12)
            plt.suptitle("Damage evolution across the loading window", fontsize=13)
            plt.savefig(os.path.join(OUTPUT_DIR, "damage_evolution_panels.png"), dpi=180)
            plt.close(fig)
            print("[INFO] damage_evolution_panels.png saved")

    # ------ Stress fields ------
    sname, sigma_flat = _get_point_array(
        last_mesh_pv, ["Stress_steel (points)", "Stress (points)"]
    )
    if sigma_flat is not None and sigma_flat.ndim == 2 and sigma_flat.shape[1] >= 9:
        # Von Mises
        vm_MPa = von_mises_from_flat(sigma_flat) / 1e6
        vm_high = max(float(np.nanpercentile(vm_MPa, 99.0)), 1.0)
        plot_field_2d(
            triang, vm_MPa,
            output_path=os.path.join(OUTPUT_DIR, "stress_vm_field.png"),
            title=(f"Von Mises equivalent stress at the final loading step\n"
                   f"colorbar clipped at the 99th percentile of |sigma_vm|"),
            cbar_label="sigma_vm (MPa)", cmap="viridis",
            vmin=0.0, vmax=vm_high,
        )
        print("[INFO] stress_vm_field.png saved")

        # sigma_yy (the Mode-I tensile driver of the crack)
        syy_MPa = sigma_flat[:, 4] / 1e6
        syy_abs = max(float(np.nanpercentile(np.abs(syy_MPa), 99.0)), 1.0)
        plot_field_2d(
            triang, syy_MPa,
            output_path=os.path.join(OUTPUT_DIR, "stress_yy_field.png"),
            title=(f"Tensile stress sigma_yy at the final loading step\n"
                   f"(positive = tensile; drives the Mode-I crack from the notch tip)"),
            cbar_label="sigma_yy (MPa)", cmap="seismic",
            vmin=-syy_abs, vmax=+syy_abs,
        )
        print("[INFO] stress_yy_field.png saved")

        # sigma_xx (transverse, used to identify crack-tip-shielding compressive zones)
        sxx_MPa = sigma_flat[:, 0] / 1e6
        sxx_abs = max(float(np.nanpercentile(np.abs(sxx_MPa), 99.0)), 1.0)
        plot_field_2d(
            triang, sxx_MPa,
            output_path=os.path.join(OUTPUT_DIR, "stress_xx_field.png"),
            title=(f"Transverse stress sigma_xx at the final loading step\n"
                   f"(compressive lobes at the crack tip indicate shielding)"),
            cbar_label="sigma_xx (MPa)", cmap="seismic",
            vmin=-sxx_abs, vmax=+sxx_abs,
        )
        print("[INFO] stress_xx_field.png saved")

    # ------ Crack driving force H ------
    hname, H_arr = _get_point_array(
        last_mesh_pv, ["CrackDrivingForce", "crack_driving_force", "H"]
    )
    if H_arr is not None:
        H_arr = H_arr.reshape(-1)
        H_high = max(float(np.nanpercentile(H_arr, 99.0)), 1e-30)
        plot_field_2d(
            triang, H_arr,
            output_path=os.path.join(OUTPUT_DIR, "crack_driving_force_field.png"),
            title=(f"Crack-driving force H at the final loading step\n"
                   f"(non-dimensional for AT2: H = (2*lc/Gc) * psi+;  "
                   f"colorbar clipped at 99th percentile)"),
            cbar_label="H (-)", cmap="magma",
            vmin=0.0, vmax=H_high,
        )
        print("[INFO] crack_driving_force_field.png saved")

except Exception as e:
    print(f"[WARN] 2D field plotting failed: {e}")


print("\n[SUMMARY] Reference Ambati Fig. 6: peak force at u_y ~ 5.5-6 um; brittle drop after.")
print(f"           Final max D = {d_max_list[-1]:.4f}  (target: 1.0)")
print(f"           Final E_frac (if logged) ~ Gc * (Lx - Dn) = {Gc * (Lx - Dn):.2f} J at full traversal.")
