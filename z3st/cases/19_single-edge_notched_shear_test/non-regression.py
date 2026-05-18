#!/usr/bin/env python3
"""Non-regression / diagnostic script for the SENS shear test.

Reproduces Ambati et al. 2015 (Comput Mech 55:383-405) §4.2 with AT2 +
hybrid: τ_xy / γ response at the notch tip, crack profile along the
expected curved-crack diagonal, global energy balance, and ParaView-style
2D field plots of D, σ_xy, σ_xx, σ_vm, and the crack-driving force H.
"""

import os
import re
from glob import glob

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
import yaml

from z3st.utils.utils_extract_vtu import extract_field

# ----- configuration --------------------------------------------------------
CASE_DIR    = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(CASE_DIR, "output")
VTU_FILES   = sorted(glob(os.path.join(OUTPUT_DIR, "fields_*.vtu")))

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
G_shear = E / (2.0 * (1.0 + nu))

# AT2/AT1 analytical sigma_c (matches spine.py reconciliation at material load).
if dmg_type == "AT2":
    sigma_c = float(np.sqrt(27.0 / 256.0 * Gc * E / lc))
elif dmg_type == "AT1":
    sigma_c = float(np.sqrt(3.0 / 8.0 * Gc * E / lc))
else:
    raise ValueError(f"Unknown damage type: {dmg_type}")

with open(os.path.join(CASE_DIR, "mesh.geo")) as f:
    Dn = float(re.search(r"Dn\s*=\s*([\d\.eE+-]+)", f.read()).group(1))

print(f"[INFO] Geometry  : Lx = {Lx*1e3:.3f} mm, Ly = {Ly*1e3:.3f} mm, Dn = {Dn*1e3:.3f} mm")
print(f"[INFO] Material  : {os.path.basename(MATERIAL_FILE)} (E={E:.2e} Pa, nu={nu}, Gc={Gc} J/m^2)")
print(f"[INFO] Phase-fld : {dmg_type}, lc = {lc*1e6:.2f} um -> sigma_c = {sigma_c/1e9:.3f} GPa")

# ----- per-step extraction at the notch tip ---------------------------------
# Sample a small disc around the notch tip (Dn, Ly/2). Shear stress σ_xy is
# index 1 of the row-major flattened 3x3 stress tensor.
x_target, y_target = Dn, Ly / 2.0
mask_tol = Lx / 80.0

shear_stress, shear_strain, u_x_top, d_max_list, h_max_list = [], [], [], [], []
for step, vtufile in enumerate(VTU_FILES):
    x_S, y_S, _, S_all = extract_field(vtufile, field_name="Stress_steel (cells)")
    mask_S = (np.abs(x_S - x_target) < mask_tol) & (np.abs(y_S - y_target) < mask_tol)
    shear_stress.append(float(np.mean(S_all[mask_S, 1])))

    x_E, y_E, _, E_all = extract_field(vtufile, field_name="Strain (cells)")
    mask_E = (np.abs(x_E - x_target) < mask_tol) & (np.abs(y_E - y_target) < mask_tol)
    shear_strain.append(float(np.mean(E_all[mask_E, 1])))

    x_u, y_u, _, u_all = extract_field(vtufile, field_name="Displacement")
    u_x_top.append(float(np.max(u_all[np.abs(y_u - Ly) < mask_tol, 0])))

    _, _, _, D_all = extract_field(vtufile, field_name="Damage")
    d_max_list.append(float(np.max(D_all)))

    _, _, _, H_all = extract_field(vtufile, field_name="CrackDrivingForce")
    h_max_list.append(float(np.max(H_all)))


# ----- shear response: τ_xy(tip) vs γ_remote --------------------------------
gamma_remote = np.array(u_x_top) / Ly
plt.figure(figsize=(7, 5))
plt.plot(gamma_remote, np.array(shear_stress) / 1e6, "C0-o", markersize=3, label=r"$\tau_{xy}$ at notch tip")
plt.plot(gamma_remote, G_shear * gamma_remote / 1e6, "k--", lw=1.5, label=r"$G\,\gamma$ (linear elastic)")
plt.xlabel(r"Remote shear strain $\gamma = u_x / L_y$")
plt.ylabel(r"$\tau_{xy}$ at notch tip (MPa)")
plt.title("SENS shear response (cf. Ambati Fig. 13)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "shear_response.png"), dpi=200)
plt.close()
print("[INFO] shear_response.png saved")


# ----- crack profile along the notch-tip -> bottom-right diagonal ------------
plt.figure(figsize=(8, 6))
band_half = 4.0 * lc
diag_x0, diag_y0 = Dn,  Ly / 2.0
diag_x1, diag_y1 = Lx,  0.0
dx, dy = diag_x1 - diag_x0, diag_y1 - diag_y0
diag_len = float(np.hypot(dx, dy))
ux_d, uy_d = dx / diag_len, dy / diag_len
nx_d, ny_d = -uy_d, ux_d

colors = plt.cm.jet(np.linspace(0, 1, len(VTU_FILES)))
for i, vtufile in enumerate(VTU_FILES):
    x_d, y_d, _, D_all = extract_field(vtufile, field_name="Damage")
    s_along = (x_d - diag_x0) * ux_d + (y_d - diag_y0) * uy_d
    s_perp  = (x_d - diag_x0) * nx_d + (y_d - diag_y0) * ny_d
    band = (np.abs(s_perp) < band_half) & (s_along >= 0) & (s_along <= diag_len)
    if not np.any(band):
        continue
    order = np.argsort(s_along[band])
    plt.plot(s_along[band][order] * 1e3, D_all[band][order],
             color=colors[i], lw=1.0, alpha=0.7,
             label=f"Step {i}" if i % max(1, len(VTU_FILES)//5) == 0 else None)
plt.xlabel("Distance along notch-tip -> bottom-right diagonal (mm)")
plt.ylabel(r"Damage $D$")
plt.title(f"Crack profile along the curved-crack diagonal ({2*band_half*1e6:.0f} um band)")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(loc="upper right", fontsize="small")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "crack_profile_diagonal.png"), dpi=200)
plt.close()
print("[INFO] crack_profile_diagonal.png saved")


# ----- global energy balance -------------------------------------------------
energy_file = os.path.join(CASE_DIR, "energies.txt")
if os.path.exists(energy_file):
    data = np.genfromtxt(energy_file, names=True, skip_header=0)
    ambati_arc = 0.55e-3   # Ambati Fig. 12d hybrid arrest length
    plt.figure(figsize=(8, 5))
    plt.plot(data["Step"], data["E_el"],   "b-o", markersize=3, label=r"Elastic $E_{el}$")
    plt.plot(data["Step"], data["E_frac"], "r-s", markersize=3, label=r"Fracture $E_{frac}$")
    plt.plot(data["Step"], data["E_tot"],  "k--", lw=1.5, label=r"Total $E_{tot}$")
    plt.axhline(Gc * ambati_arc, color="gray", ls=":", alpha=0.6,
                label=rf"Gc * 0.55 mm = {Gc*ambati_arc:.2f} J (Ambati Fig. 12d hybrid arrest)")
    plt.xlabel("Step")
    plt.ylabel("Energy (J)")
    plt.title("Global energy balance")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend(loc="upper left", fontsize="small")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "energy_balance.png"), dpi=200)
    plt.close()
    print("[INFO] energy_balance.png saved")


# ----- damage/H/τ at the notch tip ------------------------------------------
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
ax2.set_ylabel(r"$\tau_{xy}$ at notch tip (MPa)", color="tab:blue")
lns3 = ax2.plot(steps, np.array(shear_stress) / 1e6, "b-s", markersize=4, label=r"$\tau_{xy}$ (notch tip)")
ax2.axhline(sigma_c / 1e6, color="black", ls=":", alpha=0.5,
            label=rf"Analytic $\sigma_c$ ({dmg_type}) = {sigma_c/1e9:.2f} GPa")
ax2.tick_params(axis="y", labelcolor="tab:blue")
ax1.legend(lns1 + lns2 + lns3, [l.get_label() for l in lns1 + lns2 + lns3],
           loc="center left", frameon=True)
plt.title(f"Crack initiation and softening at the notch tip (x = {x_target*1e3:.3f} mm)")
fig.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "damage_evolution.png"), dpi=200)
plt.close()
print("[INFO] damage_evolution.png saved")


# ----- ParaView-style 2D field plots ----------------------------------------
def _build_triangulation(pv_mesh):
    pts = pv_mesh.points
    x, y = pts[:, 0], pts[:, 1]
    cells = pv_mesh.cells_dict
    if 5 in cells:
        tri = np.asarray(cells[5])
    elif 9 in cells:
        q = np.asarray(cells[9])
        tri = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    else:
        tri = None
    return (mtri.Triangulation(x, y, tri) if tri is not None and len(tri) > 0
            else mtri.Triangulation(x, y)), x, y


def _draw_notch(ax):
    ax.plot([0.0, Dn], [Ly / 2.0, Ly / 2.0], color="cyan", linewidth=2.5,
            label=f"Pre-crack slit (0 -- {Dn*1e3:.2f} mm at y = Ly/2)")


def plot_field(triang, field, *, output_path, title, cbar_label, cmap, vmin, vmax, mark_notch=True):
    levels = np.linspace(vmin, vmax, 41)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    cf = ax.tricontourf(triang, field, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, extend="both")
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


def von_mises(s):
    sxx, sxy, sxz = s[:, 0], s[:, 1], s[:, 2]
    syy, syz      = s[:, 4], s[:, 5]
    szz           = s[:, 8]
    return np.sqrt(0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2
                          + 6.0 * (sxy**2 + syz**2 + sxz**2)))


def _get(pv_mesh, names):
    for n in names:
        if n in pv_mesh.point_data:
            return n, np.asarray(pv_mesh.point_data[n])
    return None, None


try:
    last = pv.read(VTU_FILES[-1])
    triang, x_pts, y_pts = _build_triangulation(last)
    src = os.path.basename(VTU_FILES[-1])

    dname, D_field = _get(last, ["Damage", "damage", "D"])
    if D_field is not None:
        plot_field(triang, D_field.reshape(-1),
                   output_path=os.path.join(OUTPUT_DIR, "damage_field.png"),
                   title=f"Damage field D (Ambati Fig. 12 hybrid reproducer)  -- {src}",
                   cbar_label="Damage D", cmap="hot_r", vmin=0.0, vmax=1.0)
        print("[INFO] damage_field.png saved")

        n = len(VTU_FILES)
        if n >= 4:
            idx = sorted(set(np.linspace(0, n - 1, 4, dtype=int).tolist()))
            fig, axes = plt.subplots(1, 4, figsize=(20, 5.5))
            for ax, i in zip(axes, idx):
                m = pv.read(VTU_FILES[i])
                tri, _, _ = _build_triangulation(m)
                D = np.asarray(m.point_data[dname]).reshape(-1)
                cf = ax.tricontourf(tri, D, levels=np.linspace(0, 1, 41), cmap="hot_r",
                                    vmin=0, vmax=1, extend="both")
                ax.plot([0.0, Dn], [Ly/2.0, Ly/2.0], color="cyan", lw=2)
                ax.set_aspect("equal")
                ax.set_xlim(-0.05 * Lx, 1.05 * Lx)
                ax.set_ylim(-0.05 * Ly, 1.05 * Ly)
                ax.set_xlabel("x (m)")
                if ax is axes[0]:
                    ax.set_ylabel("y (m)")
                ax.set_title(f"Step {i}")
            fig.colorbar(cf, ax=axes.tolist(), label="Damage D",
                         orientation="horizontal", fraction=0.06, pad=0.12)
            plt.suptitle("Damage evolution across the loading window", fontsize=13)
            plt.savefig(os.path.join(OUTPUT_DIR, "damage_evolution_panels.png"), dpi=180)
            plt.close(fig)
            print("[INFO] damage_evolution_panels.png saved")

    _, S = _get(last, ["Stress_steel (points)", "Stress (points)"])
    if S is not None and S.ndim == 2 and S.shape[1] >= 9:
        vm = von_mises(S) / 1e6
        vm_hi = max(float(np.nanpercentile(vm, 99.0)), 1.0)
        plot_field(triang, vm,
                   output_path=os.path.join(OUTPUT_DIR, "stress_vm_field.png"),
                   title="Von Mises equivalent stress (99th-pct clip)",
                   cbar_label="sigma_vm (MPa)", cmap="viridis", vmin=0.0, vmax=vm_hi)
        print("[INFO] stress_vm_field.png saved")

        sxy = S[:, 1] / 1e6
        sxy_abs = max(float(np.nanpercentile(np.abs(sxy), 99.0)), 1.0)
        plot_field(triang, sxy,
                   output_path=os.path.join(OUTPUT_DIR, "stress_xy_field.png"),
                   title="Shear stress sigma_xy (Mode-II driver)",
                   cbar_label="sigma_xy (MPa)", cmap="seismic", vmin=-sxy_abs, vmax=+sxy_abs)
        print("[INFO] stress_xy_field.png saved")

        sxx = S[:, 0] / 1e6
        sxx_abs = max(float(np.nanpercentile(np.abs(sxx), 99.0)), 1.0)
        plot_field(triang, sxx,
                   output_path=os.path.join(OUTPUT_DIR, "stress_xx_field.png"),
                   title="Normal stress sigma_xx (shear-diagonal pattern)",
                   cbar_label="sigma_xx (MPa)", cmap="seismic", vmin=-sxx_abs, vmax=+sxx_abs)
        print("[INFO] stress_xx_field.png saved")

    _, H_arr = _get(last, ["CrackDrivingForce", "crack_driving_force", "H"])
    if H_arr is not None:
        H_arr = H_arr.reshape(-1)
        H_hi = max(float(np.nanpercentile(H_arr, 99.0)), 1e-30)
        plot_field(triang, H_arr,
                   output_path=os.path.join(OUTPUT_DIR, "crack_driving_force_field.png"),
                   title="Crack-driving force H (non-dim. AT2: H = (2 lc/Gc) psi+)",
                   cbar_label="H", cmap="magma", vmin=0.0, vmax=H_hi)
        print("[INFO] crack_driving_force_field.png saved")
except Exception as e:
    print(f"[WARN] 2D field plotting failed: {e}")


print("\n[SUMMARY] Ambati Fig. 12 (hybrid): curved Mode-II crack from notch tip into the")
print("          lower-right region, arresting around y ~ 0.1 mm (arc ~ 0.55 mm).")
print(f"          Final max D = {d_max_list[-1]:.4f}  (target: 1.0)")
print(f"          Gc * 0.55 mm = {Gc * 0.55e-3:.2f} J at arrest.")
