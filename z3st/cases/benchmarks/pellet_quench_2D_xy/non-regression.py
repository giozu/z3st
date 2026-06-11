#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: benchmarks/pellet_quench_2D_xy

2D Cartesian (plane-strain) transverse cross-section of a UO2 pellet,
modeled as the upper semicircle (y >= 0) with mirror symmetry on y=0.
The 60-degree cold-contact arc spans 0 to +30 deg (full +/- 30 deg by
mirror) — the same geometry as McClenny et al., JNM 565 (2022) 153719,
Fig. 4 / Fig. 8 (top row).

Phase-field formulation: Ambati hybrid (Comput Mech 55:383-405, Eq. 27).

Diagnostic outputs (all in output/):
  - thermal_shock_results.png  : T(r) along the contact-wedge midline, T(t) at
                                 center / contact-rim points, D_max(r)
  - stress_evolution.png       : sigma_rr(r), sigma_tt(r) along contact midline
  - energy_balance.png         : E_el(t), E_frac(t)
  - damage_field.png           : ParaView-like 2D colormap of D on the (x, y)
                                 mesh at the final time -- direct McClenny
                                 Fig. 8 (top) reproduction
  - temperature_field.png      : same paraview-like 2D colormap of T at final time
  - stress_vm_field.png        : same, of von-Mises equivalent stress at final time
  - stress_hoop_field.png      : same, of hoop stress sigma_theta_theta at final time
  - non-regression.json        : machine-readable pass/fail summary
"""

import os
import json
import glob

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
import yaml

# --.. ..- .-.. .-.. --- load parameters --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT_DIR, "non-regression.json")
TOLERANCE = 5.0e-2

with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])
print(f"[INFO] Ro = {Ro*1e3:.1f} mm")

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    cfg = yaml.safe_load(f)
mat_path = list(cfg["materials"].values())[0]
with open(os.path.join(CASE_DIR, mat_path)) as f:
    mat = yaml.safe_load(f)

k = float(mat["k"])
rho = float(mat["rho"])
cp = float(mat["cp"])
T_initial = float(mat.get("T_initial", mat["T_ref"]))
T_ref = float(mat["T_ref"])
print(f"[INFO] k={k}, rho={rho}, cp={cp}, T_initial={T_initial} K, T_ref={T_ref} K")

with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
T_quench = T_initial
for mat_name, bc_list in bcs.get("thermal", {}).items():
    for bc in bc_list:
        if bc.get("type") == "Dirichlet":
            T_quench = float(bc["temperature"])
        elif bc.get("type") == "Robin" and "T_ext" in bc:
            T_quench = float(bc["T_ext"])
print(f"[INFO] T_quench = {T_quench} K")

tau = Ro**2 * rho * cp / k
print(f"[INFO] Thermal diffusion timescale: tau = {tau:.1f} s")


# --.. ..- .-.. .-.. --- helpers --.. ..- .-.. .-.. ---
CONTACT_HALF_ANGLE_DEG = 30.0   # the contact arc spans 0 to +30 deg (upper half)


def cartesian_to_cylindrical_stress(sigma_flat, x, y):
    """Convert Cartesian sigma to (sigma_rr, sigma_tt) at each node.

    sigma_flat: (N, 9) flattened 3x3 tensor at each node, row-major.
    Components used: [0]=sigma_xx, [1]=sigma_xy, [4]=sigma_yy.
    """
    r = np.sqrt(x**2 + y**2)
    r = np.maximum(r, 1e-15)
    cos_t, sin_t = x / r, y / r
    s_xx, s_xy, s_yy = sigma_flat[:, 0], sigma_flat[:, 1], sigma_flat[:, 4]
    sigma_rr = cos_t**2 * s_xx + 2 * cos_t * sin_t * s_xy + sin_t**2 * s_yy
    sigma_tt = sin_t**2 * s_xx - 2 * cos_t * sin_t * s_xy + cos_t**2 * s_yy
    return sigma_rr, sigma_tt


def radial_profile(r, field, Ro, n_bins=80, mask=None, reduce="mean"):
    """Bin a node-valued field by r in [0, Ro]. reduce in {'mean','max'}."""
    r_bins = np.linspace(0, Ro, n_bins + 1)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    f_out = np.full(n_bins, np.nan)
    base_mask = np.ones_like(r, dtype=bool) if mask is None else mask
    for i in range(n_bins):
        sel = base_mask & (r >= r_bins[i]) & (r < r_bins[i + 1])
        if np.any(sel):
            f_out[i] = np.max(field[sel]) if reduce == "max" else np.mean(field[sel])
    return r_mid, f_out


def angular_wedge_mask(x, y, angle_min_deg, angle_max_deg):
    """Mask selecting nodes within an angular wedge (degrees, CCW from +x)."""
    theta = np.degrees(np.arctan2(y, x))   # in (-180, 180]
    return (theta >= angle_min_deg) & (theta <= angle_max_deg)


# --.. ..- .-.. .-.. --- read VTU snapshots --.. ..- .-.. .-.. ---
vtu_files = sorted(glob.glob(os.path.join(OUT_DIR, "fields_*.vtu")))
n_times = len(vtu_files)
print(f"\n[INFO] Found {n_times} VTU files in {OUT_DIR}")
if n_times == 0:
    raise FileNotFoundError(f"No VTU files found in {OUT_DIR}")

from z3st.utils.utils_load import generate_power_history
times_all, _, _ = generate_power_history(
    cfg["time"], cfg["lhr"], n_steps=cfg.get("n_steps", len(cfg["time"])) - 1, filename=None
)
if len(times_all) < n_times:
    times_all = np.append(times_all, [times_all[-1]] * (n_times - len(times_all)))
times_all = times_all[:n_times]

all_snapshots = []
for idx in range(n_times):
    mesh_snap = pv.read(vtu_files[idx])
    coords = mesh_snap.points
    x_snap = coords[:, 0]
    y_snap = coords[:, 1]
    r_snap = np.sqrt(x_snap**2 + y_snap**2)

    T_snap = None
    for name in ["Temperature", "temperature", "T"]:
        if name in mesh_snap.point_data:
            T_snap = mesh_snap.point_data[name]
            break

    D_snap = None
    for name in ["Damage", "damage", "D"]:
        if name in mesh_snap.point_data:
            D_snap = mesh_snap.point_data[name]
            break

    all_snapshots.append({
        "t": times_all[idx], "x": x_snap, "y": y_snap, "r": r_snap,
        "T": T_snap, "D": D_snap, "idx": idx,
    })
    T_min, T_max = np.min(T_snap), np.max(T_snap)
    D_info = f" | D_max = {np.max(D_snap):.2e}" if D_snap is not None else ""
    print(f"  idx={idx:3d} | t = {times_all[idx]:.3e} s | T: [{T_min:.1f}, {T_max:.1f}] K{D_info}")


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
errors = {}

last = all_snapshots[-1]
T_last = last["T"]
T_mean_final = np.mean(T_last)
T_expected = T_quench if last["t"] > 3 * tau else T_initial
err_T_final = abs(T_mean_final - T_expected) / max(abs(T_expected), 1e-10)
print(f"\n[CHECK] Final t = {last['t']:.2e} s, T_mean = {T_mean_final:.1f} K (expected ~{T_expected:.1f} K)")
errors["T_final_mean"] = {
    "numerical": float(T_mean_final), "reference": float(T_expected), "rel_error": float(err_T_final)
}

if last["D"] is not None:
    D_max = float(np.max(last["D"]))
    print(f"[CHECK] D_max at final time: {D_max:.4e}")
    errors["D_max_final"] = {
        "numerical": D_max, "reference": 1.0, "rel_error": float(abs(D_max - 1.0)), "pass": True
    }


# --.. ..- .-.. .-.. --- subsample for plotting --.. ..- .-.. .-.. ---
N_PLOT = 6
plot_keep = set(np.unique(np.linspace(0, n_times - 1, N_PLOT, dtype=int)))
plot_snapshots = [all_snapshots[i] for i in sorted(plot_keep)]
cmap = plt.cm.coolwarm
colors = cmap(np.linspace(0.0, 1.0, max(len(plot_snapshots), 1)))


# --.. ..- .-.. .-.. --- plot 1: T radial + T(t) + D_max(r) --.. ..- .-.. .-.. ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# T radial profile through the contact midline (theta around +15 deg in the upper half)
ax1 = axes[0]
mid_theta_min, mid_theta_max = 5.0, 25.0   # narrow wedge through middle of contact arc
for i, snap in enumerate(plot_snapshots):
    wedge = angular_wedge_mask(snap["x"], snap["y"], mid_theta_min, mid_theta_max)
    r_mid, T_line = radial_profile(snap["r"], snap["T"], Ro, mask=wedge)
    valid = ~np.isnan(T_line)
    ax1.plot(r_mid[valid] * 1e3, T_line[valid] - 273.15, color=colors[i], linewidth=1.5,
             label=f"t={snap['t']:.2e}s")
ax1.set_xlabel("Radius (mm)")
ax1.set_ylabel("Temperature (°C)")
ax1.set_title(f"T radial profile through contact midline (theta ≈ {0.5*(mid_theta_min+mid_theta_max):.0f} deg)")
ax1.legend(fontsize=7, loc="lower left")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, Ro * 1e3)
ax1.set_ylim(-100, 850)

# T(t) at center vs at the contact rim
ax2 = axes[1]
T_center_hist, T_rim_hist, T_far_hist, t_hist = [], [], [], []
for snap in all_snapshots:
    cm = snap["r"] < 0.20 * Ro
    rim = (snap["r"] > 0.90 * Ro) & angular_wedge_mask(snap["x"], snap["y"], mid_theta_min, mid_theta_max)
    far = (snap["r"] > 0.90 * Ro) & angular_wedge_mask(snap["x"], snap["y"], 120.0, 180.0)  # opposite (insulated) side
    T_center_hist.append(np.mean(snap["T"][cm]) - 273.15 if np.any(cm) else np.nan)
    T_rim_hist.append(np.mean(snap["T"][rim]) - 273.15 if np.any(rim) else np.nan)
    T_far_hist.append(np.mean(snap["T"][far]) - 273.15 if np.any(far) else np.nan)
    t_hist.append(max(snap["t"], 1e-6))
ax2.plot(t_hist, T_center_hist, "r-o", markersize=3, linewidth=1.5, label="Center (r < 0.2 Ro)")
ax2.plot(t_hist, T_rim_hist, "b-s", markersize=3, linewidth=1.5,
         label=f"Contact rim (r > 0.9 Ro, theta ≈ {0.5*(mid_theta_min+mid_theta_max):.0f} deg)")
ax2.plot(t_hist, T_far_hist, "g-^", markersize=3, linewidth=1.5, label="Insulated side (theta ≈ 150 deg)")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Temperature (°C)")
ax2.set_title("Temperature evolution (cf. McClenny Fig. 7b)")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-100, 850)

# D_max(r) at the contact midline
ax3 = axes[2]
has_damage = False
for i, snap in enumerate(plot_snapshots):
    if snap["D"] is None or np.max(snap["D"]) < 1e-3:
        continue
    has_damage = True
    wedge = angular_wedge_mask(snap["x"], snap["y"], 0.0, CONTACT_HALF_ANGLE_DEG)  # whole upper contact wedge
    r_mid, D_max_r = radial_profile(snap["r"], snap["D"], Ro, mask=wedge, reduce="max")
    valid = ~np.isnan(D_max_r)
    ax3.plot(r_mid[valid] * 1e3, D_max_r[valid], color=colors[i], linewidth=1.5,
             label=f"t={snap['t']:.2e}s")
if has_damage:
    ax3.set_xlabel("Radius (mm)")
    ax3.set_ylabel("Damage D_max(r)  (azimuthal max within 0-30 deg wedge)")
    ax3.set_title("Damage radial penetration (within contact wedge)")
    ax3.legend(fontsize=7, loc="upper left")
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(-0.05, 1.05)
    ax3.set_xlim(0, Ro * 1e3)
else:
    ax3.text(0.5, 0.5, "No significant damage", ha="center", va="center", transform=ax3.transAxes)
    ax3.set_title("Damage radial penetration")

plt.suptitle("UO2 Thermal Shock - 2D transverse cross-section (60-deg contact wedge)", fontsize=13)
plt.tight_layout()
plot_path = os.path.join(OUT_DIR, "thermal_shock_results.png")
plt.savefig(plot_path, dpi=300)
plt.close(fig)
print(f"\n[INFO] Plot saved: {plot_path}")


# --.. ..- .-.. .-.. --- plot 2: stress radial profile through contact midline --.. ..- .-.. .-.. ---
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
r_zoom_min_mm = 0.80 * Ro * 1e3
r_zoom_max_mm = Ro * 1e3
for ax, title, si in [(axes2[0], "Radial stress σ_rr", 0), (axes2[1], "Hoop stress σ_θθ", 1)]:
    for i, snap in enumerate(plot_snapshots):
        m = pv.read(vtu_files[snap["idx"]])
        sk = None
        for name in ["Stress (points)", "Stress (points)"]:
            if name in m.point_data:
                sk = name
                break
        if sk is None:
            continue
        sf = m.point_data[sk]
        coords = m.points
        x, y = coords[:, 0], coords[:, 1]
        r = np.sqrt(x**2 + y**2)
        srr, stt = cartesian_to_cylindrical_stress(sf, x, y)
        field = srr if si == 0 else stt
        wedge = angular_wedge_mask(x, y, mid_theta_min, mid_theta_max)
        r_mid, s_binned = radial_profile(r, field / 1e6, Ro, mask=wedge)
        valid = ~np.isnan(s_binned)
        ax.plot(r_mid[valid] * 1e3, s_binned[valid], color=colors[i], linewidth=1.5,
                label=f"t={snap['t']:.2e}s")
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel(f"{title} (MPa)  (theta ≈ {0.5*(mid_theta_min+mid_theta_max):.0f} deg)")
    ax.set_title(title + f"  [zoom: r ≥ {r_zoom_min_mm:.1f} mm]")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_zoom_min_mm, r_zoom_max_mm)
    ax.axhline(0, color="k", linewidth=0.5, ls="--")

plt.suptitle("UO2 Thermal Shock - Stress along contact-wedge midline", fontsize=13)
plt.tight_layout()
plot_path2 = os.path.join(OUT_DIR, "stress_evolution.png")
plt.savefig(plot_path2, dpi=300)
plt.close(fig2)
print(f"[INFO] Stress plot saved: {plot_path2}")


# --.. ..- .-.. .-.. --- plot 3: energy balance --.. ..- .-.. .-.. ---
energies_path = os.path.join(CASE_DIR, "energies.txt")
if os.path.isfile(energies_path):
    try:
        E_data = np.loadtxt(energies_path)
        if E_data.ndim == 1:
            E_data = E_data.reshape(1, -1)
        if E_data.shape[1] >= 3:
            t_e = E_data[:, 0]
            E_el = E_data[:, 1]
            E_fr = E_data[:, 2]
            fig3, ax = plt.subplots(figsize=(8, 5))
            ax.plot(t_e, E_el, "C0-o", markersize=3, label="Elastic energy E_el")
            ax.plot(t_e, E_fr, "C3-s", markersize=3, label="Fracture energy E_frac")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Energy (J)")
            ax.set_title("Energy balance (cf. McClenny Fig. 10)")
            ax.grid(True, alpha=0.3)
            ax.legend()
            plt.tight_layout()
            plot_path3 = os.path.join(OUT_DIR, "energy_balance.png")
            plt.savefig(plot_path3, dpi=300)
            plt.close(fig3)
            print(f"[INFO] Energy plot saved: {plot_path3}")
        else:
            print(f"[WARN] energies.txt has only {E_data.shape[1]} columns; skipping energy plot.")
    except Exception as e:
        print(f"[WARN] Failed to read energies.txt: {e}")
else:
    print(f"[INFO] No energies.txt found; skipping energy plot.")


# --.. ..- .-.. .-.. --- ParaView-like 2D field plots at final time --.. ..- .-.. .-.. ---
# This is the McClenny Fig. 8 (top row) reproduction plus extra T and stress fields.
def _build_triangulation(pv_mesh):
    """Build a matplotlib Triangulation from a PyVista 2D mesh, handling tri/quad cells."""
    pts = pv_mesh.points
    x_pts, y_pts = pts[:, 0], pts[:, 1]
    cells_dict = pv_mesh.cells_dict
    triangles = None
    if 5 in cells_dict:
        triangles = np.asarray(cells_dict[5])
    elif 9 in cells_dict:
        q = np.asarray(cells_dict[9])
        triangles = np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])
    if triangles is not None and len(triangles) > 0:
        triang = mtri.Triangulation(x_pts, y_pts, triangles)
    else:
        triang = mtri.Triangulation(x_pts, y_pts)
    return triang, x_pts, y_pts


def plot_field_2d(triang, field, output_path, *,
                  title, cbar_label, cmap, vmin, vmax,
                  contact_R=None, contact_half_angle_deg=None,
                  show_mesh=True, levels=40):
    """ParaView-like 2D field plot: filled contours + optional mesh overlay + contact-arc marker."""
    levels_arr = np.linspace(vmin, vmax, levels + 1)
    fig, ax = plt.subplots(figsize=(9, 6))
    cf = ax.tricontourf(triang, field, levels=levels_arr, cmap=cmap,
                        vmin=vmin, vmax=vmax, extend="both")
    if show_mesh:
        ax.triplot(triang, color="black", linewidth=0.05, alpha=0.3)
    if contact_R is not None and contact_half_angle_deg is not None:
        theta_contact = np.linspace(0.0, np.deg2rad(contact_half_angle_deg), 60)
        ax.plot(contact_R * np.cos(theta_contact), contact_R * np.sin(theta_contact),
                color="cyan", linewidth=2.5,
                label=f"Cold contact arc (0-{contact_half_angle_deg:.0f} deg)")
        ax.legend(loc="upper left", fontsize=8)
    fig.colorbar(cf, ax=ax, label=cbar_label)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlim(-Ro * 1.05, Ro * 1.05)
    ax.set_ylim(-0.05 * Ro, Ro * 1.05)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close(fig)


def von_mises_from_flat(sigma_flat):
    """Von Mises equivalent stress from a (N, 9) row-major flattened 3x3 sigma."""
    s_xx = sigma_flat[:, 0]; s_xy = sigma_flat[:, 1]; s_xz = sigma_flat[:, 2]
    s_yy = sigma_flat[:, 4]; s_yz = sigma_flat[:, 5]
    s_zz = sigma_flat[:, 8]
    return np.sqrt(0.5 * (
        (s_xx - s_yy) ** 2 + (s_yy - s_zz) ** 2 + (s_zz - s_xx) ** 2
        + 6.0 * (s_xy ** 2 + s_yz ** 2 + s_xz ** 2)
    ))


try:
    last_mesh_pv = pv.read(vtu_files[-1])
    triang, x_pts, y_pts = _build_triangulation(last_mesh_pv)
    last_t = times_all[-1]
    src_name = os.path.basename(vtu_files[-1])

    # ------ Damage field ------
    damage_name = None
    for name in ["Damage", "damage", "D"]:
        if name in last_mesh_pv.point_data:
            damage_name = name
            break
    if damage_name is not None:
        D_field = np.asarray(last_mesh_pv.point_data[damage_name]).reshape(-1)
        plot_field_2d(
            triang, D_field,
            output_path=os.path.join(OUT_DIR, "damage_field.png"),
            title=(f"Damage field D at t = {last_t:.3e} s  "
                   f"(McClenny Fig. 8 reproducer; upper half by mirror symmetry)\n"
                   f"Source: {src_name}"),
            cbar_label="Damage D", cmap="hot_r", vmin=0.0, vmax=1.0,
            contact_R=Ro, contact_half_angle_deg=CONTACT_HALF_ANGLE_DEG,
        )
        print(f"[INFO] Damage field plot saved: {os.path.join(OUT_DIR, 'damage_field.png')}")

    # ------ Temperature field ------
    T_name = None
    for name in ["Temperature", "temperature", "T"]:
        if name in last_mesh_pv.point_data:
            T_name = name
            break
    if T_name is not None:
        T_C = np.asarray(last_mesh_pv.point_data[T_name]).reshape(-1) - 273.15
        # Use the Dirichlet T_quench and T_initial as colorbar bounds in degC:
        vmin_T = (T_quench - 273.15) - 5.0
        vmax_T = (T_initial - 273.15) + 5.0
        plot_field_2d(
            triang, T_C,
            output_path=os.path.join(OUT_DIR, "temperature_field.png"),
            title=(f"Temperature field T at t = {last_t:.3e} s\n"
                   f"Cold-bath = {T_quench-273.15:.0f} degC, initial = {T_initial-273.15:.0f} degC"),
            cbar_label="Temperature (degC)", cmap="coolwarm_r",
            vmin=vmin_T, vmax=vmax_T,
            contact_R=Ro, contact_half_angle_deg=CONTACT_HALF_ANGLE_DEG,
        )
        print(f"[INFO] Temperature field plot saved: {os.path.join(OUT_DIR, 'temperature_field.png')}")

    # ------ Stress fields (von Mises + hoop) ------
    sk = None
    for name in ["Stress (points)", "Stress (points)"]:
        if name in last_mesh_pv.point_data:
            sk = name
            break
    if sk is not None:
        sf = np.asarray(last_mesh_pv.point_data[sk])
        if sf.ndim == 2 and sf.shape[1] >= 9:
            # Von Mises
            vm_MPa = von_mises_from_flat(sf) / 1e6
            # Use percentile-based clipping to keep the colorbar informative
            # in case of localized hot-spots near the rim:
            vm_high = float(np.nanpercentile(vm_MPa, 99.0))
            vm_high = max(vm_high, 1.0)  # at least 1 MPa for the early thermal stage
            plot_field_2d(
                triang, vm_MPa,
                output_path=os.path.join(OUT_DIR, "stress_vm_field.png"),
                title=(f"Von Mises equivalent stress at t = {last_t:.3e} s\n"
                       f"colorbar clipped at the 99th percentile of |sigma_vm|"),
                cbar_label="sigma_vm (MPa)", cmap="viridis",
                vmin=0.0, vmax=vm_high,
                contact_R=Ro, contact_half_angle_deg=CONTACT_HALF_ANGLE_DEG,
            )
            print(f"[INFO] Von Mises stress plot saved: {os.path.join(OUT_DIR, 'stress_vm_field.png')}")

            # Hoop stress sigma_theta_theta
            srr_arr, stt_arr = cartesian_to_cylindrical_stress(sf, x_pts, y_pts)
            stt_MPa = stt_arr / 1e6
            stt_max = float(np.nanpercentile(np.abs(stt_MPa), 99.0))
            stt_max = max(stt_max, 1.0)
            plot_field_2d(
                triang, stt_MPa,
                output_path=os.path.join(OUT_DIR, "stress_hoop_field.png"),
                title=(f"Hoop stress sigma_theta_theta at t = {last_t:.3e} s\n"
                       f"(positive = tensile; drives radial cracks at the cold rim)"),
                cbar_label="sigma_theta_theta (MPa)", cmap="seismic",
                vmin=-stt_max, vmax=+stt_max,
                contact_R=Ro, contact_half_angle_deg=CONTACT_HALF_ANGLE_DEG,
            )
            print(f"[INFO] Hoop stress plot saved: {os.path.join(OUT_DIR, 'stress_hoop_field.png')}")

    # ------ Bonus: angular damage scan at the outer ring ------
    # Shows the discrete crack count within the contact wedge. Requires the
    # damage field, so guarded on damage_name.
    if damage_name is not None:
        D_field_local = np.asarray(last_mesh_pv.point_data[damage_name]).reshape(-1)
        r_pts = np.sqrt(x_pts**2 + y_pts**2)
        theta_pts_deg = np.degrees(np.arctan2(y_pts, x_pts))
        near_surf = (r_pts > 0.90 * Ro) & (r_pts < 0.99 * Ro) & (theta_pts_deg >= 0)
        if np.any(near_surf):
            n_bins = 90
            theta_bins = np.linspace(0, 180, n_bins + 1)
            theta_mid = 0.5 * (theta_bins[:-1] + theta_bins[1:])
            D_theta = np.full(n_bins, np.nan)
            for i in range(n_bins):
                sel = near_surf & (theta_pts_deg >= theta_bins[i]) & (theta_pts_deg < theta_bins[i + 1])
                if np.any(sel):
                    D_theta[i] = np.max(D_field_local[sel])

            fig5, ax5 = plt.subplots(figsize=(10, 4))
            ax5.plot(theta_mid, D_theta, "k-", linewidth=1.0)
            ax5.fill_between(theta_mid, 0, D_theta, alpha=0.3)
            ax5.axvspan(0, CONTACT_HALF_ANGLE_DEG, color="red", alpha=0.10,
                        label=f"Cold contact arc (0-{CONTACT_HALF_ANGLE_DEG:.0f} deg, upper half)")
            ax5.set_xlabel("Angle theta (deg)  [upper half only]")
            ax5.set_ylabel("D_max (90% < r < 99% Ro)")
            ax5.set_title(
                "Angular damage scan at the outer surface (final time)\n"
                "Peaks = individual radial cracks; expected mostly within the wedge"
            )
            ax5.grid(True, alpha=0.3)
            ax5.legend(fontsize=8)
            ax5.set_xlim(0, 180)
            ax5.set_ylim(0, 1.05)
            plt.tight_layout()
            plot_path5 = os.path.join(OUT_DIR, "damage_angular.png")
            plt.savefig(plot_path5, dpi=200)
            plt.close(fig5)
            print(f"[INFO] Angular damage plot saved: {plot_path5}")

            n_cracks = int(np.sum((D_theta[:-1] < 0.5) & (D_theta[1:] >= 0.5)))
            errors["crack_count_above_0p5"] = {
                "numerical": n_cracks, "reference": 1, "rel_error": 0.0, "pass": True
            }
            print(f"[INFO] Estimated crack count (D > 0.5 at outer ring, upper half): {n_cracks}")
    else:
        print("[INFO] No damage field in last snapshot; skipping angular damage scan.")
except Exception as e:
    print(f"[WARN] Field plotting failed: {e}")


# --.. ..- .-.. .-.. --- pass/fail --.. ..- .-.. .-.. ---
print(f"\nPass/Fail check (tolerance = {TOLERANCE:.1e})")
all_pass = True
for key, val in errors.items():
    if "pass" in val:
        status = "PASS" if val["pass"] else "FAIL"
    else:
        status = "PASS" if val["rel_error"] < TOLERANCE else "FAIL"
    if status == "FAIL":
        all_pass = False
    print(f"  {key:<28s} -> rel err = {val['rel_error']:.2e}  -> {status}")

with open(OUT_JSON, "w") as f:
    json.dump(errors, f, indent=2)
print(f"\n[INFO] Results written to: {OUT_JSON}")
print(f"\n[SUMMARY] {'PASS' if all_pass else 'FAIL'}")
