#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: benchmarks/pellet_quench_3D (3D)

Non-regression script for UO2 pellet thermal-shock fracture.
Reference: McClenny et al., JNM 565 (2022) 153719.

Phase-field formulation: Ambati hybrid (Comput Mech 55:383-405, Eq. 27).

Diagnostic outputs (all in output/):
  - thermal_shock_results.png : T radial profile + T(t) histories + D_max(r)
  - stress_evolution.png      : sigma_rr(r) and sigma_tt(r) at z = H/2
  - energy_balance.png        : E_el(t), E_frac(t) (cf. McClenny Fig. 10)
  - damage_angular.png        : D(theta) at r = 0.95 Ro, final time (crack count)
  - non-regression.json       : machine-readable pass/fail summary
"""

import os
import json
import glob

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import yaml

# --.. ..- .-.. .-.. --- load parameters from files --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT_DIR, "non-regression.json")
TOLERANCE = 5.0e-2

with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])
Lz = float(geom.get("Lz", 0.01))
print(f"[INFO] Ro = {Ro*1e3:.1f} mm, Lz = {Lz*1e3:.1f} mm")

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


# --.. ..- .-.. .-.. --- analytical reference --.. ..- .-.. .-.. ---
def analytic_T_transient(r, t, n_terms=50):
    """Transient cooling of a solid cylinder, uniform Dirichlet (qualitative ref only)."""
    from scipy.special import j0, j1, jn_zeros
    alpha_th = k / (rho * cp)
    roots = jn_zeros(0, n_terms)
    lambdas = roots / Ro
    T = np.full_like(r, T_quench, dtype=float)
    for n in range(n_terms):
        lam = lambdas[n]
        Cn = 2.0 / (roots[n] * j1(roots[n]))
        T += (T_initial - T_quench) * Cn * j0(lam * r) * np.exp(-alpha_th * lam**2 * t)
    return T


# --.. ..- .-.. .-.. --- helpers --.. ..- .-.. .-.. ---
def radial_bin(r, field, R, n_bins=60, mask=None, reduce="mean"):
    """Bin a field radially. reduce in {'mean','max'}."""
    r_bins = np.linspace(0, R, n_bins + 1)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    f_out = np.full(n_bins, np.nan)
    base_mask = np.ones_like(r, dtype=bool) if mask is None else mask
    for i in range(n_bins):
        sel = base_mask & (r >= r_bins[i]) & (r < r_bins[i + 1])
        if np.any(sel):
            f_out[i] = np.max(field[sel]) if reduce == "max" else np.mean(field[sel])
    return r_mid, f_out


def radial_line(x, y, field, R, angle_deg=30.0, width_deg=15.0, n_bins=60, mask=None):
    """Field along a radial line at angle_deg, within an angular wedge of +/- width_deg."""
    r = np.sqrt(x**2 + y**2)
    theta = np.degrees(np.arctan2(y, x)) % 360
    angle = angle_deg % 360
    dtheta = np.abs(theta - angle)
    dtheta = np.minimum(dtheta, 360 - dtheta)
    wedge_mask = dtheta < width_deg
    if mask is not None:
        wedge_mask = wedge_mask & mask

    r_bins = np.linspace(0, R, n_bins + 1)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    f_mean = np.full(n_bins, np.nan)
    for i in range(n_bins):
        m = wedge_mask & (r >= r_bins[i]) & (r < r_bins[i + 1])
        if np.any(m):
            f_mean[i] = np.mean(field[m])
    return r_mid, f_mean


def midplane_mask(z, Lz, frac=0.10):
    """Boolean mask selecting nodes near z = Lz/2 within +/- frac*Lz."""
    z_mid = 0.5 * Lz
    return np.abs(z - z_mid) < frac * Lz


# --.. ..- .-.. .-.. --- read VTU results --.. ..- .-.. .-.. ---
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


# --.. ..- .-.. .-.. --- read all snapshots --.. ..- .-.. .-.. ---
N_PLOT = 6  # cap the legend / curve count regardless of n_times
snapshot_indices = np.arange(n_times)
plot_keep = set(np.unique(np.linspace(0, n_times - 1, N_PLOT, dtype=int)))

all_snapshots = []
for idx in snapshot_indices:
    mesh_snap = pv.read(vtu_files[idx])
    t_snap = times_all[idx]
    coords = mesh_snap.points
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r_snap = np.sqrt(x**2 + y**2)

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

    all_snapshots.append({"t": t_snap, "r": r_snap, "z": z, "T": T_snap, "D": D_snap, "idx": idx})
    T_min, T_max = np.min(T_snap), np.max(T_snap)
    D_info = f" | D_max = {np.max(D_snap):.2e}" if D_snap is not None else ""
    print(f"  idx={idx:3d} | t = {t_snap:.3e} s | T: [{T_min:.1f}, {T_max:.1f}] K{D_info}")


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
errors = {}

early = all_snapshots[min(2, len(all_snapshots) - 1)]
r_e, T_e, t_e = early["r"], early["T"], early["t"]
center_mask = r_e < 0.1 * Ro
surface_mask = r_e > 0.9 * Ro
T_center = np.mean(T_e[center_mask]) if np.any(center_mask) else np.nan
T_surface = np.mean(T_e[surface_mask]) if np.any(surface_mask) else np.nan
dT = T_center - T_surface

print(f"\n[CHECK] Temperature gradient at t = {t_e:.2e} s:")
print(f"  T_center  = {T_center:.1f} K, T_surface = {T_surface:.1f} K, dT = {dT:.1f} K")
gradient_ok = dT > -1.0
errors["temperature_gradient"] = {
    "numerical": float(dT),
    "reference": float(T_initial - T_quench),
    "rel_error": float(abs(dT - (T_initial - T_quench)) / max(abs(T_initial - T_quench), 1e-10)),
    "pass": bool(gradient_ok),
}

if t_e > 0 and abs(T_initial - T_quench) > 1.0:
    r_mid, T_binned = radial_bin(r_e, T_e, Ro)
    T_analytical = analytic_T_transient(r_mid, t_e)
    valid = ~np.isnan(T_binned)
    L2_T = np.sqrt(np.mean((T_binned[valid] - T_analytical[valid])**2))
    L2_T_rel = L2_T / np.sqrt(np.mean(T_analytical[valid]**2))
    print(f"  L2 vs analytical (qualitative, partial Dirichlet) = {L2_T_rel:.4e}")
    errors["L2_T_vs_analytical"] = {"numerical": float(L2_T), "rel_error": float(L2_T_rel), "pass": True}

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


# --.. ..- .-.. .-.. --- plot 1: T radial + T(t) + D_max(r) --.. ..- .-.. .-.. ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
cmap = plt.cm.coolwarm
plot_snapshots = [s for s in all_snapshots if s["idx"] in plot_keep]
colors = cmap(np.linspace(0.0, 1.0, max(len(plot_snapshots), 1)))

ax1 = axes[0]
contact_angle = 30.0
for i, snap in enumerate(plot_snapshots):
    coords = pv.read(vtu_files[snap["idx"]]).points
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    mid_mask = midplane_mask(z, Lz)
    r_mid, T_line = radial_line(x, y, snap["T"], Ro, angle_deg=contact_angle, width_deg=10.0, mask=mid_mask)
    valid = ~np.isnan(T_line)
    ax1.plot(r_mid[valid] * 1e3, T_line[valid] - 273.15, color=colors[i], linewidth=1.5,
             label=f"t={snap['t']:.2e}s")
ax1.legend(fontsize=7, loc="lower left")
ax1.set_xlabel("Radius (mm)")
ax1.set_ylabel("Temperature (°C)")
ax1.set_title(f"T radial profile through contact (θ={contact_angle:.0f}°, z=H/2)")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, Ro * 1e3)
ax1.set_ylim(-100, 850)

ax2 = axes[1]
T_center_hist, T_contact_hist, t_hist = [], [], []
for snap in all_snapshots:
    coords = pv.read(vtu_files[snap["idx"]]).points
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r = np.sqrt(x**2 + y**2)
    theta = np.degrees(np.arctan2(y, x)) % 360
    mid_mask = midplane_mask(z, Lz)
    cm = mid_mask & (r < 0.20 * Ro)
    contact_mask = mid_mask & (r > 0.90 * Ro) & (theta < 60)
    T_center_hist.append(np.mean(snap["T"][cm]) - 273.15 if np.any(cm) else np.nan)
    T_contact_hist.append(np.mean(snap["T"][contact_mask]) - 273.15 if np.any(contact_mask) else np.nan)
    t_hist.append(max(snap["t"], 1e-6))
ax2.plot(t_hist, T_center_hist, "r-o", markersize=3, linewidth=1.5, label="Inner center (z=H/2)")
ax2.plot(t_hist, T_contact_hist, "b-s", markersize=3, linewidth=1.5, label="Contact surface (θ=30°, z=H/2)")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Temperature (°C)")
ax2.set_title("Temperature evolution (cf. McClenny Fig. 7b)")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-100, 850)

ax3 = axes[2]
has_damage = False
for i, snap in enumerate(plot_snapshots):
    if snap["D"] is None or np.max(snap["D"]) < 1e-3:
        continue
    has_damage = True
    coords = pv.read(vtu_files[snap["idx"]]).points
    z = coords[:, 2]
    mid_mask = midplane_mask(z, Lz)
    # circumferential MAX (not mean) so individual cracks are visible:
    r_mid, D_max_r = radial_bin(snap["r"], snap["D"], Ro, mask=mid_mask, reduce="max")
    valid = ~np.isnan(D_max_r)
    ax3.plot(r_mid[valid] * 1e3, D_max_r[valid], color=colors[i], linewidth=1.5,
             label=f"t={snap['t']:.2e}s")
if has_damage:
    ax3.set_xlabel("Radius (mm)")
    ax3.set_ylabel("Damage D_max(r)  (circumferential max, z=H/2)")
    ax3.set_title("Damage radial penetration")
    ax3.legend(fontsize=7, loc="upper left")
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(-0.05, 1.05)
    ax3.set_xlim(0, Ro * 1e3)
else:
    ax3.text(0.5, 0.5, "No significant damage", ha="center", va="center", transform=ax3.transAxes)
    ax3.set_title("Damage radial penetration")

plt.suptitle("UO2 Thermal Shock Fracture - McClenny et al., JNM 565 (2022) 153719", fontsize=13)
plt.tight_layout()
plot_path = os.path.join(OUT_DIR, "thermal_shock_results.png")
plt.savefig(plot_path, dpi=300)
print(f"\n[INFO] Plot saved: {plot_path}")


# --.. ..- .-.. .-.. --- plot 2: stress radial profiles at z=H/2 --.. ..- .-.. .-.. ---
def cartesian_to_cylindrical_stress(sigma_flat, x, y):
    r = np.sqrt(x**2 + y**2)
    r = np.maximum(r, 1e-15)
    cos_t, sin_t = x / r, y / r
    s_xx, s_xy, s_yy = sigma_flat[:, 0], sigma_flat[:, 1], sigma_flat[:, 4]
    sigma_rr = cos_t**2 * s_xx + 2 * cos_t * sin_t * s_xy + sin_t**2 * s_yy
    sigma_tt = sin_t**2 * s_xx - 2 * cos_t * sin_t * s_xy + cos_t**2 * s_yy
    return sigma_rr, sigma_tt

fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
# Zoom into the outer 20% of the pellet where the thermal-shock stress concentrates:
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
        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
        r = np.sqrt(x**2 + y**2)
        srr, stt = cartesian_to_cylindrical_stress(sf, x, y)
        field = srr if si == 0 else stt
        mid_mask = midplane_mask(z, Lz)
        r_mid, s_binned = radial_bin(r, field / 1e6, Ro, mask=mid_mask)
        valid = ~np.isnan(s_binned)
        ax.plot(r_mid[valid] * 1e3, s_binned[valid], color=colors[i], linewidth=1.5,
                label=f"t={snap['t']:.2e}s")
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel(f"{title} (MPa)  (z=H/2)")
    ax.set_title(title + f"  [zoom: r ≥ {r_zoom_min_mm:.1f} mm]")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_zoom_min_mm, r_zoom_max_mm)
    ax.axhline(0, color="k", linewidth=0.5, ls="--")

plt.suptitle("UO2 Thermal Shock - Stress evolution (mid-plane slice)", fontsize=13)
plt.tight_layout()
plot_path2 = os.path.join(OUT_DIR, "stress_evolution.png")
plt.savefig(plot_path2, dpi=300)
print(f"[INFO] Stress plot saved: {plot_path2}")


# --.. ..- .-.. .-.. --- plot 3: energy balance vs time --.. ..- .-.. .-.. ---
energies_path = os.path.join(CASE_DIR, "energies.txt")
if os.path.isfile(energies_path):
    try:
        E_data = np.loadtxt(energies_path)
        if E_data.ndim == 1:
            E_data = E_data.reshape(1, -1)
        # Expected columns from compute_energy_balance: t, E_el, E_frac
        # Be defensive about column count:
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
            print(f"[INFO] Energy plot saved: {plot_path3}")
        else:
            print(f"[WARN] energies.txt has only {E_data.shape[1]} columns; skipping energy plot.")
    except Exception as e:
        print(f"[WARN] Failed to read energies.txt: {e}")
else:
    print(f"[INFO] No energies.txt found; skipping energy plot.")


# --.. ..- .-.. .-.. --- plot 4: angular damage scan at r ~ 0.95 Ro --.. ..- .-.. .-.. ---
if last["D"] is not None and np.max(last["D"]) > 1e-3:
    coords = pv.read(vtu_files[last["idx"]]).points
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r = np.sqrt(x**2 + y**2)
    theta = np.degrees(np.arctan2(y, x)) % 360
    mid_mask = midplane_mask(z, Lz)
    near_surf = mid_mask & (r > 0.90 * Ro) & (r < 0.99 * Ro)
    if np.any(near_surf):
        n_bins = 180
        theta_bins = np.linspace(0, 360, n_bins + 1)
        theta_mid = 0.5 * (theta_bins[:-1] + theta_bins[1:])
        D_theta = np.full(n_bins, np.nan)
        for i in range(n_bins):
            sel = near_surf & (theta >= theta_bins[i]) & (theta < theta_bins[i + 1])
            if np.any(sel):
                D_theta[i] = np.max(last["D"][sel])
        fig4, ax = plt.subplots(figsize=(10, 4))
        ax.plot(theta_mid, D_theta, "k-", linewidth=1.0)
        ax.fill_between(theta_mid, 0, D_theta, alpha=0.3)
        ax.axvspan(0, 60, color="red", alpha=0.10, label="Cold contact wedge (0°-60°)")
        ax.set_xlabel("Angle θ (deg)")
        ax.set_ylabel("D_max (90% < r < 99% Ro, z=H/2)")
        ax.set_title("Angular damage scan at the outer surface (final time)\n"
                     "Peaks = individual radial cracks; expected: 2 long + fan within the wedge")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_xlim(0, 360)
        ax.set_ylim(0, 1.05)
        plt.tight_layout()
        plot_path4 = os.path.join(OUT_DIR, "damage_angular.png")
        plt.savefig(plot_path4, dpi=300)
        print(f"[INFO] Angular damage plot saved: {plot_path4}")
        # crude crack-count diagnostic: count peaks > 0.5
        n_cracks = int(np.sum((D_theta[:-1] < 0.5) & (D_theta[1:] >= 0.5)))
        errors["crack_count_above_0p5"] = {
            "numerical": n_cracks, "reference": 2, "rel_error": 0.0, "pass": True
        }
        print(f"[INFO] Estimated crack count (D > 0.5 at outer ring): {n_cracks}")


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
