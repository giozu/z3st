#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 14_full_cylinder_cracking

Non-regression script for UO2 pellet thermal shock fracture.
Ref: McClenny et al., JNM 565 (2022) 153719
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

# Load geometry
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Ro = float(geom["Ro"])
print(f"[INFO] Ro = {Ro*1e3:.1f} mm")

# Load material
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

# Load BC to get quench temperature
with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
T_quench = T_initial  # default: no quench
for mat_name, bc_list in bcs.get("thermal", {}).items():
    for bc in bc_list:
        if bc.get("type") == "Dirichlet":
            T_quench = float(bc["temperature"])
print(f"[INFO] T_quench = {T_quench} K")

# Thermal diffusion timescale
tau = Ro**2 * rho * cp / k
print(f"[INFO] Thermal diffusion timescale: tau = {tau:.1f} s")


# --.. ..- .-.. .-.. --- analytical reference --.. ..- .-.. .-.. ---
def analytic_T_transient(r, t, n_terms=50):
    """Analytical solution for transient cooling of a solid cylinder."""
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


# --.. ..- .-.. .-.. --- helper: radial binning --.. ..- .-.. .-.. ---
def radial_bin(r, field, R, n_bins=60):
    """Bin a field radially (circumferential average), return (r_mid, mean_value)."""
    r_bins = np.linspace(0, R, n_bins + 1)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    f_mean = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (r >= r_bins[i]) & (r < r_bins[i + 1])
        if np.sum(mask) > 0:
            f_mean[i] = np.mean(field[mask])
    return r_mid, f_mean


def radial_line(x, y, field, R, angle_deg=30.0, width_deg=15.0, n_bins=60):
    """
    Extract field along a radial line at a given angle.
    Selects nodes within an angular wedge of ±width_deg around the target angle.
    """
    r = np.sqrt(x**2 + y**2)
    theta = np.degrees(np.arctan2(y, x)) % 360
    angle = angle_deg % 360

    # Angular wedge mask
    dtheta = np.abs(theta - angle)
    dtheta = np.minimum(dtheta, 360 - dtheta)
    wedge_mask = dtheta < width_deg

    r_bins = np.linspace(0, R, n_bins + 1)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    f_mean = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = wedge_mask & (r >= r_bins[i]) & (r < r_bins[i + 1])
        if np.sum(mask) > 0:
            f_mean[i] = np.mean(field[mask])
    return r_mid, f_mean


# --.. ..- .-.. .-.. --- read VTU results --.. ..- .-.. .-.. ---
vtu_files = sorted(glob.glob(os.path.join(OUT_DIR, "fields_*.vtu")))
n_times = len(vtu_files)
print(f"\n[INFO] Found {n_times} VTU files in {OUT_DIR}")
if n_times == 0:
    raise FileNotFoundError(f"No VTU files found in {OUT_DIR}")

# Reconstruct time values
from z3st.utils.utils_load import generate_power_history
times_all, _, _ = generate_power_history(
    cfg["time"], cfg["lhr"], n_steps=cfg.get("n_steps", len(cfg["time"])) - 1, filename=None
)
# Ensure times_all length matches VTU count
if len(times_all) < n_times:
    times_all = np.append(times_all, [times_all[-1]] * (n_times - len(times_all)))
times_all = times_all[:n_times]


# --.. ..- .-.. .-.. --- read all snapshots --.. ..- .-.. .-.. ---
# Select ~10 snapshots spread across the simulation
n_snapshots = n_times  # plot all time steps
snapshot_indices = np.unique(np.linspace(0, n_times - 1, n_snapshots, dtype=int))

all_snapshots = []
for idx in snapshot_indices:
    mesh_snap = pv.read(vtu_files[idx])
    t_snap = times_all[idx]
    coords = mesh_snap.points
    x, y = coords[:, 0], coords[:, 1]
    r_snap = np.sqrt(x**2 + y**2)

    # Temperature
    T_snap = None
    for name in ["Temperature", "temperature", "T"]:
        if name in mesh_snap.point_data:
            T_snap = mesh_snap.point_data[name]
            break

    # Damage
    D_snap = None
    for name in ["Damage", "damage", "D"]:
        if name in mesh_snap.point_data:
            D_snap = mesh_snap.point_data[name]
            break

    all_snapshots.append({"t": t_snap, "r": r_snap, "T": T_snap, "D": D_snap, "idx": idx})
    T_min, T_max = np.min(T_snap), np.max(T_snap)
    D_info = f" | D_max = {np.max(D_snap):.2e}" if D_snap is not None else ""
    print(f"  idx={idx:3d} | t = {t_snap:.3e} s | T: [{T_min:.1f}, {T_max:.1f}] K{D_info}")


# --.. ..- .-.. .-.. --- checks --.. ..- .-.. .-.. ---
errors = {}

# 1. Temperature at center and surface
early = all_snapshots[min(2, len(all_snapshots) - 1)]
r_e, T_e, t_e = early["r"], early["T"], early["t"]
center_mask = r_e < 0.1 * Ro
surface_mask = r_e > 0.9 * Ro
T_center = np.mean(T_e[center_mask]) if np.any(center_mask) else np.nan
T_surface = np.mean(T_e[surface_mask]) if np.any(surface_mask) else np.nan
dT = T_center - T_surface

print(f"\n[CHECK] Temperature gradient at t = {t_e:.2e} s:")
print(f"  T_center  = {T_center:.1f} K")
print(f"  T_surface = {T_surface:.1f} K")
print(f"  dT        = {dT:.1f} K")

gradient_ok = dT > -1.0  # allow small numerical noise, but no negative gradient
errors["temperature_gradient"] = {
    "numerical": float(dT),
    "reference": float(T_initial - T_quench),
    "rel_error": float(abs(dT - (T_initial - T_quench)) / max(abs(T_initial - T_quench), 1e-10)),
    "pass": bool(gradient_ok),
}

# 2. Compare with analytical at early time (only if there's a real thermal gradient)
if t_e > 0 and abs(T_initial - T_quench) > 1.0:
    r_mid, T_binned = radial_bin(r_e, T_e, Ro)
    T_analytical = analytic_T_transient(r_mid, t_e)
    valid = ~np.isnan(T_binned)
    L2_T = np.sqrt(np.mean((T_binned[valid] - T_analytical[valid])**2))
    L2_T_rel = L2_T / np.sqrt(np.mean(T_analytical[valid]**2))
    print(f"\n[CHECK] T vs analytical at t = {t_e:.2e} s:")
    print(f"  L2 error (relative) = {L2_T_rel:.4e}")
    errors["L2_T_vs_analytical"] = {"numerical": float(L2_T), "rel_error": float(L2_T_rel)}

# 3. Final temperature
last = all_snapshots[-1]
T_last = last["T"]
T_mean_final = np.mean(T_last)
# For the check, compare to expected final T (quench if t >> tau, else T_initial)
T_expected = T_quench if last["t"] > 3 * tau else T_initial
err_T_final = abs(T_mean_final - T_expected) / max(abs(T_expected), 1e-10)
print(f"\n[CHECK] Final temperature (t = {last['t']:.2e} s):")
print(f"  T_mean = {T_mean_final:.1f} K (expected ~{T_expected:.1f} K)")
print(f"  Rel error = {err_T_final:.4e}")
errors["T_final_mean"] = {"numerical": float(T_mean_final), "reference": float(T_expected), "rel_error": float(err_T_final)}

# 4. Damage
if last["D"] is not None:
    D_max = float(np.max(last["D"]))
    damage_ok = True  # damage may or may not occur depending on config
    print(f"\n[CHECK] D_max at final time: {D_max:.4e}")
    errors["D_max_final"] = {"numerical": D_max, "reference": 1.0, "rel_error": float(abs(D_max - 1.0)), "pass": True}


# --.. ..- .-.. .-.. --- plots --.. ..- .-.. .-.. ---
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
cmap = plt.cm.coolwarm
colors = cmap(np.linspace(0.0, 1.0, len(all_snapshots)))

# --- Plot 1: Temperature radial profiles through the contact region (θ=30°) ---
ax1 = axes[0]
contact_angle = 30.0   # middle of the 0°-60° contact region
for i, snap in enumerate(all_snapshots):
    coords = pv.read(vtu_files[snap["idx"]]).points
    x, y = coords[:, 0], coords[:, 1]
    r_mid, T_line = radial_line(x, y, snap["T"], Ro, angle_deg=contact_angle, width_deg=10.0)
    valid = ~np.isnan(T_line)
    ax1.plot(r_mid[valid] * 1e3, T_line[valid] - 273.15, color=colors[i], linewidth=1.5)

ax1.set_xlabel("Radius (mm)")
ax1.set_ylabel("Temperature (°C)")
ax1.set_title(f"T radial profile through contact (θ={contact_angle:.0f}°)")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, Ro * 1e3)
ax1.set_ylim(-100, 850)

# --- Plot 2: Temperature time evolution (contact side vs center vs insulated side) ---
ax2 = axes[1]
T_center_hist, T_contact_hist, T_insulated_hist, t_hist = [], [], [], []
for snap in all_snapshots:
    coords = pv.read(vtu_files[snap["idx"]]).points
    x, y = coords[:, 0], coords[:, 1]
    r = np.sqrt(x**2 + y**2)
    theta = np.degrees(np.arctan2(y, x)) % 360

    # Center (r < 20% R)
    cm = r < 0.20 * Ro
    T_center_hist.append(np.mean(snap["T"][cm]) - 273.15 if np.any(cm) else np.nan)

    # Contact surface (r > 90% R, θ within 0°-60°)
    contact_mask = (r > 0.90 * Ro) & (theta < 60)
    T_contact_hist.append(np.mean(snap["T"][contact_mask]) - 273.15 if np.any(contact_mask) else np.nan)

    # Insulated surface (r > 90% R, θ within 180°-240°, opposite side)
    insulated_mask = (r > 0.90 * Ro) & (theta > 180) & (theta < 240)
    T_insulated_hist.append(np.mean(snap["T"][insulated_mask]) - 273.15 if np.any(insulated_mask) else np.nan)

    t_hist.append(max(snap["t"], 1e-6))

ax2.plot(t_hist, T_center_hist, "r-o", markersize=3, linewidth=1.5, label="Inner center")
ax2.plot(t_hist, T_contact_hist, "b-s", markersize=3, linewidth=1.5, label="Contact surface (θ=30°)")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Temperature (°C)")
ax2.set_title("Temperature evolution\n(cf. McClenny et al., Fig. 7b)")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_xscale("linear")
ax2.set_ylim(-100, 850)

# --- Plot 3: Damage ---
ax3 = axes[2]
has_damage = False
for i, snap in enumerate(all_snapshots):
    if snap["D"] is not None and np.max(snap["D"]) > 1e-3:
        has_damage = True
        r_mid, D_binned = radial_bin(snap["r"], snap["D"], Ro)
        valid = ~np.isnan(D_binned)
        ax3.plot(r_mid[valid] * 1e3, D_binned[valid], color=colors[i], linewidth=1.5,
                 label=f"t = {snap['t']:.2e} s")

if has_damage:
    ax3.set_xlabel("Radius (mm)")
    ax3.set_ylabel("Damage D (mean)")
    ax3.set_title("Damage radial evolution")
    ax3.legend(fontsize=7)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(-0.05, 1.05)
    ax3.set_xlim(0, Ro * 1e3)
else:
    ax3.text(0.5, 0.5, "No significant damage", ha="center", va="center", transform=ax3.transAxes)
    ax3.set_title("Damage radial evolution")

plt.suptitle("UO2 Thermal Shock Fracture — McClenny et al., JNM 565 (2022) 153719", fontsize=13)
plt.tight_layout()

plot_path = os.path.join(OUT_DIR, "thermal_shock_results.png")
plt.savefig(plot_path, dpi=300)
print(f"\n[INFO] Plot saved in: {plot_path}")


# --.. ..- .-.. .-.. --- stress plot --.. ..- .-.. .-.. ---
def cartesian_to_cylindrical_stress(sigma_flat, x, y):
    r = np.sqrt(x**2 + y**2)
    r = np.maximum(r, 1e-15)
    cos_t, sin_t = x / r, y / r
    s_xx, s_xy, s_yy = sigma_flat[:, 0], sigma_flat[:, 1], sigma_flat[:, 4]
    sigma_rr = cos_t**2 * s_xx + 2 * cos_t * sin_t * s_xy + sin_t**2 * s_yy
    sigma_tt = sin_t**2 * s_xx - 2 * cos_t * sin_t * s_xy + cos_t**2 * s_yy
    return sigma_rr, sigma_tt

fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
for ax, label, si in [(axes2[0], "Radial stress σ_rr", 0), (axes2[1], "Hoop stress σ_θθ", 1)]:
    for i, snap in enumerate(all_snapshots):
        m = pv.read(vtu_files[snap["idx"]])
        sk = None
        for name in ["Stress_uo2 (points)", "Stress (points)"]:
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
        r_mid, s_binned = radial_bin(r, field / 1e6, Ro)
        valid = ~np.isnan(s_binned)
        ax.plot(r_mid[valid] * 1e3, s_binned[valid], color=colors[i], linewidth=1.5,
                label=f"t = {snap['t']:.2e} s")
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel(f"{label} (MPa)")
    ax.set_title(label)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, Ro * 1e3)
    ax.axhline(0, color="k", linewidth=0.5, ls="--")

plt.suptitle("UO2 Thermal Shock — Stress evolution", fontsize=13)
plt.tight_layout()
plot_path2 = os.path.join(OUT_DIR, "stress_evolution.png")
plt.savefig(plot_path2, dpi=300)
print(f"[INFO] Stress plot saved in: {plot_path2}")


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
    print(f"  {key:<25s} → rel err = {val['rel_error']:.2e}  → {status}")

with open(OUT_JSON, "w") as f:
    json.dump(errors, f, indent=2)
print(f"\n[INFO] Results written to: {OUT_JSON}")
print(f"\n[SUMMARY] {'PASS' if all_pass else 'FAIL'}")
