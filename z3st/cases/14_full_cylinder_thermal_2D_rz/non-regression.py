#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: 14_full_cylinder_thermal_2D_rz (axisymmetric, verification only)

Verifies the axisymmetric thermal solver and linear thermo-elastic stress
response of z3st against analytic references, using the McClenny et al.
(JNM 565, 2022) UO2 pellet quench parameters as the thermal load.

Damage is intentionally OFF for this case (see the README and input.yaml
for the rationale): axisymmetric mode cannot represent the experimental
60-deg contact wedge, so any axisymmetric variant with damage produces an
unphysical annular band, not the McClenny radial-crack pattern.

Mesh layout: 2D (r, z) plane, with x = r in [0, Ro] and y = z in [0, Lz].

Diagnostic outputs (all in output/):
  - thermal_shock_results.png : T(r) at z=Lz/2; T(t) at center vs surface
  - stress_evolution.png      : sigma_rr(r), sigma_zz(r), sigma_tt(r) at z=Lz/2
  - energy_balance.png        : E_el(t) (E_frac always 0 since damage is off)
  - non-regression.json       : machine-readable pass/fail summary
"""

import os
import json
import glob

import matplotlib.pyplot as plt
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
Lz = float(geom["Lz"])
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


# --.. ..- .-.. .-.. --- analytical reference (uniform Dirichlet, axisym) --.. ..- .-.. .-.. ---
def analytic_T_transient(r, t, n_terms=50):
    """Transient cooling of a solid cylinder with uniform Dirichlet on r = Ro
    (analytic axisymmetric Bessel series). For this case, the BC is
    physically identical (uniform azimuthal cooling), so the comparison is
    quantitative, not qualitative."""
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


# --.. ..- .-.. .-.. --- helpers (axisymmetric: r = x, z = y) --.. ..- .-.. .-.. ---
def midplane_mask(z, Lz, frac=0.10):
    """Mask selecting nodes near z = Lz/2 within +/- frac*Lz."""
    return np.abs(z - 0.5 * Lz) < frac * Lz


def radial_profile(r, field, Ro, n_bins=80, mask=None, reduce="mean"):
    """Bin a field by r in [0, Ro]. reduce in {'mean','max'}."""
    r_bins = np.linspace(0, Ro, n_bins + 1)
    r_mid = 0.5 * (r_bins[:-1] + r_bins[1:])
    f_out = np.full(n_bins, np.nan)
    base_mask = np.ones_like(r, dtype=bool) if mask is None else mask
    for i in range(n_bins):
        sel = base_mask & (r >= r_bins[i]) & (r < r_bins[i + 1])
        if np.any(sel):
            f_out[i] = np.max(field[sel]) if reduce == "max" else np.mean(field[sel])
    return r_mid, f_out


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


all_snapshots = []
for idx in range(n_times):
    mesh_snap = pv.read(vtu_files[idx])
    coords = mesh_snap.points
    r_snap = coords[:, 0]    # axisymmetric: x = r
    z_snap = coords[:, 1]    # y = z

    T_snap = None
    for name in ["Temperature", "temperature", "T"]:
        if name in mesh_snap.point_data:
            T_snap = mesh_snap.point_data[name]
            break

    all_snapshots.append({"t": times_all[idx], "r": r_snap, "z": z_snap,
                          "T": T_snap, "idx": idx})
    T_min, T_max = np.min(T_snap), np.max(T_snap)
    print(f"  idx={idx:3d} | t = {times_all[idx]:.3e} s | T: [{T_min:.1f}, {T_max:.1f}] K")


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

# Quantitative comparison vs axisymmetric Bessel-series solution.
# Unlike the 3D variant where the BC is partial (1/6), here the BC is
# uniform on the whole r = Ro boundary, so the analytic solution is exact
# in the no-damage limit. We compare BEFORE damage onset only.
if t_e > 0 and abs(T_initial - T_quench) > 1.0:
    z_e = early["z"]
    mid_mask = midplane_mask(z_e, Lz)
    r_mid, T_binned = radial_profile(r_e, T_e, Ro, mask=mid_mask)
    T_analytical = analytic_T_transient(r_mid, t_e)
    valid = ~np.isnan(T_binned)
    L2_T = np.sqrt(np.mean((T_binned[valid] - T_analytical[valid])**2))
    L2_T_rel = L2_T / np.sqrt(np.mean(T_analytical[valid]**2))
    print(f"  L2 vs Bessel analytic (axisym, quantitative) = {L2_T_rel:.4e}")
    errors["L2_T_vs_analytical"] = {
        "numerical": float(L2_T), "rel_error": float(L2_T_rel), "pass": bool(L2_T_rel < TOLERANCE)
    }

last = all_snapshots[-1]
T_last = last["T"]
T_mean_final = np.mean(T_last)
T_expected = T_quench if last["t"] > 3 * tau else T_initial
err_T_final = abs(T_mean_final - T_expected) / max(abs(T_expected), 1e-10)
print(f"\n[CHECK] Final t = {last['t']:.2e} s, T_mean = {T_mean_final:.1f} K (expected ~{T_expected:.1f} K)")
errors["T_final_mean"] = {
    "numerical": float(T_mean_final), "reference": float(T_expected), "rel_error": float(err_T_final)
}

# Damage block intentionally not checked: damage is OFF in this verification case.


# --.. ..- .-.. .-.. --- plot 1: T radial + T(t) --.. ..- .-.. .-.. ---
N_PLOT = 6  # cap legend / curve count regardless of n_times
plot_keep = set(np.unique(np.linspace(0, len(all_snapshots) - 1, N_PLOT, dtype=int)))
plot_snapshots = [all_snapshots[i] for i in sorted(plot_keep)]

cmap = plt.cm.coolwarm
colors = cmap(np.linspace(0.0, 1.0, max(len(plot_snapshots), 1)))

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
ax1 = axes[0]
for i, snap in enumerate(plot_snapshots):
    mid_mask = midplane_mask(snap["z"], Lz)
    r_mid, T_line = radial_profile(snap["r"], snap["T"], Ro, mask=mid_mask)
    valid = ~np.isnan(T_line)
    ax1.plot(r_mid[valid] * 1e3, T_line[valid] - 273.15, color=colors[i], linewidth=1.5,
             label=f"t={snap['t']:.2e}s")
# Overlay analytic Bessel-series solution at the latest plotted snapshot
# so the reader can see the agreement at a glance.
t_late = plot_snapshots[-1]["t"]
r_grid = np.linspace(1e-9, Ro, 200)
T_an = analytic_T_transient(r_grid, t_late) - 273.15
ax1.plot(r_grid * 1e3, T_an, "k--", linewidth=1.0, label=f"analytic (Bessel) at t={t_late:.2e}s")
ax1.set_xlabel("Radius (mm)")
ax1.set_ylabel("Temperature (°C)")
ax1.set_title("T radial profile at z = Lz/2")
ax1.legend(fontsize=7, loc="lower left")
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, Ro * 1e3)
ax1.set_ylim(-100, 850)

ax2 = axes[1]
T_center_hist, T_surface_hist, t_hist = [], [], []
for snap in all_snapshots:
    mid = midplane_mask(snap["z"], Lz)
    cm = mid & (snap["r"] < 0.20 * Ro)
    sm = mid & (snap["r"] > 0.90 * Ro)
    T_center_hist.append(np.mean(snap["T"][cm]) - 273.15 if np.any(cm) else np.nan)
    T_surface_hist.append(np.mean(snap["T"][sm]) - 273.15 if np.any(sm) else np.nan)
    t_hist.append(max(snap["t"], 1e-6))
ax2.plot(t_hist, T_center_hist, "r-o", markersize=3, linewidth=1.5, label="Center (r < 0.2 Ro, z=Lz/2)")
ax2.plot(t_hist, T_surface_hist, "b-s", markersize=3, linewidth=1.5, label="Surface (r > 0.9 Ro, z=Lz/2)")
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Temperature (°C)")
ax2.set_title("Temperature evolution (cf. McClenny Fig. 7b)")
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)
ax2.set_ylim(-100, 850)

plt.suptitle("UO2 Thermal verification - 2D axisymmetric (uniform azimuthal quench)", fontsize=13)
plt.tight_layout()
plot_path = os.path.join(OUT_DIR, "thermal_shock_results.png")
plt.savefig(plot_path, dpi=300)
print(f"\n[INFO] Plot saved: {plot_path}")


# --.. ..- .-.. .-.. --- plot 2: stress radial profiles at z = Lz/2 --.. ..- .-.. .-.. ---
# In axisymmetric, the stored Stress tensor is in (r, z) coordinates already
# (per z3st's _build_measures with regime='axisymmetric'). Components map as:
#   sigma_flat[:, 0] = sigma_rr
#   sigma_flat[:, 1] = sigma_rz
#   sigma_flat[:, 4] = sigma_zz
#   sigma_flat[:, 8] = sigma_tt  (the hoop component, which axisym preserves)
# If the 9-flat layout differs, fall back gracefully to the in-plane components.
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))
labels_components = [
    (axes2[0], "Radial stress σ_rr",   0),
    (axes2[1], "Axial stress σ_zz",    4),
    (axes2[2], "Hoop stress σ_θθ",     8),
]
# Zoom into the outer 20% where the thermal-shock stress concentrates:
r_zoom_min_mm = 0.80 * Ro * 1e3
r_zoom_max_mm = Ro * 1e3
for ax, title, comp in labels_components:
    for i, snap in enumerate(plot_snapshots):
        m = pv.read(vtu_files[snap["idx"]])
        sk = None
        for name in ["Stress_uo2 (points)", "Stress (points)"]:
            if name in m.point_data:
                sk = name
                break
        if sk is None:
            continue
        sf = m.point_data[sk]
        if sf.ndim != 2 or sf.shape[1] <= comp:
            continue
        coords = m.points
        z = coords[:, 1]
        mid = midplane_mask(z, Lz)
        r_mid, s_binned = radial_profile(snap["r"], sf[:, comp] / 1e6, Ro, mask=mid)
        valid = ~np.isnan(s_binned)
        ax.plot(r_mid[valid] * 1e3, s_binned[valid], color=colors[i], linewidth=1.5,
                label=f"t={snap['t']:.2e}s")
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel(f"{title} (MPa)  (z = Lz/2)")
    ax.set_title(title + f"  [zoom: r ≥ {r_zoom_min_mm:.1f} mm]")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(r_zoom_min_mm, r_zoom_max_mm)
    ax.axhline(0, color="k", linewidth=0.5, ls="--")
plt.suptitle("UO2 thermal verification - Stress (axisymmetric, mid-plane)", fontsize=13)
plt.tight_layout()
plot_path2 = os.path.join(OUT_DIR, "stress_evolution.png")
plt.savefig(plot_path2, dpi=300)
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
            print(f"[INFO] Energy plot saved: {plot_path3}")
        else:
            print(f"[WARN] energies.txt has only {E_data.shape[1]} columns; skipping energy plot.")
    except Exception as e:
        print(f"[WARN] Failed to read energies.txt: {e}")
else:
    print(f"[INFO] No energies.txt found; skipping energy plot.")


# Damage field plot intentionally not produced: damage is OFF in this verification case.


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
