#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import glob
import os
import warnings

import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pyvista as pv

from z3st.utils.utils_load import generate_power_history

pv.OFF_SCREEN = True
# start_xvfb is deprecated in recent PyVista but is the working headless path
# here; silence the deprecation notice rather than change the render backend.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        pv.start_xvfb()
    except Exception:
        pass

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")

R_PELLET = 0.0045


def _imposed_lhr_kwm():
    """imposed LHR per step (kW/m)."""
    with open(os.path.join(HERE, "input.yaml")) as f:
        inp = yaml.safe_load(f)
    raw = inp["n_steps"]
    n_increments = raw if isinstance(raw, (list, tuple)) else int(raw) - 1
    _, lhrs, _ = generate_power_history(inp["time"], inp["lhr"],
                                        n_steps=n_increments, filename=None)
    return np.array(lhrs) / 1e3


LHR = _imposed_lhr_kwm()

# ----------------------------------------------------------------------
# 0. LHR vs time
# ----------------------------------------------------------------------
def plot_history():
    with open(os.path.join(HERE, "input.yaml")) as f:
        inp = yaml.safe_load(f)
    t_points = inp["time"]
    lhr_points = inp["lhr"]
    raw = inp["n_steps"]
    n_increments = raw if isinstance(raw, (list, tuple)) else int(raw) - 1
    times, lhrs, _ = generate_power_history(t_points, lhr_points,
                                            n_steps=n_increments, filename=None)

    # display time in hours if the ramp is long, else seconds
    span = max(t_points) - min(t_points)
    if span >= 3600:
        tscale, tunit = 1.0 / 3600.0, "h"
    else:
        tscale, tunit = 1.0, "s"

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(np.array(t_points) * tscale, np.array(lhr_points) / 1e3, "-",
            color="#C44E52", lw=1.5, alpha=0.6, label="imposed ramp")
    ax.plot(np.array(times) * tscale, np.array(lhrs) / 1e3, "o",
            color="#C44E52", ms=6, label="solve steps")
    ax.set_xlabel(f"time ({tunit})")
    ax.set_ylabel("linear heat rate (kW/m)")
    ax.set_title("Imposed power history (transient)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "power_history.png"), dpi=150)
    plt.close(fig)
    print(f"  wrote output/power_history.png  ({len(times)} steps: "
          f"{', '.join(f'{q/1e3:.0f}' for q in lhrs)} kW/m)")


# ----------------------------------------------------------------------
# 2. OPERATING HISTORY: burnup and temperature vs time
# ----------------------------------------------------------------------
def plot_operating_history():
    """Fuel-average burnup and min/max temperature vs time, from
    output/history.csv (written by diagnostics.py; format-independent)."""
    d = np.genfromtxt(os.path.join(OUT, "history.csv"), delimiter=",", names=True)
    bu = np.atleast_1d(d["burnup_avg_MWdkgU"])
    days = np.atleast_1d(d["time_days"])
    tmax = np.atleast_1d(d["T_max_K"])
    tmin = np.atleast_1d(d["T_min_K"])

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.set_title("Operating history")
    ax.plot(days, bu, "o-", color="#55A868", label="fuel-avg burnup")
    ax.set_xlabel("time (days)")
    ax.set_ylabel("burnup (MWd/kgU)", color="#55A868")
    ax.tick_params(axis="y", labelcolor="#55A868")
    ax2 = ax.twinx()
    ax2.plot(days, tmax, "^--", color="#911B1F", label="max temperature")
    ax2.plot(days, tmin, "^--", color="#BE7375", label="min temperature")
    ax2.set_ylabel("temperature (K)", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "operating_history.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/operating_history.png  (from history.csv)")
    print(f"    final: bu={bu[-1]:.1f} MWd/kgU, T_max={tmax[-1]:.0f} K")


# ----------------------------------------------------------------------
# 2b. STRESS PROFILE (radial + hoop) at mid-height, last step
# ----------------------------------------------------------------------
def _stress_profile_data():
    """(r, sigma_rr, sigma_theta) in MPa on a mid-height band, last step.

    Prefers nodal VTU stress; falls back to the cell-wise XDMF field
    (piecewise constant, sampled at cell centres). The stress tensor is
    stored in cylindrical order (r, theta, z), so in the flattened
    9-component array column 0 is sigma_rr and column 4 is sigma_theta.
    """
    z_mid, dz = 0.005, 3.0e-4
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    if files:
        g = pv.read(files[-1])
        if "Stress (points)" not in g.point_data:
            return None
        r, z = g.points[:, 0], g.points[:, 1]
        s = np.asarray(g.point_data["Stress (points)"]).reshape(-1, 9)
    else:
        # the solver keeps fields.h5 locked while running; read-only access to
        # flushed steps is safe, so disable HDF5 locking before h5py loads
        os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
        from z3st.utils.utils_extract_xdmf import extract_field_xdmf
        xdmf = os.path.join(OUT, "fields.xdmf")
        if not os.path.exists(xdmf):
            return None
        r, z, _, s = extract_field_xdmf(xdmf, "Stress", step_index=-1)
        s = np.asarray(s).reshape(len(r), 9)
    band = np.abs(z - z_mid) < dz
    return r[band], s[band, 0] / 1e6, s[band, 4] / 1e6


def plot_stress_profile():
    data = _stress_profile_data()
    if data is None:
        print("  [plots] no stress field in output: skipping stress profile")
        return
    r, s_rr, s_tt = data

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    sel = r <= R_PELLET + 1e-9
    o = np.argsort(r[sel])
    ax.plot(r[sel][o] * 1e3, s_rr[sel][o], "o-", color="#4C72B0", ms=3,
            lw=1.8, label=r"$\sigma_{rr}$ (radial)")
    ax.plot(r[sel][o] * 1e3, s_tt[sel][o], "s-", color="#C44E52", ms=3,
            lw=1.8, label=r"$\sigma_{\theta\theta}$ (hoop)")

    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("stress (MPa)")
    ax.set_title("Radial and hoop stress at mid-height (last step)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/stress_profile.png")

    # console summary at the fuel centre and surface
    def at(rt):
        i = np.argmin(np.abs(r - rt))
        return s_rr[i], s_tt[i]
    print("    mid-height stresses (MPa):")
    for name, rt in [("fuel centre", 0.0), ("fuel surface", R_PELLET)]:
        srr, stt = at(rt)
        print(f"      {name:12s} r={rt * 1e3:6.3f} mm   "
              f"sigma_rr={srr:9.1f}   sigma_theta={stt:9.1f}")


# ----------------------------------------------------------------------
# 3. FINAL-STEP FIELDS: temperature, radial displacement, clad von Mises
# ----------------------------------------------------------------------
def _bar_args(title):
    return dict(title=title, vertical=True, position_x=0.80, position_y=0.20,
                height=0.6, width=0.10, title_font_size=18, label_font_size=15,
                n_labels=4, fmt="%.0f", color="black")


def plot_fields():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    grid = pv.read(files[-1])

    # 3D displacement vector for warping (u_r, u_z, 0)
    u2 = grid.point_data["Displacement"]
    u3 = np.zeros((u2.shape[0], 3))
    u3[:, 0] = u2[:, 0]
    u3[:, 1] = u2[:, 1]
    grid.point_data["u3"] = u3
    grid.point_data["u_r_um"] = u2[:, 0] * 1e6
    grid.point_data["T"] = grid.point_data["Temperature"]

    # ---- (a) temperature ----
    p = pv.Plotter(off_screen=True, window_size=(560, 1050))
    p.add_mesh(grid, scalars="T", cmap="inferno", show_scalar_bar=True,
               scalar_bar_args=_bar_args("T (K)"))
    p.add_text("Temperature (centreline 1706 K -> coolant 580 K)",
               font_size=11, position=(0.03, 0.965), viewport=True, color="black")
    p.view_xy(); p.camera.zoom(0.92)
    p.screenshot(os.path.join(OUT, "field_temperature.png")); p.close()
    print("  wrote output/field_temperature.png")

    # ---- (b) radial displacement field u_r on the (undeformed) mesh ----
    p = pv.Plotter(off_screen=True, window_size=(560, 1050))
    p.add_mesh(grid, scalars="u_r_um", cmap="viridis", show_edges=False,
               scalar_bar_args=_bar_args("u_r (um)"))
    p.add_text("Radial displacement (thermal expansion + swelling)",
               font_size=11, position=(0.03, 0.965), viewport=True, color="black")
    p.view_xy(); p.camera.zoom(0.92)
    p.screenshot(os.path.join(OUT, "field_radial_disp.png")); p.close()
    print("  wrote output/field_radial_disp.png")


# ----------------------------------------------------------------------
# 4. FISSION GAS (SCIANTIX coupling): total Xe+Kr in each state
# ----------------------------------------------------------------------
# (vtu point_data key, legend label, colour) for the four gas states the
# coupling writes (spine.fg_fields / SciantixField.gas_concentrations).
_FG_KEYS = [
    ("FG produced (at_m3)",          "produced",       "#4C72B0"),
    ("FG in grain (at_m3)",          "in grain",       "#55A868"),
    ("FG at grain boundary (at_m3)", "grain boundary", "#C44E52"),
    ("FG released (at_m3)",          "released",       "#8172B3"),
]


def _have_fg(g):
    return "FG produced (at_m3)" in g.point_data


def plot_fg_radial():
    """Radial profiles of total fission gas (Xe+Kr) at mid-height, last step."""
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    if not files:
        return
    g = pv.read(files[-1])
    if not _have_fg(g):
        print("  [plots] no FG fields in output: skipping FG radial profile")
        return
    z_mid, dz = 0.005, 3.0e-4
    r, z = g.points[:, 0], g.points[:, 1]
    band = (np.abs(z - z_mid) < dz) & (r <= R_PELLET + 1e-9)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for key, label, color in _FG_KEYS:
        v = g.point_data[key]
        o = np.argsort(r[band])
        ax.plot(r[band][o] * 1e3, v[band][o], "-", color=color, lw=1.8, label=label)
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("fission gas Xe+Kr (at/m³)")
    ax.set_title("Radial fission-gas profile — total Xe+Kr per state\n"
                 "(mid-height, last step)", fontsize=10)
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fg_radial_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/fg_radial_profile.png")


def plot_fg_axial():
    """Axial profiles of total fission gas at the centreline and the pellet surface,
    last step."""
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    if not files:
        return
    g = pv.read(files[-1])
    if not _have_fg(g):
        print("  [plots] no FG fields in output: skipping FG axial profile")
        return
    r, z = g.points[:, 0], g.points[:, 1]
    fuel = r <= R_PELLET + 1e-9
    rf = r[fuel]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    for ax, r_target, rname in [(axes[0], 0.0, "centreline"),
                                (axes[1], R_PELLET, "pellet surface")]:
        # snap to the single nearest radial node-column: a fixed-r band would
        # straddle two columns where the rim peak varies steeply (sawtooth).
        r_near = rf[np.argmin(np.abs(rf - r_target))]
        sel = fuel & (np.abs(r - r_near) < 1.0e-6)
        for key, label, color in _FG_KEYS:
            v = g.point_data[key]
            o = np.argsort(z[sel])
            ax.plot(z[sel][o] * 1e3, v[sel][o], "-", color=color, lw=1.6, label=label)
        ax.set_xlabel("axial position z (mm)")
        ax.set_title(f"{rname}  (r = {r_near * 1e3:.2f} mm)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("fission gas Xe+Kr (at/m³)")
    axes[0].legend(fontsize=8)
    fig.suptitle("Axial fission-gas profile (last step) — "
                 "near-flat under the uniform axial power")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fg_axial_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/fg_axial_profile.png")


def plot_fg_average():
    """Fuel-average fission-gas inventory and fractional release vs burnup, from
    output/history.csv (written by diagnostics.py; format-independent)."""
    csv = os.path.join(OUT, "history.csv")
    if not os.path.exists(csv):
        return
    d = np.genfromtxt(csv, delimiter=",", names=True)
    if d.dtype.names is None or "fg_produced_avg" not in d.dtype.names:
        print("  [plots] history.csv has no FG columns: skipping FG average plot")
        return
    bu = np.atleast_1d(d["burnup_avg_MWdkgU"])
    fgr = np.atleast_1d(d["fgr_frac"]) * 100.0
    series = [
        ("fg_produced_avg", "produced", "#4C72B0"),
        ("fg_ingrain_avg", "in grain", "#55A868"),
        ("fg_gb_avg", "grain boundary", "#C44E52"),
        ("fg_released_avg", "released", "#8172B3"),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for col, label, color in series:
        ax1.plot(bu, np.atleast_1d(d[col]), "o-", color=color, ms=3, label=label)
    ax1.set_xlabel("fuel-average burnup (MWd/kgU)")
    ax1.set_ylabel("fuel-average fission gas Xe+Kr (at/m³)")
    ax1.set_title("Fission-gas inventory vs burnup")
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.plot(bu, fgr, "o-", color="#8172B3", ms=3)
    ax2.set_xlabel("fuel-average burnup (MWd/kgU)")
    ax2.set_ylabel("fractional fission-gas release (%)")
    ax2.set_title("Fractional fission-gas release vs burnup")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fg_average.png"), dpi=150)
    plt.close(fig)
    print(f"  wrote output/fg_average.png  (final FGR = {fgr[-1]:.2f} %, "
          f"released = {np.atleast_1d(d['fg_released_avg'])[-1]:.3e} at/m³)")


if __name__ == "__main__":
    have_vtu = bool(glob.glob(os.path.join(OUT, "*.vtu")))
    have_csv = os.path.exists(os.path.join(OUT, "history.csv"))
    if not have_vtu and not have_csv:
        raise SystemExit("No output found - run the case first (Allrun).")
    print("[plots] generating figures...")
    plot_history()
    plot_operating_history()     # CSV-preferred; works for XDMF-only runs too
    plot_stress_profile()        # VTU-preferred; falls back to the XDMF h5
    plot_fg_average()            # from history.csv (works regardless of field format)
    if have_vtu:
        plot_fields()
        plot_fg_radial()
        plot_fg_axial()
    else:
        print("  [plots] XDMF-only run: skipping field renders (pyvista cannot "
              "read dolfinx XDMF). Open output/fields.xdmf in ParaView, or re-run "
              "with 'output: {format: vtu}' for those figures.")
    print("[plots] done.")
