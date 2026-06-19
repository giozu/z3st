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

# parameters
R_PELLET = 0.0045
R_CLAD_I = 0.004565
R_CLAD_O = 0.005315
G0 = R_CLAD_I - R_PELLET
K_PEN = 5.0e13


def burnup_per_step(files):
    """Average burnup (MWd/kgU), from the Burnup VTU field"""
    bu = []
    for f in files:
        g = pv.read(f)
        if "Burnup" not in g.point_data:
            bu.append(0.0)
            continue
        r = g.points[:, 0]
        b = g.point_data["Burnup"]
        fuel = r <= R_PELLET + 1e-6
        bu.append(float(b[fuel].mean()))
    return np.array(bu)


def _imposed_lhr_kwm():
    """LHR (kW/m) """
    with open(os.path.join(HERE, "input.yaml")) as f:
        inp = yaml.safe_load(f)
    raw = inp["n_steps"]
    n_increments = raw if isinstance(raw, (list, tuple)) else int(raw) - 1
    _, lhrs, _ = generate_power_history(inp["time"], inp["lhr"],
                                        n_steps=n_increments, filename=None)
    return np.array(lhrs) / 1e3


LHR = _imposed_lhr_kwm()

TOL = 2.0e-5


def surface_mean_ur(grid, r_target):
    r = grid.points[:, 0]
    ur = grid.point_data["Displacement"][:, 0]
    sel = np.abs(r - r_target) < TOL
    return ur[sel].mean(), ur[sel]


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
# 1. MESH PLOT
# ----------------------------------------------------------------------
def plot_mesh():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    grid = pv.read(files[0])
    # region by cell-centre radius: pellet vs cladding (gap is unmeshed)
    cc = grid.cell_centers().points[:, 0]
    region = np.where(cc <= R_PELLET + 1e-6, 0, 1)  # 0 pellet, 1 clad
    grid.cell_data["region"] = region

    p = pv.Plotter(off_screen=True, window_size=(700, 1100))
    p.add_mesh(grid, scalars="region", show_edges=True, line_width=1,
               cmap=["#4C72B0", "#DD8452"], show_scalar_bar=False)
    p.add_text("UO2 pellet | gap 65 um | Zircaloy clad",
               font_size=12, position=(0.04, 0.965), viewport=True, color="black")
    p.view_xy()
    p.camera.zoom(0.95)
    p.add_axes(xlabel="r", ylabel="z")
    p.screenshot(os.path.join(OUT, "mesh.png"))
    p.close()
    print("  wrote output/mesh.png")


# ----------------------------------------------------------------------
# 2. PCMI CURVES (gap, contact pressure, load transfer) vs LHR
# ----------------------------------------------------------------------
def plot_pcmi_curves():
    """Prefers output/history.csv (written by diagnostics.py)
    It works for VTU *or* XDMF runs. Falls back to per-step VTU reads when no CSV
    is present."""
    if os.path.exists(os.path.join(OUT, "history.csv")):
        _pcmi_from_csv()
    else:
        _pcmi_from_vtu()


def _pcmi_from_csv():
    d = np.genfromtxt(os.path.join(OUT, "history.csv"), delimiter=",", names=True)
    bu = np.atleast_1d(d["burnup_avg_MWdkgU"])
    days = np.atleast_1d(d["time_days"])
    gaps = np.atleast_1d(d["gap_um"])
    pres = np.atleast_1d(d["contact_pressure_MPa"])
    tmax = np.atleast_1d(d["T_max_K"])
    tmin = np.atleast_1d(d["T_min_K"])

    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(12, 4.6))

    ax1.set_title("Burnup-driven gap closure and PCMI contact pressure")
    ax1.plot(bu, gaps, "o-", color="#4C72B0", label="gap")
    ax1.axhline(0, color="grey", lw=0.8, ls=":")
    ax1.set_xlabel("fuel-average burnup (MWd/kgU)")
    ax1.set_ylabel("gap (um)", color="#4C72B0")
    ax1.tick_params(axis="y", labelcolor="#4C72B0")
    ax2 = ax1.twinx()
    ax2.plot(bu, pres, "s--", color="#C44E52", label="contact pressure")
    ax2.set_ylabel("contact pressure (MPa)", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")

    ax3.set_title("Operating history")
    ax3.plot(days, bu, "o-", color="#55A868", label="fuel-avg burnup")
    ax3.set_xlabel("time (days)")
    ax3.set_ylabel("burnup (MWd/kgU)", color="#55A868")
    ax3.tick_params(axis="y", labelcolor="#55A868")
    ax4 = ax3.twinx()
    ax4.plot(days, tmax, "^--", color="#911B1F", label="max temperature")
    ax4.plot(days, tmin, "^--", color="#BE7375", label="min temperature")
    ax4.set_ylabel("temperature (K)", color="#C44E52")
    ax4.tick_params(axis="y", labelcolor="#C44E52")
    ax3.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pcmi_curves.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/pcmi_curves.png  (from history.csv)")
    onset = next((bu[i] for i in range(len(bu)) if pres[i] > 0.0), None)
    onset_s = f"{onset:.1f} MWd/kgU" if onset is not None else "not reached"
    print(f"    PCMI onset at burnup = {onset_s}")
    print(f"    final: bu={bu[-1]:.1f} MWd/kgU, gap={gaps[-1]:+.2f} um, "
          f"p={pres[-1]:.1f} MPa, T_max={tmax[-1]:.0f} K")


def _pcmi_from_vtu():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    bu = burnup_per_step(files)        # rod-average burnup (MWd/kgU), x-axis
    gaps, pres, ur_p, ur_c, tmax = [], [], [], [], []
    for f in files:
        g = pv.read(f)
        urp, _ = surface_mean_ur(g, R_PELLET)
        urc, _ = surface_mean_ur(g, R_CLAD_I)
        gap = G0 + urc - urp
        gaps.append(gap * 1e6)             # um
        pres.append(K_PEN * max(0.0, -gap) / 1e6)  # MPa
        ur_p.append(urp * 1e6)
        ur_c.append(urc * 1e6)
        tmax.append(float(g.point_data["Temperature"].max()))

    fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(12, 4.6))

    # left: gap closure + contact pressure vs burnup
    ax1.plot(bu, gaps, "o-", color="#4C72B0", label="gap")
    ax1.axhline(0, color="grey", lw=0.8, ls=":")
    ax1.set_xlabel("rod-average burnup (MWd/kgU)")
    ax1.set_ylabel("gap (um)", color="#4C72B0")
    ax1.tick_params(axis="y", labelcolor="#4C72B0")
    ax2 = ax1.twinx()
    ax2.plot(bu, pres, "s--", color="#C44E52", label="contact pressure")
    ax2.set_ylabel("contact pressure (MPa)", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")
    ax1.set_title("Burnup-driven gap closure and PCMI contact pressure")

    # right: load transfer (radial displacement of the two surfaces) vs burnup
    ax3.plot(bu, ur_p, "o-", color="#4C72B0", label="pellet outer (swelling + thermal)")
    ax3.plot(bu, ur_c, "s-", color="#DD8452", label="clad inner (pushed out on contact)")
    ax3.set_xlabel("rod-average burnup (MWd/kgU)")
    ax3.set_ylabel("radial displacement (um)")
    ax3.set_title("Load transfer (clad pushed out on contact)")
    ax3.legend(loc="upper left")
    ax3.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pcmi_curves.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/pcmi_curves.png")
    # report PCMI onset (first step where the gap closes)
    onset = next((bu[i] for i in range(len(bu)) if pres[i] > 0.0), None)
    onset_str = f"{onset:.1f} MWd/kgU" if onset is not None else "not reached"
    print(f"    PCMI onset at burnup = {onset_str}")
    print(f"    final: bu={bu[-1]:.1f} MWd/kgU, gap={gaps[-1]:+.2f} um, "
          f"p={pres[-1]:.1f} MPa, T_max={tmax[-1]:.0f} K")


# ----------------------------------------------------------------------
# 2b. RADIAL PROFILE at mid-height
# ----------------------------------------------------------------------
def plot_radial_profile():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    z_mid = 0.005          # mid-height (m)
    dz = 3.0e-4            # band half-width to pick a horizontal cut

    bu = burnup_per_step(files)        # MWd/kgU per step, for the legend
    n = len(files)
    label_every = max(1, n // 6)       # 6 labelled curves, rest unlabelled
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    cmap = plt.cm.viridis
    for i, f in enumerate(files):
        g = pv.read(f)
        r = g.points[:, 0]
        z = g.points[:, 1]
        ur = g.point_data["Displacement"][:, 0] * 1e6  # um
        band = np.abs(z - z_mid) < dz
        rb, ub = r[band], ur[band]
        color = cmap(i / (len(files) - 1))
        for lo, hi in [(0.0, R_PELLET + 1e-9), (R_CLAD_I - 1e-9, R_CLAD_O + 1e-9)]:
            sel = (rb >= lo) & (rb <= hi)
            o = np.argsort(rb[sel])
            ax.plot(rb[sel][o] * 1e3, ub[sel][o], "-", color=color, lw=1.8)
        if i % label_every == 0 or i == n - 1:
            ax.plot([], [], "-", color=color, lw=1.8, label=f"{bu[i]:.0f} MWd/kgU")

    # shade the gap region
    ax.axvspan(R_PELLET * 1e3, R_CLAD_I * 1e3, color="0.85", alpha=0.7)
    ax.text((R_PELLET + R_CLAD_I) / 2 * 1e3, ax.get_ylim()[1] * 0.05,
            "gap", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("radial displacement u_r (um)")
    ax.set_title("Radial displacement profile (mid-height)\n"
                 "pellet swells with burnup, gap closes, clad pushed out")
    ax.legend(title="burnup", fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    ax.text(0.97, 0.06, "fresh (0 MWd/kgU) curve is non-zero: thermal\n"
            "expansion from T_ref (293 K) to 580 K coolant",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right", color="#444",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "radial_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/radial_profile.png")


# ----------------------------------------------------------------------
# 2c. STRESS PROFILE (radial + hoop) at mid-height, last step
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
    first = True
    for lo, hi in [(0.0, R_PELLET + 1e-9), (R_CLAD_I - 1e-9, R_CLAD_O + 1e-9)]:
        sel = (r >= lo) & (r <= hi)
        o = np.argsort(r[sel])
        ax.plot(r[sel][o] * 1e3, s_rr[sel][o], "o-", color="#4C72B0", ms=3,
                lw=1.8, label=r"$\sigma_{rr}$ (radial)" if first else None)
        ax.plot(r[sel][o] * 1e3, s_tt[sel][o], "s-", color="#C44E52", ms=3,
                lw=1.8, label=r"$\sigma_{\theta\theta}$ (hoop)" if first else None)
        first = False

    ax.axvspan(R_PELLET * 1e3, R_CLAD_I * 1e3, color="0.85", alpha=0.7)
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.text((R_PELLET + R_CLAD_I) / 2 * 1e3, ax.get_ylim()[1] * 0.05,
            "gap", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("stress (MPa)")
    ax.set_title("Radial and hoop stress at mid-height (last step)")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "stress_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/stress_profile.png")

    # console summary at the four radii of interest
    def at(rt):
        i = np.argmin(np.abs(r - rt))
        return s_rr[i], s_tt[i]
    print("    mid-height stresses (MPa):")
    for name, rt in [("fuel centre", 0.0), ("fuel surface", R_PELLET),
                     ("clad inner", R_CLAD_I), ("clad outer", R_CLAD_O)]:
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

    # cell-centre radius to split pellet / clad robustly
    cc = grid.cell_centers().points[:, 0]
    pellet = grid.extract_cells(np.where(cc <= R_PELLET + 1e-6)[0])
    clad = grid.extract_cells(np.where(cc > R_CLAD_I - 1e-6)[0])

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
    # NB: a warped mesh is misleading for contact -- each body warps by its own
    # displacement, so the pellet (large u_r) visually overruns the clad. The
    # quantitative gap-closure story is told by plot_radial_profile() instead.
    p = pv.Plotter(off_screen=True, window_size=(560, 1050))
    p.add_mesh(grid, scalars="u_r_um", cmap="viridis", show_edges=False,
               scalar_bar_args=_bar_args("u_r (um)"))
    p.add_text("Radial displacement (pellet expands into the clad)",
               font_size=11, position=(0.03, 0.965), viewport=True, color="black")
    p.view_xy(); p.camera.zoom(0.92)
    p.screenshot(os.path.join(OUT, "field_radial_disp.png")); p.close()
    print("  wrote output/field_radial_disp.png")

    # ---- (c) cladding von Mises ----
    vm_key = "VonMises (points)" if "VonMises (points)" in grid.point_data else None
    if vm_key is not None:
        vm = np.asarray(clad.point_data[vm_key]) / 1e6  # MPa
        clad.point_data["vM_MPa"] = vm
        finite = vm[np.isfinite(vm)]
        clim = [float(finite.min()), float(finite.max())]
        p = pv.Plotter(off_screen=True, window_size=(560, 1050))
        p.add_mesh(clad, scalars="vM_MPa", cmap="plasma", clim=clim,
                   show_edges=True, line_width=0.4,
                   scalar_bar_args=_bar_args("vM (MPa)"))
        p.add_text(f"Cladding von Mises ({clim[0]:.0f}-{clim[1]:.0f} MPa)",
                   font_size=11, position=(0.03, 0.965), viewport=True, color="black")
        p.view_xy(); p.camera.zoom(0.92)
        p.screenshot(os.path.join(OUT, "field_clad_vonmises.png")); p.close()
        print("  wrote output/field_clad_vonmises.png")


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
    plot_pcmi_curves()           # CSV-preferred; works for XDMF-only runs too
    plot_stress_profile()        # VTU-preferred; falls back to the XDMF h5
    plot_fg_average()            # from history.csv (works regardless of field format)
    if have_vtu:
        plot_mesh()
        plot_radial_profile()
        plot_fields()
        plot_fg_radial()
        plot_fg_axial()
    else:
        print("  [plots] XDMF-only run: skipping mesh / radial-profile / field "
              "renders (pyvista cannot read dolfinx XDMF). Open output/fields.xdmf "
              "in ParaView, or re-run with 'output: {format: vtu}' for those figures.")
    print("[plots] done.")
