#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST case U_coaxial_contact_2D : PCMI penalty-contact plots
#
# Produces, into output/:
#   mesh.png          2D axisymmetric (r-z) mesh: pellet + gap + cladding
#   pcmi_curves.png   gap closure, contact pressure, and load transfer vs LHR
#   fields.png        final-step temperature, deformed shape, clad von Mises
#
# Run (with the z3st conda env active), after a solve has written output/:
#   python3 plots.py
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

# --- geometry / contact parameters (keep in sync with input.yaml) ---
R_PELLET = 0.0041      # pellet outer radius (m)
R_CLAD_I = 0.00413     # clad inner radius (m)
R_CLAD_O = 0.00475     # clad outer radius (m)
G0 = R_CLAD_I - R_PELLET       # initial gap (m)
K_PEN = 5.0e13                 # penalty stiffness (Pa/m)


def _imposed_lhr_kwm():
    """LHR (kW/m) per solve step, derived from input.yaml (any ramp/length)."""
    with open(os.path.join(HERE, "input.yaml")) as f:
        inp = yaml.safe_load(f)
    raw = inp["n_steps"]  # int: legacy total points; list: intervals per segment
    n_increments = raw if isinstance(raw, (list, tuple)) else int(raw) - 1
    _, lhrs, _ = generate_power_history(inp["time"], inp["lhr"],
                                        n_steps=n_increments, filename=None)
    return np.array(lhrs) / 1e3


LHR = _imposed_lhr_kwm()       # kW/m per step (matches the actual run)

TOL = 2.0e-5  # radial tolerance to pick a surface band


def surface_mean_ur(grid, r_target):
    r = grid.points[:, 0]
    ur = grid.point_data["Displacement"][:, 0]
    sel = np.abs(r - r_target) < TOL
    return ur[sel].mean(), ur[sel]


# ----------------------------------------------------------------------
# 0. IMPOSED POWER HISTORY (LHR vs time)
# ----------------------------------------------------------------------
def plot_history():
    with open(os.path.join(HERE, "input.yaml")) as f:
        inp = yaml.safe_load(f)
    t_points = inp["time"]
    lhr_points = inp["lhr"]
    raw = inp["n_steps"]
    # match __main__: int → total points − 1; list → intervals per segment
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
    # underlying ramp (the input control points, linearly interpolated)
    ax.plot(np.array(t_points) * tscale, np.array(lhr_points) / 1e3, "-",
            color="#C44E52", lw=1.5, alpha=0.6, label="imposed ramp")
    # the discrete solve steps actually evaluated
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
    p.add_text("UO2 pellet | gap 30 um | Zircaloy clad",
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
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
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

    # left: gap closure + contact pressure
    ax1.plot(LHR, gaps, "o-", color="#4C72B0", label="gap")
    ax1.axhline(0, color="grey", lw=0.8, ls=":")
    ax1.set_xlabel("linear heat rate (kW/m)")
    ax1.set_ylabel("gap (um)", color="#4C72B0")
    ax1.tick_params(axis="y", labelcolor="#4C72B0")
    ax2 = ax1.twinx()
    ax2.plot(LHR, pres, "s--", color="#C44E52", label="contact pressure")
    ax2.set_ylabel("contact pressure (MPa)", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")
    ax1.set_title("Gap closure and contact pressure")

    # right: load transfer (radial displacement of the two surfaces)
    ax3.plot(LHR, ur_p, "o-", color="#4C72B0", label="pellet outer")
    ax3.plot(LHR, ur_c, "s-", color="#DD8452", label="clad inner")
    ax3.set_xlabel("linear heat rate (kW/m)")
    ax3.set_ylabel("radial displacement (um)")
    ax3.set_title("Load transfer (clad pushed out on contact)")
    ax3.legend(loc="upper left")
    ax3.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "pcmi_curves.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/pcmi_curves.png")
    print(f"    final: gap={gaps[-1]:+.2f} um, p={pres[-1]:.1f} MPa, T_max={tmax[-1]:.0f} K")


# ----------------------------------------------------------------------
# 2b. RADIAL DISPLACEMENT PROFILE u_r(r) at mid-height, all power levels
# ----------------------------------------------------------------------
def plot_radial_profile():
    files = sorted(glob.glob(os.path.join(OUT, "*.vtu")))
    z_mid = 0.005          # mid-height (m)
    dz = 3.0e-4            # band half-width to pick a horizontal cut

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
        # pellet segment and clad segment, drawn separately (gap between)
        for lo, hi in [(0.0, R_PELLET + 1e-9), (R_CLAD_I - 1e-9, R_CLAD_O + 1e-9)]:
            sel = (rb >= lo) & (rb <= hi)
            o = np.argsort(rb[sel])
            ax.plot(rb[sel][o] * 1e3, ub[sel][o], "-", color=color, lw=1.8)
        ax.plot([], [], "-", color=color, lw=1.8, label=f"{LHR[i]:.0f} kW/m")

    # shade the gap region
    ax.axvspan(R_PELLET * 1e3, R_CLAD_I * 1e3, color="0.85", alpha=0.7)
    ax.text((R_PELLET + R_CLAD_I) / 2 * 1e3, ax.get_ylim()[1] * 0.05,
            "gap", ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("radial displacement u_r (um)")
    ax.set_title("Radial displacement profile (mid-height)\n"
                 "pellet expands, gap closes, clad pushed out")
    ax.legend(title="linear heat rate", fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    # the 0 kW/m curve is already non-zero: stress-free reference is T_ref
    # (293/300 K) while the rod sits in 580 K coolant -> baseline expansion.
    ax.text(0.97, 0.06, "0 kW/m curve is non-zero: thermal expansion\n"
            "from T_ref (293 K) to 580 K coolant",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right", color="#444",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "radial_profile.png"), dpi=150)
    plt.close(fig)
    print("  wrote output/radial_profile.png")


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
    # Single merged von Mises field; the clad is selected by region (the `clad`
    # extract above), so reading it on those cells gives the cladding stress.
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


if __name__ == "__main__":
    if not glob.glob(os.path.join(OUT, "*.vtu")):
        raise SystemExit("No output/*.vtu found - run the case first (Allrun).")
    print("[plots] generating figures...")
    plot_history()
    plot_mesh()
    plot_pcmi_curves()
    plot_radial_profile()
    plot_fields()
    print("[plots] done.")
