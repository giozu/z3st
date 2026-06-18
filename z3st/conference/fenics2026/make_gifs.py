#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Presentation GIFs for the FEniCS 2026 talk.

Two animations:

  1. pcmi_evolution.gif   (regression/pwr_rod_2D)
     fuel + cladding warped by displacement (TRUE scale, factor 1.0), coloured
     by temperature, beside the three tracked PCMI quantities -- gap width,
     contact pressure, temperature -- with a marker sweeping the time axis so
     the mesh and the curves advance together.

     NOTE on the warp factor: the unmeshed 65 um pellet-clad gap is baked into
     the mesh coordinates and is NOT scaled by warp_by_vector (only the
     displacement is). Exaggerating the warp therefore makes the pellet
     visually overrun the clad by ~(factor-1)x65 um, which is unphysical. The
     faithful animation is factor = 1.0; the quantitative gap-closure / contact
     story is carried by the side-panel curves.

  2. damage_evolution.gif (benchmarks/damage/pellet_quench_2D_xy)
     AT1 phase-field damage field over the thermal-shock transient, cracks
     nucleating and propagating from the quenched surface inward.

Outputs land in the oral deck's figures/ directory, both as a .gif (web /
README) and as a numbered PNG frame sequence under figures/anim_<name>/ for
LaTeX \animategraphics in slides.tex.

Run with the project env:
    export PATH="/home/giovanni/miniconda3/envs/z3st/bin:$PATH"
    export HDF5_USE_FILE_LOCKING=FALSE
    python3 make_gifs.py
"""

import glob
import os
import warnings

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image
import pyvista as pv

pv.OFF_SCREEN = True
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        pv.start_xvfb()
    except Exception:
        pass

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.normpath(os.path.join(HERE, "..", "..", "cases"))
FIG = os.path.join(HERE, "oral", "figures")

PCMI_OUT = os.path.join(CASES, "regression", "pwr_rod_2D", "output")
QUENCH_OUT = os.path.join(CASES, "benchmarks", "damage", "pellet_quench_2D_xy", "output")


def _fig_to_array(fig):
    """Render a Matplotlib figure to an (H, W, 3) uint8 array."""
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return buf[:, :, :3].copy()


def _save_gif(frames, path, duration_ms, hold_last=8):
    """Write a list of RGB arrays as a looping GIF; hold the final frame."""
    frames = list(frames) + [frames[-1]] * hold_last
    imgs = [Image.fromarray(a) for a in frames]
    imgs[0].save(path, save_all=True, append_images=imgs[1:],
                 duration=duration_ms, loop=0, optimize=True)
    print(f"  wrote {os.path.relpath(path, HERE)}  ({len(imgs)} frames)")


def _dump_frames(frames, name, stride=1):
    """Write numbered PNG frames for LaTeX \\animategraphics."""
    d = os.path.join(FIG, f"anim_{name}")
    os.makedirs(d, exist_ok=True)
    for f in glob.glob(os.path.join(d, "frame_*.png")):
        os.remove(f)
    # unpadded names: LaTeX \animategraphics builds filenames as frame_<n>.png
    kept = frames[::stride]
    for i, a in enumerate(kept):
        Image.fromarray(a).save(os.path.join(d, f"frame_{i}.png"))
    print(f"  wrote {len(kept)} frames -> {os.path.relpath(d, HERE)}/frame_<n>.png  "
          f"(0..{len(kept) - 1})")
    return len(kept)


# ======================================================================
# 1. PCMI evolution
# ======================================================================
def _pcmi_grid(h5):
    """Quad UnstructuredGrid from the dolfinx-written mesh (identity node
    ordering matches VTK_QUAD here)."""
    geom = np.array(h5["Mesh/mesh/geometry"])
    topo = np.array(h5["Mesh/mesh/topology"])
    ncell = topo.shape[0]
    cells = np.hstack([np.full((ncell, 1), 4), topo]).astype(np.int64).ravel()
    ctypes = np.full(ncell, pv.CellType.QUAD, np.uint8)
    return pv.UnstructuredGrid(cells, ctypes, geom.copy()), geom


def _step_key(s):
    try:
        return (0, float(s.replace("_", ".")))
    except ValueError:
        return (1, s)


def make_pcmi_gif():
    h5path = os.path.join(PCMI_OUT, "fields.h5")
    csvpath = os.path.join(PCMI_OUT, "history.csv")
    if not (os.path.exists(h5path) and os.path.exists(csvpath)):
        print("  [pcmi] missing fields.h5 or history.csv -- run the case first")
        return

    d = np.genfromtxt(csvpath, delimiter=",", names=True)
    days = np.atleast_1d(d["time_days"])
    gap = np.atleast_1d(d["gap_um"])
    pres = np.atleast_1d(d["contact_pressure_MPa"])
    tmax = np.atleast_1d(d["T_max_K"])
    tmin = np.atleast_1d(d["T_min_K"])
    contact_day = next((days[i] for i in range(len(gap)) if gap[i] <= 0.0), None)

    h5 = h5py.File(h5path, "r")
    grid, geom = _pcmi_grid(h5)
    steps = sorted(h5["Function/Displacement"].keys(), key=_step_key)
    n = min(len(steps), len(days))

    # global temperature colour limits across the transient
    tlo, thi = float(tmin.min()), float(tmax.max())

    frames = []
    for i in range(n):
        disp = np.array(h5["Function/Displacement"][steps[i]])
        u3 = np.zeros((geom.shape[0], 3))
        u3[:, 0], u3[:, 1] = disp[:, 0], disp[:, 1]
        grid.point_data["u"] = u3
        grid.point_data["T"] = np.array(h5["Function/Temperature"][steps[i]])
        warped = grid.warp_by_vector("u", factor=1.0)   # TRUE scale -- see header

        pl = pv.Plotter(off_screen=True, window_size=(470, 820))
        pl.add_mesh(warped, scalars="T", cmap="inferno", clim=[tlo, thi],
                    show_edges=True, line_width=0.25, edge_color="dimgray",
                    show_scalar_bar=True,
                    scalar_bar_args=dict(title="T (K)", vertical=True,
                                         position_x=0.88, position_y=0.18,
                                         height=0.62, width=0.05,
                                         title_font_size=18, label_font_size=15,
                                         n_labels=4, fmt="%.0f", color="black"))
        pl.add_text("UO2 pellet  |  gap  |  Zr clad",
                    font_size=10, position=(0.02, 0.02), viewport=True,
                    color="black")
        pl.view_xy()
        pl.camera.zoom(1.25)
        mesh_img = pl.screenshot(return_img=True)
        pl.close()

        fig = plt.figure(figsize=(11.0, 5.2), dpi=110)
        gs = GridSpec(3, 2, width_ratios=[1.05, 1.55], hspace=0.12, wspace=0.18,
                      left=0.02, right=0.96, top=0.90, bottom=0.11)
        axm = fig.add_subplot(gs[:, 0])
        axm.imshow(mesh_img)
        axm.axis("off")
        axm.set_title("fuel-clad deformation (true displacement)", fontsize=11)

        cur_d = days[i]
        common = dict(lw=1.6)
        ax1 = fig.add_subplot(gs[0, 1])
        ax1.plot(days, gap, color="#4C72B0", **common)
        ax1.axhline(0, color="grey", lw=0.8, ls=":")
        ax1.scatter([cur_d], [gap[i]], color="#4C72B0", zorder=5, s=36)
        ax1.set_ylabel("gap (um)", color="#4C72B0")
        ax1.tick_params(labelbottom=False)
        ax1.grid(alpha=0.3)

        ax2 = fig.add_subplot(gs[1, 1], sharex=ax1)
        ax2.plot(days, pres, color="#C44E52", **common)
        ax2.scatter([cur_d], [pres[i]], color="#C44E52", zorder=5, s=36)
        ax2.set_ylabel("contact p (MPa)", color="#C44E52")
        ax2.tick_params(labelbottom=False)
        ax2.grid(alpha=0.3)

        ax3 = fig.add_subplot(gs[2, 1], sharex=ax1)
        ax3.plot(days, tmax, color="#911B1F", **common, label="T max")
        ax3.plot(days, tmin, color="#DD8452", **common, label="T min")
        ax3.scatter([cur_d, cur_d], [tmax[i], tmin[i]], color="black",
                    zorder=5, s=24)
        ax3.set_ylabel("T (K)")
        ax3.set_xlabel("time (days)")
        ax3.legend(fontsize=8, loc="center right")
        ax3.grid(alpha=0.3)

        for ax in (ax1, ax2, ax3):
            ax.axvline(cur_d, color="0.6", lw=0.8, ls="--")
            if contact_day is not None:
                ax.axvline(contact_day, color="green", lw=0.8, alpha=0.5)

        state = "CONTACT" if gap[i] <= 0.0 else "open gap"
        fig.suptitle(f"PCMI in a PWR rod  -  t = {cur_d:6.1f} d   "
                     f"|   gap = {gap[i]:+5.1f} um   |   p = {pres[i]:4.1f} MPa   "
                     f"|   [{state}]", fontsize=12.5)
        frames.append(_fig_to_array(fig))
        plt.close(fig)

    h5.close()
    _save_gif(frames, os.path.join(FIG, "pcmi_evolution.gif"), duration_ms=140)
    _dump_frames(frames, "pcmi")
    onset = f"{contact_day:.0f} d" if contact_day is not None else "not reached"
    print(f"    PCMI contact onset at t = {onset}; "
          f"final gap {gap[-1]:+.1f} um, p {pres[-1]:.1f} MPa")


# ======================================================================
# 2. Damage evolution
# ======================================================================
def make_damage_gif():
    files = sorted(glob.glob(os.path.join(QUENCH_OUT, "*.vtu")))
    if not files:
        print("  [damage] no VTU output in pellet_quench_2D_xy -- run the case first")
        return

    t = np.linspace(1.0e-4, 0.1, len(files))   # transient time (s), per input.yaml

    frames = []
    stride = 2                                  # 100 steps -> 50 frames
    idx = list(range(0, len(files), stride))
    if idx[-1] != len(files) - 1:
        idx.append(len(files) - 1)

    for i in idx:
        g = pv.read(files[i])
        dmg = "Damage" if "Damage" in g.point_data else None
        scalars = dmg if dmg else list(g.point_data.keys())[0]

        pl = pv.Plotter(off_screen=True, window_size=(900, 520))
        pl.add_mesh(g, scalars=scalars, cmap="inferno", clim=[0.0, 1.0],
                    show_edges=False, show_scalar_bar=True,
                    scalar_bar_args=dict(title="damage  d", vertical=True,
                                         position_x=0.88, position_y=0.20,
                                         height=0.6, width=0.04,
                                         title_font_size=22, label_font_size=18,
                                         n_labels=3, fmt="%.1f", color="black"))
        pl.add_text(f"AT1 phase-field damage  -  thermal shock\n"
                    f"t = {t[i] * 1e3:6.2f} ms   (step {i + 1}/{len(files)})",
                    font_size=13, position=(0.02, 0.90), viewport=True,
                    color="black")
        pl.view_xy()
        pl.camera.zoom(1.3)
        frames.append(pl.screenshot(return_img=True))
        pl.close()

    _save_gif(frames, os.path.join(FIG, "damage_evolution.gif"), duration_ms=110)
    _dump_frames(frames, "damage")
    print(f"    {len(frames)} frames from {len(files)} steps")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    print("[gifs] PCMI evolution ...")
    make_pcmi_gif()
    print("[gifs] damage evolution ...")
    make_damage_gif()
    print("[gifs] done.")
