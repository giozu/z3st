#!/usr/bin/env python3
"""
Animation of the VP1 thermal bow: the gradient is applied, held and removed.

    python animate_vp1.py                 # output/vp1_bow.gif
    python animate_vp1.py --case flats    # flats/output/vp1_bow.gif
    python animate_vp1.py --warp 30 --out bow.gif
    python animate_vp1.py --show          # interactive window instead of a file

The problem is linear thermoelastic, so the state at load fraction s is
T = T_iso + s (T_vp1 - T_iso),  u = u_iso + s (u_vp1 - u_iso),  with T_iso,
u_iso the solution at uniform 400 °C. With T_ref = 400 °C (steel_iaea.yaml)
u_iso = 0; otherwise (flats_1515Ti, T_ref = 300 K) it is read from
<case>/iso/output/<fields>_iso.h5, the run at uniform 400 °C made by Allrun. Only the bow
s (u_vp1 - u_iso) is drawn, magnified by --warp.
"""

import argparse
import os
import sys

import numpy as np
import pyvista as pv

from view_vp1 import FIELDS, build_grid, h5_path

HERE = os.path.dirname(os.path.abspath(__file__))


def load_fraction(n_rest=12, n_ramp=36, n_hold=24):
    """s(frame): rest, ramp up, hold, ramp down, rest."""
    up = 0.5 - 0.5 * np.cos(np.linspace(0, np.pi, n_ramp))
    return np.concatenate([np.zeros(n_rest), up, np.ones(n_hold), up[::-1], np.zeros(n_rest)])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="", choices=list(FIELDS))
    ap.add_argument("--warp", type=float, default=20.0)
    ap.add_argument("--out")
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    a.out = a.out or os.path.join(HERE, a.case, "output", "vp1_bow.gif")

    full = build_grid("1", h5_path(a.case))
    if os.path.exists(h5_path(a.case, "_iso")):
        iso = build_grid("1", h5_path(a.case, "_iso"))
        T_iso, u_iso = iso["T (°C)"], iso["Displacement"]
    else:                                   # T_ref = 400 °C: no strain at 400 °C
        u = np.abs(full["Displacement"])
        if u[full.points[:, 2] < 1.0].max() > 1e-3 * u.max():
            sys.exit(f"[ERROR] {h5_path(a.case, '_iso')} missing and the duct "
                     "below the core moves: T_ref is not 400 °C")
        T_iso, u_iso = np.full(full.n_points, 400.0), 0.0
    dT = full["T (°C)"] - T_iso
    bow = full["Displacement"] - u_iso

    z = full.points[:, 2]
    tip = z > z.max() - 1e-6
    tip_bow = 1e3 * bow[tip, 1].mean()       # mm, lateral (y)
    s_all = load_fraction()
    t = np.arange(len(s_all))

    grid = full.copy()
    grid.point_data.clear()
    grid.cell_data.clear()
    grid["T (°C)"] = T_iso.copy()
    ref = grid.copy()                       # undeformed outline

    p = pv.Plotter(shape=(1, 2), window_size=(1200, 800), off_screen=not a.show)
    p.subplot(0, 0)
    p.add_mesh(ref.extract_feature_edges(), color="gray", opacity=0.5)
    p.add_mesh(grid, scalars="T (°C)", cmap="coolwarm",
               clim=[T_iso.min(), (T_iso + dT).max()],
               scalar_bar_args={"title": "T (°C)", "vertical": True})
    label = p.add_text("", font_size=11)
    p.add_text(f"bow x{a.warp:g}", position="lower_left", font_size=9)
    grid.points = ref.points + a.warp * bow   # frame the fully bowed shape
    p.view_yz()
    p.reset_camera()
    grid.points = ref.points.copy()

    p.subplot(0, 1)
    chart = pv.Chart2D(x_label="frame", y_label="tip deflection (mm)")
    chart.line(t, s_all * tip_bow, color="gray", style="--")
    trace = chart.line([0], [0], color="tab:red", width=3)
    marker = chart.scatter([0], [0], color="tab:red", size=12)
    chart.y_range = [min(tip_bow, 0) * 1.15, max(tip_bow, 0) * 1.15 + 0.5]
    p.add_chart(chart)

    if not a.show:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        p.open_gif(a.out, fps=24)
    else:
        p.show(interactive_update=True)

    for i, s in enumerate(s_all):
        grid.points = ref.points + a.warp * s * bow
        grid["T (°C)"] = T_iso + s * dT
        label.SetText(2, f"Delta T fraction = {s:4.2f}\n"
                         f"tip deflection = {s * tip_bow:6.2f} mm")
        trace.update(t[: i + 1], s_all[: i + 1] * tip_bow)
        marker.update([t[i]], [s * tip_bow])
        if a.show:
            p.update()
        else:
            p.write_frame()

    print(f"[INFO] peak tip bow {tip_bow:.2f} mm, residual after unloading "
          f"{s_all[-1] * tip_bow:.2e} mm")
    if not a.show:
        print(f"[INFO] {a.out}")
    p.close()
