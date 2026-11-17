#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Profiles and 3D views of a step, read only from ``output/fields.xdmf``.

Geometry and fields come from ``gui_fields.load_grid``, so the plots show the
stored run whatever the case YAML files currently say. Each flat is drawn on
its own, because the thermal BCs may differ from flat to flat.

Local frame of a flat: xi through the wall (0 at mid-wall, outward), s along
the flat (0 at mid-flat), z along the axis; theta is the angle around the axis
from the first corner.

  1. xi (through the wall), at mid-flat and mid-height, one row per flat
     with T, stresses (nn, ss, zz, vm) and strains (nn, ss, zz);
  2. theta (around the shell), at mid-height, inner and outer layers;
  3. z (along the axis), at mid-flat and mid-wall, one curve per flat.

The theta and z figures have eight panels: T, sigma_vm, sigma_nn, sigma_ss,
sigma_zz, eps_nn, eps_ss, eps_zz.
Figures are written to ``output/profile_{xi,theta,z}.png``; 3D views of the step to
``output/3d_{mesh,original,deformed,stress,vonmises,strain}.png``.

Usage: ``python plots.py`` (last step).
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

HERE = os.path.dirname(os.path.abspath(__file__))
# [TBC] rule 4 — drop the chdir comment? loses: gui_fields reads cwd at import,
# so chdir must precede the import.
os.chdir(HERE)  # gui_fields reads the case in the working directory
sys.path.insert(0, os.path.join(HERE, "..", "..", ".."))  # U_hex_shell/gui_fields.py
from gui_fields import A0, H, flat_frame, load_grid, plot_3d  # noqa: E402

OUT = os.path.join(HERE, "output")

# physical group of each flat, by sector index (see mesh.geo)
FLATS = ("nord-est", "nord", "nord-west", "sud-west", "sud", "sud-est")
COLORS = plt.cm.tab10(np.arange(len(FLATS)))
PANELS = (("T (K)", "T (K)"),
          ("sigma_vm (MPa)", r"$\sigma_{vm}$ (MPa)"),
          ("sigma_nn (MPa)", r"$\sigma_{nn}$ (MPa)"),
          ("sigma_ss (MPa)", r"$\sigma_{ss}$ (MPa)"),
          ("sigma_zz (MPa)", r"$\sigma_{zz}$ (MPa)"),
          ("eps_nn (-)", r"$\varepsilon_{nn}$ (-)"),
          ("eps_ss (-)", r"$\varepsilon_{ss}$ (-)"),
          ("eps_zz (-)", r"$\varepsilon_{zz}$ (-)"))


# --. local frame --..
def local_coords(p):
    """Sector, apothem a, s, s_frac, z and theta (deg) of points ``p`` (N x 3);
    s_frac = s / half-flat length at the point's apothem, so mid-flat rows line up
    across the wall layers of the trapezoidal sectors."""
    n, tang = flat_frame(p[:, 0], p[:, 1])
    phi_n = np.arctan2(n[:, 1], n[:, 0])
    sector = np.round((phi_n - A0 - np.pi / 6) / (np.pi / 3)).astype(int) % 6
    theta = np.degrees(np.arctan2(p[:, 1], p[:, 0]) - A0) % 360
    a, s = (p * n).sum(1), (p * tang).sum(1)
    return {"sector": sector, "a": a, "s": s, "s_frac": s / (a * np.tan(np.pi / 6)),
            "z": p[:, 2], "theta": theta}


def local_data(step_index=-1):
    """Nodal T and cell stresses of a step, with their local coordinates."""
    grid = load_grid(step_index)
    nodes = local_coords(grid.points)
    cells = local_coords(grid.cell_centers().points)
    a_mid = 0.5 * (nodes["a"].min() + nodes["a"].max())   # mid-wall apothem
    for f in (nodes, cells):
        f["xi"] = f.pop("a") - a_mid
    nodes.update(grid.point_data)
    cells.update(grid.cell_data)
    return grid, nodes, cells


def nearest(v, target):
    """The row(s) of ``v`` closest to ``target`` (both rows when equidistant)."""
    d = np.abs(v - target)
    return d < d.min() + 1e-6 * (np.ptp(v) or 1.0)


def profile(f, coord, field, mask):
    """Mean of ``field`` over the points of ``mask`` sharing a ``coord`` value."""
    c = np.round(f[coord][mask], 9)
    axis, inverse = np.unique(c, return_inverse=True)
    values = np.bincount(inverse, weights=f[field][mask]) / np.bincount(inverse)
    return axis, values


def _figure(title):
    """Two-column grid for the PANELS, row by row; a spare slot is hidden."""
    rows = -(-len(PANELS) // 2)
    fig, axes = plt.subplots(rows, 2, figsize=(12, 3.6 * rows), sharex=True)
    fig.suptitle(title)
    axes = axes.ravel()
    for ax in axes[len(PANELS):]:
        ax.set_visible(False)
    return fig, axes[:len(PANELS)]


def _finish(fig, axes, xlabel):
    for ax, (_, label) in zip(axes, PANELS):
        ax.set_ylabel(label)
        ax.grid(ls=":", alpha=0.6)
        ax.ticklabel_format(axis="y", useOffset=False)
    for ax in axes[-2:]:   # bottom of each column
        ax.set_xlabel(xlabel)
        ax.xaxis.set_tick_params(labelbottom=True)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    return fig


def _source(nodes, cells, key):
    return nodes if key in nodes else cells


# --. 1. through the wall (xi), one row per flat --..
XI_GROUPS = (
    ("Temperature", "T (K)", (("T (K)", "T", "o", "tab:red"),)),
    ("Stress", "stress (MPa)", (("sigma_nn (MPa)", r"$\sigma_{nn}$", "x", "tab:gray"),
                                ("sigma_ss (MPa)", r"$\sigma_{ss}$", "s", "tab:green"),
                                ("sigma_zz (MPa)", r"$\sigma_{zz}$", "^", "tab:blue"),
                                ("sigma_vm (MPa)", r"$\sigma_{vm}$", "o", "k"))),
    ("Strain", "strain (-)", (("eps_nn (-)", r"$\varepsilon_{nn}$", "x", "tab:gray"),
                              ("eps_ss (-)", r"$\varepsilon_{ss}$", "s", "tab:green"),
                              ("eps_zz (-)", r"$\varepsilon_{zz}$", "^", "tab:blue"))),
)


# row order of plot_xi, top to bottom: nord on top, sud at the bottom
XI_ROWS = ("nord", "nord-west", "nord-est", "sud-west", "sud-est", "sud")


def plot_xi(nodes, cells):
    """T, stresses and strains through the wall, mid-flat, mid-height: one row
    per flat (``XI_ROWS``), one column per group. T and stress share their y
    limits down the column; strain shares the y span only, each row centred
    on its own flat."""
    fig, axes = plt.subplots(len(XI_ROWS), 3, figsize=(15, 3 * len(XI_ROWS)), sharex=True)
    extent = np.empty((len(XI_ROWS), 3, 2))   # data min/max of each panel
    fig.suptitle("Through the wall, mid-flat, mid-height")
    for i, (row, name) in enumerate(zip(axes, XI_ROWS)):
        m = FLATS.index(name)
        for j, (ax, (title, ylabel, curves)) in enumerate(zip(row, XI_GROUPS)):
            values = []
            for key, label, marker, color in curves:
                f = _source(nodes, cells, key)
                line = nearest(f["s_frac"], 0.0) & nearest(f["z"], H / 2) & (f["sector"] == m)
                c, v = profile(f, "xi", key, line)
                ax.plot(c * 1e3, v, marker + "-", color=color, ms=5, mfc="none", label=label)
                values.append(v)
            values = np.concatenate(values + ([[0.0]] if title == "Stress" else []))
            extent[i, j] = values.min(), values.max()
            if title == "Stress":
                ax.axhline(0, color="gray", lw=0.6)
            ax.set_ylabel(f"{name}\n{ylabel}" if ax is row[0] else ylabel)
            ax.grid(ls=":", alpha=0.6)
            ax.ticklabel_format(axis="y", useOffset=False)
    for ax, (title, _, _) in zip(axes[0], XI_GROUPS):
        ax.set_title(title)
        ax.legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\xi$ (mm)   inner $\rightarrow$ outer")
    # T needs at least 1 K span: a uniform T is not zoomed into round-off
    for j, (title, _, _) in enumerate(XI_GROUPS):
        lo, hi = extent[:, j, 0], extent[:, j, 1]
        if title == "Strain":   # common span, each row centred on its own data
            mid, half = 0.5 * (lo + hi), np.full(len(lo), 0.5 * (hi - lo).max())
        else:                   # common limits
            mid, half = np.full(len(lo), 0.5 * (lo.min() + hi.max())), np.full(len(lo), 0.5 * (hi.max() - lo.min()))
        if title == "Temperature":
            half = np.maximum(half, 0.5)
        for ax, c, h in zip(axes[:, j], mid, half * 1.1):
            ax.set_ylim(c - h, c + h)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    return fig


# --. 2. around the shell (theta) --..
def plot_theta(nodes, cells):
    fig, axes = _figure("Around the shell, mid-height")
    for ax, (key, _) in zip(axes, PANELS):
        f = _source(nodes, cells, key)
        ring = nearest(f["z"], H / 2)
        for side, color in (("inner", "tab:blue"), ("outer", "tab:red")):
            layer = nearest(f["xi"], f["xi"].min() if side == "inner" else f["xi"].max())
            c, v = profile(f, "theta", key, ring & layer)
            where = "face" if f is nodes else "layer"
            ax.plot(c, v, "-", color=color, lw=1.2, label=f"{side} {where}")
        for m, name in enumerate(FLATS):
            ax.axvline(60 * m, color="gray", lw=0.6)
            if ax is axes[0]:
                ax.text(60 * m + 30, 1.02, name, transform=ax.get_xaxis_transform(),
                        ha="center", fontsize=8)
    axes[0].set_xticks(np.arange(0, 361, 60))
    return _finish(fig, axes, r"$\theta$ (deg) from the first corner")


# --. 3. along the axis (z) --..
def plot_z(nodes, cells):
    fig, axes = _figure("Along the axis, mid-flat, mid-wall")
    for ax, (key, _) in zip(axes, PANELS):
        f = _source(nodes, cells, key)
        line = nearest(f["s_frac"], 0.0) & nearest(f["xi"], 0.0)
        for m, name in enumerate(FLATS):
            c, v = profile(f, "z", key, line & (f["sector"] == m))
            ax.plot(c, v, "-", color=COLORS[m], lw=1.2, label=name)
    return _finish(fig, axes, r"$z$ (m)")


if __name__ == "__main__":
    grid, nodes, cells = local_data()
    for name, plot in (("xi", plot_xi), ("theta", plot_theta), ("z", plot_z)):
        path = os.path.join(OUT, f"profile_{name}.png")
        fig = plot(nodes, cells)
        fig.savefig(path, dpi=150)
        plt.close(fig)
        print(f"[INFO] {path}")
    plot_3d(grid, OUT)
