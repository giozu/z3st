#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Analytic slab solution vs z3st along the three lines of the local frame of a
flat, each with the other two coordinates fixed:

  1. xi (through the wall), at mid-flat and mid-height;
  2. s (along the flat), at mid-height;
  3. z (along the axis), at mid-flat.

Fields are read from the last step of ``output/fields.xdmf`` and mapped into
the local frame by ``case_params.local_fields``; the six flats collapse onto
one curve. Figures are written to ``output/profile_{xi,s,z}.png``.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from case_params import (
    DT, H, L_MID, NT, T_I, T_WALL, XI_CELLS,
    along_axis, along_flat, check_consistency, local_fields, on_layer, profile,
    sigma_wall, temperature, through_wall,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")

# cell layer -> (colour, label), shared by the s and z plots
_LAYER_STYLE = {0: ("tab:blue", "inner"), NT - 1: ("tab:red", "outer")}


def _overlay_layers(ax, f, coord, mask, scale):
    """sigma_ss and sigma_zz of the inner and outer cell layers along ``coord``."""
    for j, (color, label) in _LAYER_STYLE.items():
        layer = mask & on_layer(f, XI_CELLS[j])
        c, ss = profile(f, coord, "ss", layer)
        _, zz = profile(f, coord, "zz", layer)
        ax.plot(c * scale, ss / 1e6, "s", color=color, ms=4, mfc="none",
                label=rf"z3st $\sigma_{{ss}}$, {label} layer")
        ax.plot(c * scale, zz / 1e6, "^", color=color, ms=4, mfc="none",
                label=rf"z3st $\sigma_{{zz}}$, {label} layer")
        ax.axhline(sigma_wall(XI_CELLS[j]) / 1e6, color=color, lw=1.5,
                   label=f"analytic, {label} layer")


# ----------------------------------------------------------------------
# 1. THROUGH THE WALL (xi)
# ----------------------------------------------------------------------
def plot_xi(nodes, cells):
    xi = np.linspace(-T_WALL / 2, T_WALL / 2, 201)
    fig, (axT, axS) = plt.subplots(1, 2, figsize=(12, 4.5))

    c, T = profile(nodes, "xi", "T", through_wall(nodes))
    axT.plot(xi * 1e3, temperature(xi), "k-", lw=2, label="analytic (Kirchhoff)")
    axT.plot(c * 1e3, T, "o", color="tab:red", mfc="none", label="z3st, nodes")
    axT.set_xlabel(r"$\xi$ (mm)   inner $\rightarrow$ outer")
    axT.set_ylabel("T (K)")
    axT.set_title(r"Temperature, mid-flat, mid-height")

    m = through_wall(cells)
    axS.plot(xi * 1e3, sigma_wall(xi) / 1e6, "k-", lw=2,
             label=r"analytic $\sigma_{ss} = \sigma_{zz}$")
    for comp, mk, color in (("ss", "s", "tab:green"), ("zz", "^", "tab:blue"),
                            ("nn", "x", "tab:gray")):
        c, v = profile(cells, "xi", comp, m)
        axS.plot(c * 1e3, v / 1e6, mk, color=color, ms=7, mfc="none",
                 label=rf"z3st $\sigma_{{{comp}}}$")
    axS.axhline(0, color="gray", lw=0.8)
    axS.set_xlabel(r"$\xi$ (mm)   inner $\rightarrow$ outer")
    axS.set_ylabel("stress (MPa)")
    axS.set_title("Stresses, mid-flat, mid-height")

    for ax in (axT, axS):
        ax.grid(ls=":", alpha=0.6)
        ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------
# 2. ALONG THE FLAT (s)
# ----------------------------------------------------------------------
def plot_s(nodes, cells):
    fig, (axT, axS) = plt.subplots(1, 2, figsize=(12, 4.5))

    c, T = profile(nodes, "s", "T", along_flat(nodes) & on_layer(nodes, -T_WALL / 2))
    axT.plot(c * 1e3, T, "o", color="tab:red", ms=4, mfc="none", label=r"z3st, inner face")
    axT.axhline(T_I, color="k", lw=2, label=r"analytic $T_i$")
    axT.set_ylabel("T (K)")
    axT.set_title("Inner-face temperature, mid-height")

    _overlay_layers(axS, cells, "s", along_flat(cells), 1e3)
    axS.set_ylabel("stress (MPa)")
    axS.set_title("Stresses in the cell layers, mid-height")

    for ax in (axT, axS):
        ax.axvspan(-L_MID / 4 * 1e3, L_MID / 4 * 1e3, color="gray", alpha=0.1,
                   label="checked window")
        ax.set_xlabel(r"$s$ (mm)   corner $\leftarrow$ mid-flat $\rightarrow$ corner")
        ax.grid(ls=":", alpha=0.6)
        ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------
# 3. ALONG THE AXIS (z)
# ----------------------------------------------------------------------
def plot_z(nodes, cells):
    fig, (axT, axS) = plt.subplots(1, 2, figsize=(12, 4.5))

    c, T = profile(nodes, "z", "T", along_axis(nodes) & on_layer(nodes, -T_WALL / 2))
    axT.plot(c, T, "o", color="tab:red", ms=4, mfc="none", label="z3st, inner face")
    axT.axhline(T_I, color="k", lw=2, label=r"analytic $T_i$")
    axT.set_ylabel("T (K)")
    axT.set_title("Inner-face temperature, mid-flat")
    axT.set_ylim(T_I - 0.05 * DT, T_I + 0.05 * DT)  # uniform to round-off
    axT.ticklabel_format(axis="y", useOffset=False)

    _overlay_layers(axS, cells, "z", along_axis(cells), 1.0)
    axS.set_ylabel("stress (MPa)")
    axS.set_title("Stresses in the cell layers, mid-flat")

    for ax in (axT, axS):
        ax.axvspan(H / 4, 3 * H / 4, color="gray", alpha=0.1, label="checked window")
        ax.set_xlabel(r"$z$ (m)")
        ax.grid(ls=":", alpha=0.6)
        ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    nodes, cells = local_fields()
    problems = check_consistency(cells)
    if problems:
        raise RuntimeError("; ".join(problems))
    for name, plot in (("xi", plot_xi), ("s", plot_s), ("z", plot_z)):
        path = os.path.join(OUT, f"profile_{name}.png")
        plot(nodes, cells).savefig(path, dpi=150)
        print(f"[INFO] {path}")
