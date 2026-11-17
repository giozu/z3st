#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Legend of the stress and strain components in the figures of ``results/``.

Usage::

    python results_components.py

Writes ``results/components.png``, a schematic that reads no case output:

  - section of the duct, for orientation flat (V1 at 0°: inner_q_outer_T,
    T_asymmetries) and vertical (V1 at 30°: IAEA VP1, VP2): vertices V1..V6,
    flats, the local frame (n, s, z) on one flat, the perimeter coordinate p
    and the hoop direction on every flat;
  - wall element: sigma_ss, sigma_zz, sigma_nn and the shear sigma_sz on its faces;
  - table: definition, meaning, plotted unit and figures of every component.

Frame as in hex_post.frame: n outward normal of the flat, t = z x n along
the flat (s, 0 at mid-flat), z the duct axis; (n, s, z) is right-handed.
p runs at mid-wall from V1 counter-clockwise: the x-axis of the along-s
figures.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hex_post as hp  # noqa: E402

COL = {"n": "tab:purple", "s": "tab:green", "z": "tab:blue", "shear": "0.45"}

# (component, definition, meaning, plotted as, figures <case>_*.png)
ROWS = (
    ("σ_ss", "t·σ·t", "hoop stress, along the flat (s)", "MPa",
     "stress_z, along_s_stress, along_y, corners, max_stress"),
    ("σ_zz", "σ_zz (global z)", "axial stress, along the duct axis", "MPa",
     "along_y, corners, max_stress"),
    ("σ_nn", "n·σ·n", "through the wall; ≈ 0 off the clamped base (free faces)", "—",
     "not plotted"),
    ("σ_vm", "√(3/2 s:s), s = σ − (tr σ / 3) I", "von Mises, all six components", "MPa",
     "stress_z, along_s_stress, corners, max_stress"),
    ("ε_ss", "t·ε·t", "hoop strain, stretch along s", "-  (×10³ in along_y)",
     "along_s_strain_sz, along_y, corners"),
    ("ε_zz", "∂u_z/∂z", "axial strain, stretch along the axis", "-  (×10³ in along_y)",
     "along_s_strain_sz, along_y, corners"),
    ("ε_nn", "n·ε·n = ∂u_n/∂ξ", "through-wall strain, change of thickness", "-",
     "along_s_normal"),
    ("u_n", "(u − ū_section)·n", "outward wall displacement, ū = mean nodal u of the section", "mm",
     "along_s_normal"),
    ("ε_th", "α (T − T_ref)", "free thermal strain (reference curve)", "-  (×10³ in along_y)",
     "along_y, corners"),
)

NOTES = (
    "ε is the total strain sym∇u (z3st field Strain); σ = ℂ : (ε − ε*), ε* = α (T − T_ref) I "
    "(+ the card eigenstrain / creep strain). Where the wall expands freely ε = ε_th and σ = 0.",
    "Coordinates: along-s figures plot against p (perimeter at mid-wall from V1, ticks at V1..V6); max_stress gives s "
    "from mid-flat and ξ from mid-wall. Vertical meshes (IAEA): flats named by the angle of n, not by compass.",
    "curvature (default 'rotation'): κ = −d/dz of the slope of the plane fitted to u_z on each node plane, not "
    "from the output ε_zz.",
    "Sign: tension > 0, compression < 0. Cell values (DG0) for σ and ε, nodal values for u_n.",
    "From the global x, y components (ParaView Stress_XX, ...), flat with outward normal at angle φ: "
    "σ_ss = σ_xx sin²φ − 2 σ_xy sinφ cosφ + σ_yy cos²φ,  σ_nn = σ_xx cos²φ + 2 σ_xy sinφ cosφ + σ_yy sin²φ.",
    "along_s_zmid / along_s_zmaxvm (no along_s_z in post.yaml): one panel each for σ_vm, σ_ss, ε_ss, ε_zz, "
    "ε_nn and u_n; along_s_stress / _strain_sz / _normal split the same panels in pairs.",
)


# --. section --..
def _arrow(ax, x0, y0, dx, dy, color, lw=2.0, style="-|>"):
    ax.annotate("", (x0 + dx, y0 + dy), (x0, y0), arrowprops=dict(arrowstyle=style, color=color, lw=lw))


def fig_section(ax, a0, title):
    """Section with V1 at angle a0 (rad), named as hex_post.frame names the flats."""
    R, t = 1.0, 0.12                                        # vertex radius of mid-wall, wall thickness (schematic)
    a_mid = R * np.cos(np.pi / 6)
    ang = a0 + np.arange(7) * np.pi / 3
    for xi in (-t / 2, t / 2):
        r = (a_mid + xi) / np.cos(np.pi / 6)
        ax.plot(r * np.cos(ang), r * np.sin(ang), color="k", lw=1.0)
    for m in range(6):
        ax.text(1.2 * R * np.cos(ang[m]), 1.2 * R * np.sin(ang[m]), f"V{m + 1}",
                ha="center", va="center", fontsize=10, fontweight="bold")
        phi = ang[m] + np.pi / 6
        deg = int(round(np.degrees(phi))) % 360
        ax.text(0.66 * a_mid * np.cos(phi), 0.66 * a_mid * np.sin(phi), hp.COMPASS.get(deg, f"normal {deg}°"),
                ha="center", va="center", fontsize=7.5, color="0.3", rotation=np.degrees(phi) % 180 - 90)
        # hoop direction on every flat
        n = np.array([np.cos(phi), np.sin(phi)])
        tv = np.array([-n[1], n[0]])
        c = a_mid * n - 0.12 * tv
        if m:                                               # the first flat carries the local frame
            _arrow(ax, *c, *(0.24 * tv), COL["s"], lw=1.2)
    # local frame at mid-flat (s = 0) of the first flat, V1-V2
    phi = a0 + np.pi / 6
    n = np.array([np.cos(phi), np.sin(phi)])
    tv = np.array([-n[1], n[0]])
    c = a_mid * n
    _arrow(ax, *c, *(0.4 * n), COL["n"])
    ax.text(*(c + 0.56 * n), "n, ξ", color=COL["n"], fontsize=10, ha="center", va="center")
    _arrow(ax, *c, *(0.3 * tv), COL["s"])
    ax.text(*(c + 0.3 * tv - 0.16 * n), "s", color=COL["s"], fontsize=10, ha="center", va="center")
    ax.plot(*c, "o", ms=10, mfc="white", mec=COL["z"], mew=1.6)
    ax.plot(*c, ".", ms=4, color=COL["z"])
    ax.text(*(c - 0.13 * tv - 0.1 * n), "z", color=COL["z"], fontsize=10, ha="center", va="center")
    # perimeter coordinate p from V1, outside the wall
    v1 = R * np.array([np.cos(a0), np.sin(a0)])
    _arrow(ax, *(v1 + 0.1 * n - 0.1 * tv), *(0.35 * tv), "0.2", lw=1.2, style="->")
    ax.text(*(v1 + 0.24 * n + 0.08 * tv), "p", fontsize=10, color="0.2", ha="center", va="center")
    _arrow(ax, 0, 0, *(0.38 * n), "0.5", lw=1.0, style="->")
    ax.plot([0, 0.3], [0, 0], color="0.5", lw=0.8, ls=":")
    ax.text(0.22, 0.03 + 0.1 * np.sin(phi), "φ", fontsize=9, color="0.4")
    ax.text(0.0, -0.12, "z ⊙ out of the page", ha="center", va="top", fontsize=8, color="0.3")
    ax.set_aspect("equal")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.4, 1.45)
    ax.axis("off")
    ax.set_title(title, fontsize=10)


# --. wall element --..
def fig_element(ax):
    # plot axes: x = n (through the wall), y = s (along the flat), z = z
    d = np.array([0.5, 1.0, 1.0])                           # thin in n
    for i in range(3):
        for j in (0, 1):
            for k in (0, 1):
                p0, p1 = np.zeros(3), np.zeros(3)
                o = [m for m in range(3) if m != i]
                p0[o[0]] = p1[o[0]] = j * d[o[0]]
                p0[o[1]] = p1[o[1]] = k * d[o[1]]
                p1[i] = d[i]
                ax.plot(*zip(p0, p1), color="k", lw=0.8)

    def arr(p, v, color, lab, ls="-", k=1.25):
        ax.quiver(*p, *v, color=color, lw=2.0, arrow_length_ratio=0.25, linestyle=ls)
        ax.text(*(np.array(p) + k * np.array(v)), lab, color=color, fontsize=11, ha="center")

    arr((d[0] / 2, d[1], d[2] / 2), (0, 0.45, 0), COL["s"], "σ_ss")
    arr((d[0] / 2, 0, d[2] / 2), (0, -0.45, 0), COL["s"], "")
    arr((d[0] / 2, d[1] / 2, d[2]), (0, 0, 0.45), COL["z"], "σ_zz")
    arr((d[0] / 2, d[1] / 2, 0), (0, 0, -0.45), COL["z"], "")
    arr((d[0], d[1] / 2, d[2] / 2), (0.45, 0, 0), COL["n"], "σ_nn ≈ 0", k=2.0)
    ax.quiver(d[0] / 2, d[1] * 0.2, d[2], 0, 0.4, 0, color=COL["shear"], lw=2.0, arrow_length_ratio=0.25,
              linestyle="--")
    ax.text(d[0] / 2, -0.2, d[2] + 0.12, "σ_sz (shear,\nnot plotted)", color=COL["shear"], fontsize=9)
    for v, lab, col in (((0.35, 0, 0), "n (outer face)", COL["n"]), ((0, 0.35, 0), "s", COL["s"]),
                        ((0, 0, 0.35), "z", COL["z"])):
        ax.quiver(1.0, -0.6, -0.3, *v, color=col, lw=1.2, arrow_length_ratio=0.3)
        ax.text(*(np.array((1.0, -0.6, -0.3)) + 1.3 * np.array(v)), lab, color=col, fontsize=8)
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylim(-0.8, 1.6)
    ax.set_zlim(-0.3, 1.5)
    ax.set_box_aspect((1.6, 2.2, 1.8))
    ax.view_init(elev=20, azim=-35)
    ax.set_axis_off()
    ax.set_title("Wall element: tension > 0 shown\n"
                 "σ_ij acts on the face normal to i, along j", fontsize=10)


# --. table --..
def fig_table(ax):
    ax.axis("off")
    cols = ("component", "definition", "meaning", "plotted as", "figures  <case>_*.png")
    tab = ax.table(cellText=[list(r) for r in ROWS], colLabels=cols, bbox=(0, 0.3, 1, 0.7), cellLoc="left",
                   colLoc="left", colWidths=(0.08, 0.17, 0.31, 0.12, 0.32))
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor("0.75")
        if r == 0:
            cell.set_text_props(fontweight="bold")
            cell.set_facecolor("0.92")
        elif c == 0:
            sym = ROWS[r - 1][0]
            color = COL["s"] if "ss" in sym else COL["z"] if "zz" in sym else COL["n"] if sym[-2:] in ("nn", "_n") \
                else "k"
            cell.set_text_props(fontweight="bold", color=color)
    ax.text(0.0, 0.26, "\n".join(f"•  {n}" for n in NOTES), transform=ax.transAxes, va="top", fontsize=8.8,
            wrap=True)


def main():
    fig = plt.figure(figsize=(18, 12.5))
    gs = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.8), hspace=0.08, wspace=0.05)
    fig_section(fig.add_subplot(gs[0, 0]), 0.0,
                "orientation flat (V1 at 0°)\ninner_q_flux, bowing_flat-flat_T_asymmetries")
    fig_section(fig.add_subplot(gs[0, 1]), np.pi / 6,
                "orientation vertical (V1 at 30°)\nbowing_corner-corner_IAEA_VP1, VP2")
    fig_element(fig.add_subplot(gs[0, 2], projection="3d"))
    fig_table(fig.add_subplot(gs[1, :]))
    fig.suptitle("U_hex_shell results: stress σ_ij and strain ε_ij components in the local frame (n, s, z)",
                 fontweight="bold", fontsize=13)
    path = os.path.join(hp.RESULTS, "components.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] {os.path.relpath(path, HERE)}")


if __name__ == "__main__":
    main()
