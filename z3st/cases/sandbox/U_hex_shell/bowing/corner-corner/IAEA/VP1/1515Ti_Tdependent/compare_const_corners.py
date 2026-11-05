#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 (z-gradient): stress along z at the hot (T2) and cold (T1) corner,
z3st with the T-dependent 15-15Ti card, Sobolev (this folder) and Di Marcello
(DiMarcello/) sets, and with the constant card (../), same mesh, BCs and T.
Layout of ../vp1_corners.py.

  * sigma_zz, hoop sigma_ss and von Mises of the two DG0 cells next to the
    mid-wall corner (one per flat, averaged), per element layer;
  * last two panels: error of the constant card, constant - law, of sigma_zz
    and sigma_vm, with max |constant - law| in the legend.

Runs loaded as in compare_const.py (T checked equal on the nodes).

Usage (z3st conda env)::

    python compare_const_corners.py

Output: ``U_hex_shell/results/bowing_corner-corner_IAEA_VP1/1515Ti_Tdep/vs_const/IAEA_VP1_1515Ti_Tdep_vs_const_corners_z.png``.
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import compare_const as cc  # noqa: E402
import vp1_temperature as vt  # noqa: E402

TITLE = "IAEA_VP1 (z-gradient): 1515Ti_Tdep Sobolev and Di Marcello vs constant card, hot and cold corner along z"
CORNERS = {"hot corner (T2)": ("tab:red", vt.PHI_HOT), "cold corner (T1)": ("tab:blue", vt.PHI_HOT + np.pi)}
RUNS = {"Sobolev": "-", "DiMarcello": "-.", "constant": "--"}
STRESSES = (("s_zz", "σ_zz"), ("s_ss", "σ_ss (hoop)"), ("s_vm", "σ_vm"))


def corner_stress(res):
    """Per corner: z of the cell layers and sigma_zz, sigma_ss, sigma_vm (Pa)
    of the two cells next to the mid-wall corner."""
    C = res["frame"]["cells"]["z"]
    d = res["data"]
    Cxy = d["X"][d["cells"]].mean(axis=1)[:, :2]
    r_mid = res["frame"]["a_mid"] / np.cos(np.pi / 6)   # mid-wall circumradius
    out = {}
    for name, (_, phi) in CORNERS.items():
        dc = np.hypot(*(Cxy - r_mid * np.array([np.cos(phi), np.sin(phi)])).T)
        sel = dc < dc.min() * (1 + 1e-6) + 1e-12
        zc, inv = np.unique(np.round(C[sel], 9), return_inverse=True)
        mean = lambda v: np.bincount(inv, weights=v[sel]) / np.bincount(inv)
        c = r_mid * np.array([np.cos(phi), np.sin(phi)])
        T = np.ravel(vt.T_vp1(np.array([np.full_like(zc, c[0]), np.full_like(zc, c[1]), zc])))
        out[name] = {"z": zc, "s_zz": mean(d["S"][:, 2, 2]), "s_ss": mean(res["ss"]), "s_vm": mean(res["vm"]),
                     "T": T}
    hot, cold = (out[n]["T"].max() for n in CORNERS)
    if hot <= cold:
        sys.exit(f"hot corner at {hot:.1f} K, cold at {cold:.1f} K: PHI_HOT does not point to the hot corner")
    return out


def figure(p):
    fig, axes = plt.subplots(1, 5, figsize=(25, 8.5), sharey=True)
    for ax, (key, lab) in zip(axes, STRESSES):
        for run, ls in RUNS.items():
            for name, (col, _) in CORNERS.items():
                c = p[run][name]
                ax.plot(c[key] * 1e-6, c["z"] * 1e3, ls, color=col, lw=1.5, label=f"{run}, {name}")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel(f"{lab} (MPa)")
        ax.set_title(f"{lab} (corner cells, DG0)")

    for ax, (key, lab) in zip(axes[3:], (STRESSES[0], STRESSES[2])):
        for run in list(RUNS)[:-1]:
            for name, (col, _) in CORNERS.items():
                t, c = p[run][name], p["constant"][name]
                if not np.allclose(t["z"], c["z"]):
                    sys.exit(f"{name}: cell layers of {run} and constant differ")
                d = c[key] - t[key]   # Pa
                ax.plot(d * 1e-6, t["z"] * 1e3, RUNS[run], color=col, lw=1.5,
                        label=f"constant − {run}, {name}: max |Δ| = {np.abs(d).max() * 1e-6:.3f} MPa "
                              f"({100 * np.abs(d).max() / np.abs(t[key]).max():.1f} %)")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel(f"Δ{lab}, constant − T-dependent (MPa)")
        ax.set_title(f"Error of the constant card: {lab}")

    axes[0].set_ylabel("height above support z (mm)")
    for ax in axes:
        for zc in (1500, 2500):
            ax.axhline(zc, color="0.6", lw=0.8, ls=":")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8 if ax in axes[:3] else 6.5, loc="upper left")
    axes[0].text(0.02, 1540, "core bottom", fontsize=8, color="0.4", transform=axes[0].get_yaxis_transform())
    axes[0].text(0.02, 2540, "core top", fontsize=8, color="0.4", transform=axes[0].get_yaxis_transform())
    hot, cold = (p["constant"][n]["T"].max() - 273.15 for n in CORNERS)
    return cc.hp._finish(fig, TITLE, f"corners at mid-wall; T above the core: hot {hot:.1f} C, cold {cold:.1f} C; "
                         "solid: Sobolev, dash-dot: Di Marcello, dashed: constant card; same mesh, BCs and T; "
                         "bottom clamped\n(%): max |Δ| over max |T-dependent| of that corner along z", top=0.91)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    loaded = cc.load_runs()
    if loaded is None:
        return
    runs = loaded[0]
    p = {run: corner_stress(res) for run, res in runs}

    os.makedirs(cc.OUT, exist_ok=True)
    path = os.path.join(cc.OUT, f"{cc.PREFIX}corners_z.png")
    fig = figure(p)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")
    for run in RUNS:
        for name, c in p[run].items():
            print(f"{run}, {name}: max|σ_zz| = {np.abs(c['s_zz']).max() * 1e-6:.2f} MPa, "
                  f"max σ_vm = {c['s_vm'].max() * 1e-6:.2f} MPa at z = {c['z'][np.argmax(c['s_vm'])] * 1e3:.0f} mm")


if __name__ == "__main__":
    main()
