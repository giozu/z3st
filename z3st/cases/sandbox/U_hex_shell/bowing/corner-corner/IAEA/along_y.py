#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 / VP2: T, stress and strain at mid-wall against the global y, on
the half perimeter x >= 0, at the heights ``along_s_z`` of post.yaml.

With the corners on +-y ('vertical') T depends on y only and the section is
symmetric about x = 0: on the half x >= 0 each y is one point of the wall,
from the cold corner 1 (y = +R) to the hot corner 4 (y = -R), R the
mid-wall circumradius. The components stay in the local frame of hex_post
(s along the flat): hoop sigma_ss, sigma_zz, eps_ss, eps_zz, with the free
thermal strain alpha (T - T_ref).

Usage (z3st conda env)::

    python along_y.py <case> [<case> ...]     # e.g. VP1 VP2 VP1/z_uniform

Output: ``U_hex_shell/results/IAEA_<case>_along_y.png``.
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 3))
sys.path.insert(0, U_HEX)
import hex_post as hp  # noqa: E402

C2K = 273.15


def profiles(case_dir):
    """Per height of along_s_z: y, T, sigma_ss, sigma_zz, eps_ss, eps_zz and
    alpha (T - T_ref) at mid-wall on x >= 0, sorted by y."""
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    path = lambda p: os.path.normpath(os.path.join(case_dir, p))
    mat = hp.load_material(path(cfg["material"]))
    d = hp.load_z3st(path(cfg["fields"]), temperature=cfg["temperature"], case_dir=case_dir)
    res = hp.analyse(d, mat, cfg.get("support", "clamped"))
    module, func = cfg["temperature"].rsplit(".", 1)
    m = hp._import(module, case_dir)
    if abs(np.cos(m.PHI_HOT)) > 1e-9:
        sys.exit(f"{case_dir}: corners not on +-y (orientation 'vertical' expected)")
    T_of = getattr(m, func)

    fr = res["frame"]
    c = fr["cells"]
    x_mid, y_mid, _ = hp.perimeter_point(fr, c["p"], 0.0)
    out = {}
    for name, z in cfg["along_s_z"].items():
        sel = hp._nearest(c["z"], z) & (x_mid > 0)
        pu, inv = np.unique(np.round(c["p"][sel], 9), return_inverse=True)
        mean = lambda v: np.bincount(inv, weights=v[sel]) / np.bincount(inv)   # through the wall
        x, y, _ = hp.perimeter_point(fr, pu, 0.0)
        T = np.ravel(T_of(np.array([x, y, np.full_like(x, z)])))
        k = np.argsort(y)
        out[name] = {"z": z, "y": y[k], "T": T[k], "s_ss": mean(res["ss"])[k], "s_zz": mean(res["data"]["S"][:, 2, 2])[k],
                     "e_ss": mean(res["eps_ss"])[k], "e_zz": mean(res["eps_zz"])[k],
                     "e_th": mat["alpha"] * (T[k] - mat["T_ref"])}
    R = fr["a_mid"] / np.cos(np.pi / 6)
    corners = {"1": R, "2, 6": R / 2, "3, 5": -R / 2, "4": -R}
    if np.sin(m.PHI_HOT) > 0:   # hot corner 4 at +y
        corners = {k: -v for k, v in corners.items()}
    return out, corners


def figure(case, prof, corners):
    fig, axes = plt.subplots(3, 1, figsize=(12, 13), sharex=True)
    for (name, p), col in zip(prof.items(), hp.LEVEL_COLORS):
        y = p["y"] * 1e3
        axes[0].plot(y, p["T"] - C2K, ".-", color=col, lw=1.8, label=f"{name}, z = {p['z'] * 1e3:.0f} mm")
        axes[1].plot(y, p["s_ss"] * 1e-6, ".-", color=col, lw=1.8)
        axes[1].plot(y, p["s_zz"] * 1e-6, ".--", color=col, lw=1.4)
        axes[2].plot(y, p["e_ss"] * 1e3, ".-", color=col, lw=1.8)
        axes[2].plot(y, p["e_zz"] * 1e3, ".--", color=col, lw=1.4)
        axes[2].plot(y, p["e_th"] * 1e3, ":", color=col, lw=2.0)
    # line styles (colours = heights, legend of the T panel)
    for ls, lab in (("-", "σ_ss (hoop)"), ("--", "σ_zz")):
        axes[1].plot([], [], ls, color="k", label=lab)
    for ls, lab in (("-", "ε_ss (hoop)"), ("--", "ε_zz"), (":", "α (T - T_ref), free thermal")):
        axes[2].plot([], [], ls, color="k", label=lab)
    axes[0].set_ylabel("T (C)")
    axes[1].set_ylabel("stress (MPa)")
    axes[1].axhline(0, color="k", lw=0.8)
    axes[2].set_ylabel("total strain (10⁻³)")
    axes[2].set_xlabel("y at mid-wall (mm), half perimeter x >= 0")
    for ax in axes:
        for yc in corners.values():
            ax.axvline(yc * 1e3, color="k", lw=0.8, alpha=0.5)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="best")
    for lab, yc in corners.items():
        axes[0].text(yc * 1e3, 1.01, f"corner {lab}", ha="center", va="bottom", fontsize=9,
                     transform=axes[0].get_xaxis_transform())
    return hp._finish(fig, f"IAEA {case}: along y at mid-wall",
                      "cell centres (dots), cells through the wall averaged; components in the local frame (s along the flat); "
                      "vertical lines = corners", top=0.93)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    if len(sys.argv) < 2:
        sys.exit("usage: python along_y.py <case> [<case> ...]")
    os.makedirs(hp.RESULTS, exist_ok=True)
    for arg in sys.argv[1:]:
        case_dir = os.path.abspath(arg)
        case = os.path.relpath(case_dir, HERE).replace(os.sep, "_")
        prof, corners = profiles(case_dir)
        fig = figure(case, prof, corners)
        path = os.path.join(hp.RESULTS, f"IAEA_{case}_along_y.png")
        fig.savefig(path, dpi=130)
        plt.close(fig)
        print(f"[INFO] {path}")
        for name, p in prof.items():
            print(f"{case} {name}: max|σ_ss| = {np.abs(p['s_ss']).max() * 1e-6:.1f} MPa, "
                  f"max|σ_zz| = {np.abs(p['s_zz']).max() * 1e-6:.1f} MPa")


if __name__ == "__main__":
    main()
