#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 / VP2: displacement u_y, stress and strain along z at the two corners
of extreme temperature on the bowing axis, the hot corner (T2, direction
PHI_HOT of the case temperature module) and the cold corner opposite (T1).

  * u_y: mid-wall corner (mean of the inner and outer corner node), with the
    z3st section mean and the closed-form bow w(z);
  * stress: sigma_zz, hoop sigma_ss and von Mises of the two DG0 cells next to
    the corner (one per flat, averaged), per element layer;
  * strain: eps_zz and eps_ss of the same cells, with the free thermal strain
    alpha (T_corner - T_ref). VP1 (T linear in x, y): the section is
    stress-free and eps = alpha theta I away from the clamp and the core.
    VP2 (T not linear in x, y): eps departs from alpha theta I.

Usage (z3st conda env), from the case folder::

    python ../corners.py .

The case post.yaml must have a vp1: or vp2: closed form (T1, T2, PHI_HOT of its
module). Output: ``U_hex_shell/results/IAEA_<case>_corners_z.png``, <case> the folder name.
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, *[os.pardir] * 3)))
import hex_post as hp  # noqa: E402

COLORS = {"hot": "tab:red", "cold": "tab:blue"}


def corner_profiles(case_dir):
    """Per corner: z of the nodes and u_y; z of the cell layers and
    sigma_zz, sigma_ss, sigma_vm, eps_zz, eps_ss, eps_th; plus the section
    mean u_y and the closed-form bow along y."""
    case_dir = os.path.abspath(case_dir)
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    mat = hp.load_material(os.path.join(case_dir, cfg["material"]))
    d = hp.load_z3st(os.path.join(case_dir, cfg["fields"]), temperature=cfg["temperature"], case_dir=case_dir)
    fr = hp.frame(d)
    cf = hp.closed_form(cfg, mat, case_dir, fr)
    res = hp.analyse(d, cfg.get("support", "clamped"), cf)
    X, u, S, E = d["X"], d["u"], d["S"], d["E"]
    C = X[d["cells"]].mean(axis=1)
    r_mid = fr["a_mid"] / np.cos(np.pi / 6)            # mid-wall circumradius

    out = {"case": os.path.basename(case_dir), "z_nodes": res["z_nodes"], "u_y_mean": res["u_bar"][:, 1],
           "u_y_cf": cf.w(res["z_nodes"]) * cf.e[1], "corners": {}}
    for name, phi, T_corner in (("hot", cf.m.PHI_HOT, cf.m.T2), ("cold", cf.m.PHI_HOT + np.pi, cf.m.T1)):
        c = r_mid * np.array([np.cos(phi), np.sin(phi)])
        # nodes: inner and outer corner node of each section
        dn = np.hypot(*(X[:, :2] - c).T)
        sel = dn < 1.01 * np.sort(np.unique(np.round(dn, 9)))[0] + 1e-9
        zn, u_y = hp.level_mean(X[sel, 2], u[sel, 1])
        # cells: the two cells next to the corner in each layer
        dc = np.hypot(*(C[:, :2] - c).T)
        selc = dc < dc.min() * (1 + 1e-6) + 1e-12
        mean = lambda v: hp.level_mean(C[selc, 2], v[selc])[1]
        zc = hp.levels(C[selc, 2])[0]
        T = T_corner(zc)
        out["corners"][name] = {
            "z_nodes": zn, "u_y": u_y, "z_cells": zc,
            "s_zz": mean(S[:, 2, 2]), "s_ss": mean(res["ss"]), "s_vm": mean(res["vm"]),
            "e_zz": mean(E[:, 2, 2]), "e_ss": mean(res["eps_ss"]),
            "e_th": mat["alpha"] * (T - mat["T_ref"]), "T": T,
        }
    return out


def figure_corners(p):
    fig, axes = plt.subplots(1, 3, figsize=(16, 8.5), sharey=True)
    label = {"hot": "hot corner (T2)", "cold": "cold corner (T1)"}

    ax = axes[0]
    for name, c in p["corners"].items():
        ax.plot(c["u_y"] * 1e3, c["z_nodes"] * 1e3, "-", color=COLORS[name], lw=2, label=f"z3st, {label[name]}")
    ax.plot(p["u_y_mean"] * 1e3, p["z_nodes"] * 1e3, "-", color="k", lw=1, alpha=0.6, label="z3st, section mean")
    ax.plot(p["u_y_cf"] * 1e3, p["z_nodes"] * 1e3, "--", color="0.3", lw=1, label="closed form bow w")
    ax.set_xlabel("u_y (mm)")
    ax.set_ylabel("height above support z (mm)")
    ax.set_title("Displacement u_y")

    ax = axes[1]
    for name, c in p["corners"].items():
        z, col = c["z_cells"] * 1e3, COLORS[name]
        ax.plot(c["s_zz"] * 1e-6, z, "-", color=col, lw=1.8, label=f"σ_zz, {label[name]}")
        ax.plot(c["s_ss"] * 1e-6, z, "--", color=col, lw=1.2, label=f"σ_ss (hoop), {label[name]}")
        ax.plot(c["s_vm"] * 1e-6, z, ":", color=col, lw=1.5, label=f"σ_vm, {label[name]}")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("stress (MPa)")
    ax.set_title("Stress (corner cells, DG0)")

    ax = axes[2]
    for name, c in p["corners"].items():
        z, col = c["z_cells"] * 1e3, COLORS[name]
        ax.plot(c["e_zz"] * 1e3, z, "-", color=col, lw=1.8, label=f"ε_zz, {label[name]}")
        ax.plot(c["e_ss"] * 1e3, z, "--", color=col, lw=1.2, label=f"ε_ss (hoop), {label[name]}")
        ax.plot(c["e_th"] * 1e3, z, ":", color="k", lw=1.2,
                label="α (T - T_ref), free thermal" if name == "hot" else None)
    ax.set_xlabel("total strain (10⁻³)")
    ax.set_title("Strain (corner cells, DG0)")

    for ax in axes:
        for zc in (1500, 2500):
            ax.axhline(zc, color="0.6", lw=0.8, ls=":")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper left" if ax is not axes[0] else "lower right")
    axes[0].text(0.02, 1540, "core bottom", fontsize=8, color="0.4", transform=axes[0].get_yaxis_transform())
    axes[0].text(0.02, 2540, "core top", fontsize=8, color="0.4", transform=axes[0].get_yaxis_transform())
    hot = p["corners"]["hot"]["T"].max() - 273.15
    cold = p["corners"]["cold"]["T"].max() - 273.15
    return hp.finish(fig, f"IAEA {p['case']}: hot and cold corner along z",
                     f"corners at mid-wall; T above the core: hot {hot:.1f} C, cold {cold:.1f} C, "
                     "400 C below the core (= T_ref); bottom clamped", top=0.91)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    p = corner_profiles(sys.argv[1] if len(sys.argv) > 1 else ".")
    os.makedirs(hp.RESULTS, exist_ok=True)
    path = os.path.join(hp.RESULTS, f"IAEA_{p['case']}_corners_z.png")
    fig = figure_corners(p)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")
    for name, c in p["corners"].items():
        name = {"hot": "hot corner (T2)", "cold": "cold corner (T1)"}[name]
        print(f"{name}: u_y(top) = {c['u_y'][-1] * 1e3:.3f} mm, max|σ_zz| = {np.abs(c['s_zz']).max() * 1e-6:.1f} MPa "
              f"at z = {c['z_cells'][np.argmax(np.abs(c['s_zz']))] * 1e3:.0f} mm, "
              f"max σ_vm = {c['s_vm'].max() * 1e-6:.1f} MPa at z = {c['z_cells'][np.argmax(c['s_vm'])] * 1e3:.0f} mm")


if __name__ == "__main__":
    main()
