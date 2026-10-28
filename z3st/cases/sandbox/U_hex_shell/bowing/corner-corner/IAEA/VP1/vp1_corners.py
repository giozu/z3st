#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1: displacement u_y, stress and strain along z at the two corners of
extreme temperature, the hot one (T2, direction PHI_HOT of
vp1_temperature) and the cold one opposite (T1).

  * u_y: mid-wall corner (mean of the inner and outer corner node), with the
    z3st section mean and the closed-form bow w(z);
  * stress: sigma_zz, hoop sigma_ss and von Mises of the two DG0 cells next to
    the corner (one per flat, averaged), per element layer;
  * strain: eps_zz and eps_ss of the same cells, with the free thermal strain
    alpha (T_corner - T_ref): for T linear in (x, y) the free solution is
    stress-free, so eps = alpha theta I away from the clamp and the core.

Usage (z3st conda env)::

    python vp1_corners.py [--xdmf output/<run>.xdmf]

Output: ``U_hex_shell/results/IAEA_VP1_corners_z.png``.
"""

import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 4))
sys.path.insert(0, U_HEX)
sys.path.insert(0, HERE)
import hex_post as hp  # noqa: E402
import vp1_temperature as vt  # noqa: E402

TITLE = "IAEA VP1: hot and cold corner along z"
CORNERS = {"hot corner (T2)": ("tab:red", vt.PHI_HOT), "cold corner (T1)": ("tab:blue", vt.PHI_HOT + np.pi)}


def corner_profiles(xdmf):
    """Per corner: z of the nodes and u_y; z of the cell layers and
    sigma_zz, sigma_ss, sigma_vm, eps_zz, eps_ss, eps_th; plus the section
    mean u_y and the closed-form bow along y."""
    with open(os.path.join(HERE, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    mat = hp.load_material(os.path.join(HERE, cfg["material"]))
    d = hp.load_z3st(xdmf, temperature=cfg["temperature"], case_dir=HERE)
    fr = hp.frame(d)
    cf = hp.closed_form(cfg, mat, HERE, fr)
    res = hp.analyse(d, mat, cfg.get("support", "clamped"), cf)
    X, u, S, E = d["X"], d["u"], d["S"], d["E"]
    C = X[d["cells"]].mean(axis=1)
    r_mid = fr["a_mid"] / np.cos(np.pi / 6)            # mid-wall circumradius

    out = {"z_nodes": res["z_nodes"], "u_y_mean": res["u_bar"][:, 1],
           "u_y_cf": cf.w(res["z_nodes"]) * cf.e[1], "corners": {}}
    for name, (_, phi) in CORNERS.items():
        c = r_mid * np.array([np.cos(phi), np.sin(phi)])
        # nodes: inner and outer corner node of each section
        dn = np.hypot(*(X[:, :2] - c).T)
        sel = dn < 1.01 * np.sort(np.unique(np.round(dn, 9)))[0] + 1e-9
        zn, inv = np.unique(np.round(X[sel, 2], 9), return_inverse=True)
        u_y = np.bincount(inv, weights=u[sel, 1]) / np.bincount(inv)
        # cells: the two cells next to the corner in each layer
        dc = np.hypot(*(C[:, :2] - c).T)
        selc = dc < dc.min() * (1 + 1e-6) + 1e-12
        zc, invc = np.unique(np.round(C[selc, 2], 9), return_inverse=True)
        mean = lambda v: np.bincount(invc, weights=v[selc]) / np.bincount(invc)
        T_corner = np.ravel(vt.T_vp1(np.array([np.full_like(zc, c[0]), np.full_like(zc, c[1]), zc])))
        out["corners"][name] = {
            "z_nodes": zn, "u_y": u_y, "z_cells": zc,
            "s_zz": mean(S[:, 2, 2]), "s_ss": mean(res["ss"]), "s_vm": mean(res["vm"]),
            "e_zz": mean(E[:, 2, 2]), "e_ss": mean(res["eps_ss"]),
            "e_th": mat["alpha"] * (T_corner - mat["T_ref"]), "T": T_corner,
        }
    return out


def figure_corners(p):
    fig, axes = plt.subplots(1, 3, figsize=(16, 8.5), sharey=True)

    ax = axes[0]
    for name, (col, _) in CORNERS.items():
        c = p["corners"][name]
        ax.plot(c["u_y"] * 1e3, c["z_nodes"] * 1e3, "-", color=col, lw=2, label=f"z3st, {name}")
    ax.plot(p["u_y_mean"] * 1e3, p["z_nodes"] * 1e3, "-", color="k", lw=1, alpha=0.6, label="z3st, section mean")
    ax.plot(p["u_y_cf"] * 1e3, p["z_nodes"] * 1e3, "--", color="0.3", lw=1, label="closed form bow w")
    ax.set_xlabel("u_y (mm)")
    ax.set_ylabel("height above support z (mm)")
    ax.set_title("Displacement u_y")

    ax = axes[1]
    for name, (col, _) in CORNERS.items():
        c = p["corners"][name]
        z = c["z_cells"] * 1e3
        ax.plot(c["s_zz"] * 1e-6, z, "-", color=col, lw=1.8, label=f"σ_zz, {name}")
        ax.plot(c["s_ss"] * 1e-6, z, "--", color=col, lw=1.2, label=f"σ_ss (hoop), {name}")
        ax.plot(c["s_vm"] * 1e-6, z, ":", color=col, lw=1.5, label=f"σ_vm, {name}")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("stress (MPa)")
    ax.set_title("Stress (corner cells, DG0)")

    ax = axes[2]
    for name, (col, _) in CORNERS.items():
        c = p["corners"][name]
        z = c["z_cells"] * 1e3
        ax.plot(c["e_zz"] * 1e3, z, "-", color=col, lw=1.8, label=f"ε_zz, {name}")
        ax.plot(c["e_ss"] * 1e3, z, "--", color=col, lw=1.2, label=f"ε_ss (hoop), {name}")
        ax.plot(c["e_th"] * 1e3, z, ":", color="k", lw=1.2,
                label="α (T - T_ref), free thermal" if col == "tab:red" else None)
    ax.set_xlabel("total strain (10⁻³)")
    ax.set_title("Strain (corner cells, DG0)")

    for ax in axes:
        for zc in (1500, 2500):
            ax.axhline(zc, color="0.6", lw=0.8, ls=":")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="upper left" if ax is not axes[0] else "lower right")
    axes[0].text(0.02, 1540, "core bottom", fontsize=8, color="0.4", transform=axes[0].get_yaxis_transform())
    axes[0].text(0.02, 2540, "core top", fontsize=8, color="0.4", transform=axes[0].get_yaxis_transform())
    hot = p["corners"]["hot corner (T2)"]["T"].max() - 273.15
    cold = p["corners"]["cold corner (T1)"]["T"].max() - 273.15
    return hp._finish(fig, TITLE, f"corners at mid-wall; T above the core: hot {hot:.1f} C, cold {cold:.1f} C, "
                      "400 C below the core (= T_ref); bottom clamped", top=0.91)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    with open(os.path.join(HERE, "post.yaml")) as fh:
        default = os.path.join(HERE, yaml.safe_load(fh)["fields"])
    ap.add_argument("--xdmf", default=default)
    ap.add_argument("--out", default=hp.RESULTS)
    args = ap.parse_args()
    p = corner_profiles(args.xdmf)
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, "IAEA_VP1_corners_z.png")
    fig = figure_corners(p)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")
    for name, c in p["corners"].items():
        print(f"{name}: u_y(top) = {c['u_y'][-1] * 1e3:.3f} mm, max|σ_zz| = {np.abs(c['s_zz']).max() * 1e-6:.1f} MPa "
              f"at z = {c['z_cells'][np.argmax(np.abs(c['s_zz']))] * 1e3:.0f} mm, "
              f"max σ_vm = {c['s_vm'].max() * 1e-6:.1f} MPa at z = {c['z_cells'][np.argmax(c['s_vm'])] * 1e3:.0f} mm")


if __name__ == "__main__":
    main()
