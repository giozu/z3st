#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 (z-gradient): z3st with the T-dependent 15-15Ti card (this folder)
against z3st with the constant card (../), same mesh, BCs and T, with the
outputs of ``U_hex_shell/hex_post.py`` overlaid. Each run is compared with its
own beam: kappa = alpha dT/D for the constant card, the E(T)-weighted fit of
eta(T) for this one (hex_post VP1 / VP1Tdep).

This folder from ``post.yaml``; ../ from its ``input.yaml`` (output and card),
the rest as here. Thermal model on in both: T from the fields file, checked
equal on the nodes (T_TOL).

Usage (z3st conda env)::

    python compare_const.py

Outputs in ``U_hex_shell/results/bowing_corner-corner_IAEA_VP1/1515Ti_Tdep_vs_const/``,
named ``IAEA_VP1_Tdep_vs_const_<figure>``: the hex_post figures (along s with the
hot/cold band of results_organize.py), a table Tdep vs constant
(``differences``) and the summary (.json).
"""

import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
CONST = os.path.dirname(HERE)                                   # VP1/, constant card
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 5))   # U_hex_shell, home of hex_post.py
sys.path.insert(0, U_HEX)
import hex_post as hp  # noqa: E402
import results_organize as ro  # noqa: E402

OUT = os.path.join(hp.RESULTS, "bowing_corner-corner_IAEA_VP1", "1515Ti_Tdep_vs_const")
PREFIX = "IAEA_VP1_Tdep_vs_const_"
TITLE = "IAEA_VP1 (z-gradient): 1515Ti_Tdep vs constant card"
T_TOL = 1e-6   # (K) max |T Tdep - T constant| on the nodes: above it the runs are not the same case
ABOVE = (2.6, 3.9)   # (m) above the core, away from the kink and the top: stress-free with the constant card


def run(case_dir, cfg):
    mat = hp.load_material(os.path.join(case_dir, cfg["material"]))
    d = hp.load_z3st(os.path.join(case_dir, cfg["fields"]), temperature=cfg["temperature"], case_dir=case_dir)
    return hp.analyse(d, mat, cfg["support"], hp.closed_form(cfg, mat, case_dir, hp.frame(d)))


def diff_rows(runs, zs):
    """Tdep and constant side by side: bow, curvature, u_z, largest stresses."""
    (lt, t), (lc, c) = runs
    rows = []

    def row(name, get, fmt="{:.4g}"):
        a, b = get(t), get(c)
        rows.append([name, fmt.format(b), fmt.format(a), fmt.format(a - b), f"{a / b:.4f}" if b else "-"])

    for z in zs:
        row(f"w at z = {z * 1e3:.0f} mm (mm)", lambda r, z=z: np.interp(z, r["z_nodes"], r["w"]) * 1e3)
    row("κ at z = 3000 mm (1/m)", lambda r: np.interp(3.0, r["z_sec"], r["kappa"]))
    row("u_z at the top (mm)", lambda r: r["u_bar"][-1, 2] * 1e3)
    row("max σ_vm (MPa)", lambda r: r["vm"].max() * 1e-6)
    zc = t["frame"]["cells"]["z"]
    above = (zc > ABOVE[0]) & (zc < ABOVE[1])
    row(f"max σ_vm, {ABOVE[0] * 1e3:.0f} < z < {ABOVE[1] * 1e3:.0f} mm (MPa)", lambda r: r["vm"][above].max() * 1e-6)
    row("max |σ_zz| (MPa)", lambda r: np.abs(r["data"]["S"][:, 2, 2]).max() * 1e-6)
    row("max |σ_ss| (MPa)", lambda r: np.abs(r["ss"]).max() * 1e-6)
    return rows


def fig_differences(runs, zs, title):
    rows = diff_rows(runs, zs)
    fig, ax = plt.subplots(figsize=(13, 1.3 + 0.36 * (len(rows) + 1)))
    hp._table(ax, ["quantity", runs[1][0], runs[0][0], f"{runs[0][0]} − {runs[1][0]}",
                   f"{runs[0][0]} / {runs[1][0]}"], rows, fill=True)
    return hp._finish(fig, title, "Same mesh, BCs and T; bow w and u_z: section mean of the nodes; "
                      "stresses: cell values (DG0)", top=0.86)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    with open(os.path.join(HERE, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    with open(os.path.join(CONST, "input.yaml")) as fh:
        inp = yaml.safe_load(fh)
    cfg_c = dict(cfg, fields=f"output/{inp['output']['filename']}.xdmf", material=inp["materials"]["steel"],
                 closed_form="vp1:vp1_temperature")
    runs = [("1515Ti_Tdep", run(HERE, cfg)), ("constant", run(CONST, cfg_c))]   # first: along-s closed form
    dT = np.abs(runs[0][1]["data"]["T"] - runs[1][1]["data"]["T"]).max()
    if dT > T_TOL:
        sys.exit(f"T of {HERE} and {CONST} differ by up to {dT:.2e} K (> {T_TOL} K): not the same case")

    levels = cfg["along_s_z"]
    summary = hp.save_all(runs, OUT, PREFIX, TITLE, cfg.get("stations_z"), vtb_3d=False, levels=levels)
    res = runs[0][1]
    notes = "\n".join(ro.hot_cold_text(res, z, ro.thermal_sides(res, z), tag) for tag, z in levels.items())
    for group, keys in hp.ALONG_S_GROUPS:
        ro.annotate(os.path.join(OUT, f"{PREFIX}along_s_{group}.png"), notes,
                    [f"{tag}: {x}" for tag, z in levels.items() for x in ro.along_s_flags(res, z, keys)])
    zs = [2.5, 3.0, 4.0]   # core top, ACLP, duct top
    path = os.path.join(OUT, f"{PREFIX}differences.png")
    fig = fig_differences(runs, zs, TITLE)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")
    summary["sources"] = {lab: os.path.join(d, c["fields"]) for (lab, _), d, c in zip(runs, (HERE, CONST), (cfg, cfg_c))}
    summary["differences"] = diff_rows(runs, zs)
    with open(os.path.join(OUT, f"{PREFIX}summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
