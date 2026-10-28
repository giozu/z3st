#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 (z-gradient), T-dependent 15-15Ti card: z3st with the Di Marcello set
(DiMarcello/) against z3st with the Sobolev set (this folder), same mesh, BCs
and T, with the outputs of ``U_hex_shell/hex_post.py`` overlaid. Each run is
compared with its own beam, the E(T)-weighted fit of eta(T) of its card
(hex_post VP1Tdep). Functions of compare_const.py.

Both folders from their ``post.yaml``; T from the fields files, checked equal
on the nodes (T_TOL).

Usage (z3st conda env)::

    python compare_laws.py

Outputs in ``U_hex_shell/results/steel_Tdependent/``, named
``IAEA_VP1_DiMarcello_vs_Sobolev_<figure>``: the hex_post figures (along s with
the hot/cold band of results_organize.py), a table Di Marcello vs Sobolev
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
sys.path.insert(0, HERE)
from compare_const import T_TOL, diff_rows, fig_differences, hp, ro, run  # noqa: E402

DIM = os.path.join(HERE, "DiMarcello")                          # Di Marcello set
OUT = os.path.join(hp.RESULTS, "steel_Tdependent")
PREFIX = "IAEA_VP1_DiMarcello_vs_Sobolev_"
TITLE = "IAEA_VP1 (z-gradient), 1515Ti_Tdep: Di Marcello vs Sobolev"


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    cfgs = {}
    for d in (DIM, HERE):
        with open(os.path.join(d, "post.yaml")) as fh:
            cfgs[d] = yaml.safe_load(fh)
    runs = [("DiMarcello", run(DIM, cfgs[DIM])), ("Sobolev", run(HERE, cfgs[HERE]))]
    dT = np.abs(runs[0][1]["data"]["T"] - runs[1][1]["data"]["T"]).max()
    if dT > T_TOL:
        sys.exit(f"T of {DIM} and {HERE} differ by up to {dT:.2e} K (> {T_TOL} K): not the same case")

    os.makedirs(OUT, exist_ok=True)
    cfg = cfgs[HERE]
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
    summary["sources"] = {lab: os.path.join(d, cfgs[d]["fields"]) for (lab, _), d in zip(runs, (DIM, HERE))}
    summary["differences"] = diff_rows(runs, zs)
    with open(os.path.join(OUT, f"{PREFIX}summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
