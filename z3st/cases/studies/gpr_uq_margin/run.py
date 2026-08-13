#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: propagation of the GPR conductivity uncertainty to the melting margin
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Solve the same high-rated MA-MOX pin at several offsets of the GPR posterior.

The GPR conductivity hook keeps the posterior standard deviation of its log
residual, so a scenario solve can be run at k = k_Magni * exp(mean + xi*sigma).
This script sweeps xi and reports the resulting centre temperature and the
margin to the MOX solidus, i.e. what the uncertainty of a data-driven material
law does to a safety-relevant coupled outcome.

Run:  python3 run.py     (after ../magni_gpr_conductivity/make_synthetic_gpr.py)
"""

import glob
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.plotstyle import apply as _apply_plotstyle
_apply_plotstyle()
import numpy as np
import pyvista as pv
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
XI = [-2.0, -1.0, 0.0, 1.0, 2.0]

# MOX solidus, Adamson et al., J. Nucl. Mater. 130 (1985) 349:
#   T_sol = 3120 - 655.3 y + 336.4 y^2 - 99.9 y^3   (y = Pu/(U+Pu))
def solidus(y):
    return 3120.0 - 655.3 * y + 336.4 * y ** 2 - 99.9 * y ** 3


def run_one(xi):
    card = yaml.safe_load(open(os.path.join(HERE, "fuel.yaml")))
    card["k"]["mode"] = "affine"
    card["k"]["xi"] = xi
    with open(os.path.join(HERE, "fuel.yaml"), "w") as f:
        yaml.safe_dump(card, f, sort_keys=False)

    for f in glob.glob(os.path.join(OUT, "*.vtu")):
        os.remove(f)
    subprocess.run([sys.executable, "-m", "z3st"], cwd=HERE, check=True,
                   stdout=open(os.path.join(HERE, "log_z3st.md"), "w"),
                   stderr=subprocess.STDOUT)

    grid = pv.read(sorted(glob.glob(os.path.join(OUT, "*.vtu")))[-1])
    T = grid.point_data["Temperature"]
    return float(T.max()), card["Pu"]


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    temps = []
    for xi in XI:
        Tc, Pu = run_one(xi)
        temps.append(Tc)
        print(f"  xi = {xi:+.1f} sigma   T_centre = {Tc:8.1f} K")
    temps = np.array(temps)
    Tsol = solidus(Pu)
    margin = Tsol - temps
    print(f"\n  MOX solidus (Pu = {Pu:.2f}): {Tsol:.0f} K")
    for xi, Tc, m in zip(XI, temps, margin):
        print(f"  xi = {xi:+.1f}  T_centre = {Tc:7.1f} K   margin = {m:6.1f} K")
    span = margin[0] - margin[-1]
    print(f"\n  margin spread over +/-2 sigma: {abs(span):.0f} K "
          f"({abs(span)/margin[XI.index(0.0)]*100:.0f}% of the nominal margin)")

    np.savez(os.path.join(OUT, "uq_margin.npz"), xi=XI, T=temps,
             solidus=Tsol, margin=margin)

    plt.rcParams.update({"font.size": 13, "axes.labelsize": 14,
                         "xtick.labelsize": 12, "ytick.labelsize": 12,
                         "legend.fontsize": 11})
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(XI, temps, "o-", color="#0072B2", lw=2, label="centre temperature")
    ax.axhline(Tsol, color="#D55E00", ls="--", lw=2,
               label=f"MOX solidus, {Tsol:.0f} K")
    ax.fill_between(XI, temps, Tsol, color="#0072B2", alpha=0.10)
    ax.set_xlabel(r"GPR posterior offset $\xi$ (standard deviations)")
    ax.set_ylabel("temperature (K)")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "uq_margin.png"), dpi=150)
    print("  wrote output/uq_margin.png")
