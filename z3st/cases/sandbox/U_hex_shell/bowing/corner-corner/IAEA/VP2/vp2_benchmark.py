#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP2 (benchmark N1, example 2, non-linear thermal gradient, corner to
corner): z3st against the published results.

Reference data (SOURCES), IWGFR/75 (Part I) and IWGFR/76 (Part II):
  * ../iwgfr.py: IWGFR/75 linear theory (DT1BOW, DT1LEN), the AT1BOW
    AVERAGE(10) at the URP, Y of the 11 participants (DT1BOW), the outer
    across-flats increase (SISCO) and the sensitivity studies of Part II
    (ARKAS shell, HIBEACON standard: other T);
  * S1 (CEA, SOLO), Part II: Table 2 "S1 deflection", Y and Z at core top,
    LRP and URP;
  * BelgoNucleaire (STRAW), Part II: Table III "Results of example 2 -
    reference case" (Y, Z and the mid-wall dimensions across corners OY and
    across flats OX), Table IV "zero thickness" (Y, Z) and "non adiabatic
    wrapper" (Y, Z, dimensions). The non adiabatic wrapper has another
    temperature field: listed for reference, not the z3st case.

No bow profile along z is published for example 2 (the S1 Fig. 2 used for
VP1 is example 1 only): the bow is compared at the three stations.

The axial values (Z) and the dimensions depend on the stress-free
temperature: STRAW and the linear theory are given for T_ref = 0 / 20 C, S1
Table 2 matches about 20 C, while the case uses T_ref = 400 C. The z3st values
are rebased to the T_ref of each source by superposing the free dilatation
alpha (T_ref,case - T_ref,source), valid away from the clamp: u_z + alpha dT z
and L + alpha dT L0. The bow does not depend on T_ref.

Usage (z3st conda env)::

    python vp2_benchmark.py

Outputs in ``U_hex_shell/results/``, named ``IAEA_VP2_vs_benchmark*``: the
overview figure (.png), the full comparison table (_tables.png) and the
summary (.json). A point is flagged when |z3st - benchmark| exceeds both
benchmark.TOL_MM and TOL_REL of the benchmark value.
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import benchmark as bc  # noqa: E402
import iwgfr  # noqa: E402
from benchmark import Z_ABOVE, Z_BELOW  # noqa: E402
hp = bc.hp

PREFIX, TITLE = "IAEA_VP2_vs_benchmark", "IAEA VP2: z3st vs IAEA benchmark N1, example 2"

# Per source: T_ref (C) of its Z and dimensions, Y / Z per station (mm),
# dimensions per level (across corners OY, across flats OX, mid-wall, mm),
# plot style and whether it is the same problem as the z3st case.
SOURCES = {
    "IWGFR/75 theory": {
        "label": "IWGFR/75 AT1BOW, linear theory", "T_ref_C": 20.0, "same_case": True, "marker": "s",
        "color": "tab:red",
        "Y": {"Core top": 1.333, "LRP": 4.333, "URP": 16.333},
        "Z": {}, "dims": {},
    },
    "S1": {
        "label": "S1 (CEA), Table 2", "T_ref_C": 20.0, "same_case": True, "marker": "D", "color": "tab:green",
        "Y": {"Core top": 1.333, "LRP": 4.338, "URP": 16.334},
        "Z": {"Core top": 18.338, "LRP": 22.882, "URP": 31.973},
        "dims": {},
    },
    "Table III": {
        "label": "STRAW Table III, reference", "T_ref_C": 0.0, "same_case": True, "marker": "o",
        "color": "tab:purple",
        "Y": {"Core top": 1.333, "LRP": 4.332, "URP": 16.328},
        "Z": {"Core top": 19.125, "LRP": 23.850, "URP": 33.300},
        "dims": {"hot, below core bottom": (151.076, 130.840),
                 "hot, above core top": (151.412, 131.157),
                 "cold": (150.000, 129.904)},
    },
    "Table IV, zero thickness": {
        "label": "STRAW Table IV, zero thickness", "T_ref_C": 0.0, "same_case": True, "marker": "^",
        "color": "tab:brown",
        "Y": {"Core top": 1.330, "LRP": 4.333, "URP": 16.333},
        "Z": {"Core top": 19.125, "LRP": 23.850, "URP": 33.300},
        "dims": {},
    },
    "Table IV, non adiabatic": {
        "label": "STRAW Table IV, non adiabatic", "T_ref_C": 0.0, "same_case": False, "marker": "v",
        "color": "tab:gray",
        "Y": {"Core top": 1.334, "LRP": 4.336, "URP": 16.345},
        "Z": {"Core top": 19.125, "LRP": 23.850, "URP": 33.300},
        "dims": {"hot, below core bottom": (151.076, 130.840),
                 "hot, above core top": (151.415, 131.080),
                 "cold": (150.000, 129.904)},
    },
}

# IWGFR/75 theory and declared average, the participants (DT1BOW, AT1LEN) and the
# sensitivity studies of Part II: data in ../iwgfr.py, example 2
EX = "ex2"
SOURCES["IWGFR/75 theory"]["Z"] = dict(zip(iwgfr.STATIONS, iwgfr.AXIAL_THEORY))   # DT1LEN, same for both examples
SOURCES["IWGFR/75 AVERAGE(10)"] = {
    "label": "IWGFR/75 AT1BOW, AVERAGE(10)", "T_ref_C": 20.0, "same_case": True, "marker": "*", "color": "tab:orange",
    "Y": {"URP": iwgfr.DECLARED_URP[EX]["AVERAGE(10)"][0]}, "Z": {}, "dims": {},
}
for _name, _c in {**iwgfr.CODES, **iwgfr.SENSITIVITY}.items():
    _code = _name in iwgfr.CODES
    SOURCES[_name] = {
        "label": "IWGFR/75 DT1BOW, participants" if _code else "IWGFR/76, sensitivity studies",
        "T_ref_C": 20.0, "same_case": _c.get("same_case", True), "marker": "x" if _code else "+", "color": "0.45",
        "Y": dict(zip(iwgfr.STATIONS, _c[EX])), "Z": {}, "dims": {},
        "dH": {"hot, above core top": iwgfr.LATERAL[EX][_name]} if _name in iwgfr.LATERAL[EX] else {},
    }

# Interpretation of a benchmark value, shown next to it
NOTES = {
    ("IWGFR/75 theory", None, None): "beam theory, thin wall: 4/3 of VP1 (VP2/VP1 ratio 1.333 in AT1BOW)",
    ("S1", "Z", None): "S1 Z matches T_ref = 20 C up to a constant, as in VP1",
    ("Table IV, non adiabatic", None, None): "other T field (sensitivity study), not the z3st case",
}
NOTES.update({(_n, None, None): _c["note"] for _n, _c in iwgfr.SENSITIVITY.items()})
NOTES[(iwgfr.REJECTED, None, None)] = "out of the IAEA AVERAGE(10) (IWGFR/75 para 19)"
NOTES[("IWGFR/75 AVERAGE(10)", None, None)] = f"IAEA declared, spread {iwgfr.DECLARED_URP[EX]['AVERAGE(10)'][1]} %"
NOTES[("IWGFR/75 theory", "Z", None)] = "DT1LEN, average temperature inlet -> station"


def note(src, qty, station):
    return bc.note(NOTES, src, qty, station)


def z3st_profiles(xdmf=None):
    """Section-mean bow, u_z, dimensions and closed forms of this run (benchmark.profiles)."""
    return bc.profiles(HERE, xdmf)


def shift_per_m(p, T_ref_C):
    return bc.shift_per_m(p, T_ref_C)


def compare(p, tol_mm=bc.TOL_MM, tol_rel=bc.TOL_REL):
    """benchmark.compare of SOURCES."""
    return bc.compare(SOURCES, NOTES, p, tol_mm, tol_rel)


def figure(p, rows, tol_mm, tol_rel):
    n = max(sum(r["flagged"] for r in rows), 1)
    fig = plt.figure(figsize=(16, 9 + 0.3 * n))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.2, 0.5 + 0.3 * n], width_ratios=[1.25, 1, 1])
    zmm = p["z"] * 1e3

    # bow
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(p["w"] * 1e3, zmm, "-", color="tab:blue", lw=2, label="z3st")
    ax.plot(p["w_cf"] * 1e3, zmm, "--", color="k", lw=1,
            label=rf"closed form $\kappa={p['C']:.3f}\,\alpha\,\Delta T/D$")
    bc.source_markers(ax, SOURCES, "Y")
    for zc, lab in ((1500, "core bottom"), (2500, "core top")):
        ax.axhline(zc, color="0.6", lw=0.8, ls=":")
        ax.text(0.2, zc + 30, lab, fontsize=8, color="0.4")
    ax.set_xlabel("bow w towards the cold corner (mm)")
    ax.set_ylabel("height above support (mm)")
    ax.set_title("Thermal bow (Y)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")

    # bow residual
    ax = fig.add_subplot(gs[0, 1], sharey=fig.axes[0])
    ax.plot((p["w"] - p["w_cf"]) * 1e3, zmm, "-", color="k", lw=1, label="z3st - closed form")
    bc.source_markers(ax, SOURCES, "Y", p, residual=True)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("z3st - benchmark (mm)")
    ax.set_title("Bow: z3st - benchmark")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower right")

    # axial
    bc.axial_panel(fig.add_subplot(gs[0, 2], sharey=fig.axes[0]), p, SOURCES)

    # flagged points
    bc.flagged_table(fig.add_subplot(gs[1, :]), rows, tol_mm, tol_rel)
    return hp.finish(fig, TITLE, "bottom clamped at z = 0, core 1500-2500 mm, corners 1, 2, 6 / 3, 4, 5 at "
                      "500 / 550 C above the core; Z rebased to the T_ref of each source (S1 20 C, STRAW 0 C)",
                      top=0.91)


def figure_tables(rows, tol_mm, tol_rel):
    """Every benchmark value against z3st, flagged rows highlighted."""
    fig = plt.figure(figsize=(17, 1.6 + 0.32 * (len(rows) + 1)))
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.84])
    bc.put_table(ax, rows, highlight=[r["flagged"] for r in rows])
    return hp.finish(fig, TITLE, f"IWGFR/75 (../iwgfr.py), S1 Table 2 and STRAW Tables III / IV against z3st; "
                      f"highlighted: |diff| > {tol_mm} mm and > {100 * tol_rel:.0f} %. Dimensions at mid-wall, "
                      f"z = {Z_BELOW * 1e3:.0f} / {Z_ABOVE * 1e3:.0f} mm, rebased to T_ref = 0 C; "
                      "cold = undeformed mesh", top=0.88)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    p = z3st_profiles()
    rows = compare(p)
    print(f"[INFO] section factor C = {p['C']:.5f} (thin-wall theory 4/3)")
    bc.write(PREFIX, p, rows, {"": figure(p, rows, bc.TOL_MM, bc.TOL_REL),
                               "_tables": figure_tables(rows, bc.TOL_MM, bc.TOL_REL)}, bc.TOL_MM, bc.TOL_REL)


if __name__ == "__main__":
    main()
