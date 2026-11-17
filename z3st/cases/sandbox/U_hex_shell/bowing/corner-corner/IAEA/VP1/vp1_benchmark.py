#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 (benchmark N1, example 1, corner to corner): z3st against the
published results.

Reference data (SOURCES):
  * S1: Fig. 2 "S1 thermal bow", deflection DT (mm) against the height above
    the support (mm), digitized from the scanned page (about +-0.15 mm), and
    Table 1 "S1 deflection", Y (lateral bow) and Z (axial) at core top, LRP
    and URP;
  * Table I "Results of example 1 - reference case": Y and Z deflections, and
    the mid-wall dimensions across corners (OY) and across flats (OX) below
    the core bottom, above the core top and cold;
  * Table II "Sensitivity studies": zero thickness (Y, Z) and non adiabatic
    wrapper (Y, Z, dimensions). The non adiabatic wrapper has another
    temperature field: listed for reference, not the z3st case;
  * ../iwgfr.py: IWGFR/75 linear theory (DT1BOW, DT1LEN), the AT1BOW
    AVERAGE(10) at the URP, Y of the 11 participants (DT1BOW), their Z at the
    URP (AT1LEN), the outer across-flats increase (SISCO, SOLO) and the
    sensitivity studies of Part II (ARKAS shell, HIBEACON standard: other T).

The axial values (Z) and the dimensions depend on the stress-free
temperature: Table I and II match T_ref = 0 C, Table 1 of S1 about 20 C,
while the case uses T_ref = 400 C (VP1 stress_free_temperature). The z3st
values are rebased to the T_ref of each source by superposing the free
dilatation alpha (T_ref,case - T_ref,source), valid away from the clamp:
u_z + alpha dT z and L + alpha dT L0.

The dimensions are measured on the z3st mesh at Z_BELOW (below the core,
away from the clamp) and Z_ABOVE (above the core top): mid-wall point = mean
of the inner and outer node of the corner (OY) or of the mid-flat (OX).

Usage (z3st conda env)::

    python vp1_benchmark.py

Outputs in ``U_hex_shell/results/``, named ``IAEA_VP1_vs_benchmark*``: the
overview figure (.png), the bow along z in the layout of Fig. 2 (_fig2.png),
the full comparison with all the sources (_tables.png) and the summary (.json).
A point is flagged when |z3st - benchmark| exceeds both benchmark.TOL_MM and
TOL_REL of the benchmark value.
"""

import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import benchmark as bc  # noqa: E402
import iwgfr  # noqa: E402
from benchmark import Z_ABOVE, Z_BELOW  # noqa: E402
hp = bc.hp

PREFIX, TITLE = "IAEA_VP1_vs_benchmark", "IAEA VP1: z3st vs IAEA benchmark N1, example 1"

# Fig. 2 (S1), digitized: (z above support (mm), DT (mm))
FIG2 = np.array([
    (1500, 0.00), (1600, 0.00), (1700, 0.00), (1800, 0.02), (1900, 0.05),
    (2000, 0.10), (2100, 0.18), (2200, 0.30), (2300, 0.50), (2400, 0.70),
    (2500, 0.95), (2700, 1.65), (2900, 2.60), (3100, 3.85), (3300, 5.35),
    (3500, 7.10), (3700, 8.95), (3900, 11.20), (4000, 12.35),
])
FIG2_ERR = 0.15   # (mm) reading uncertainty of the scan

# Per source: T_ref (C) of its Z and dimensions, Y / Z per station (mm),
# dimensions per level (across corners OY, across flats OX, mid-wall, mm),
# plot style and whether it is the same problem as the z3st case.
SOURCES = {
    "S1": {
        "label": "S1, Table 1", "T_ref_C": 20.0, "same_case": True, "marker": "D", "color": "tab:green",
        "Y": {"Core top": 0.125, "LRP": 3.125, "URP": 12.250},
        "Z": {"Core top": 18.338, "LRP": 22.882, "URP": 31.973},
        "dims": {},
    },
    "Table I": {
        "label": "Table I, reference", "T_ref_C": 0.0, "same_case": True, "marker": "o", "color": "tab:purple",
        "Y": {"Core top": 1.000, "LRP": 3.249, "URP": 12.245},
        "Z": {"Core top": 19.125, "LRP": 23.850, "URP": 33.300},
        "dims": {"hot, below core bottom": (151.076, 130.840),
                 "hot, above core top": (151.412, 131.132),
                 "cold": (150.000, 129.904)},
    },
    "Table II, zero thickness": {
        "label": "Table II, zero thickness", "T_ref_C": 0.0, "same_case": True, "marker": "^", "color": "tab:brown",
        "Y": {"Core top": 1.000, "LRP": 3.250, "URP": 12.250},
        "Z": {"Core top": 19.125, "LRP": 23.850, "URP": 33.300},
        "dims": {},
    },
    "Table II, non adiabatic": {
        "label": "Table II, non adiabatic", "T_ref_C": 0.0, "same_case": False, "marker": "v", "color": "tab:gray",
        "Y": {"Core top": 1.002, "LRP": 3.257, "URP": 12.275},
        "Z": {"Core top": 19.125, "LRP": 23.850, "URP": 33.300},
        "dims": {"hot, below core bottom": (151.076, 130.840),
                 "hot, above core top": (151.400, 131.146),
                 "cold": (150.000, 129.904)},
    },
}

# IWGFR/75 theory and declared average, the participants (DT1BOW, AT1LEN) and the
# sensitivity studies of Part II: data in ../iwgfr.py, example 1
EX = "ex1"
SOURCES["IWGFR/75 theory"] = {
    "label": "IWGFR/75 DT1BOW / DT1LEN, linear theory", "T_ref_C": 20.0, "same_case": True, "marker": "s",
    "color": "tab:red", "Y": dict(zip(iwgfr.STATIONS, iwgfr.THEORY[EX])),
    "Z": dict(zip(iwgfr.STATIONS, iwgfr.AXIAL_THEORY)), "dims": {},
}
SOURCES["IWGFR/75 AVERAGE(10)"] = {
    "label": "IWGFR/75 AT1BOW, AVERAGE(10)", "T_ref_C": 20.0, "same_case": True, "marker": "*", "color": "tab:orange",
    "Y": {"URP": iwgfr.DECLARED_URP[EX]["AVERAGE(10)"][0]}, "Z": {}, "dims": {},
}
for _name, _c in {**iwgfr.CODES, **iwgfr.SENSITIVITY}.items():
    _code = _name in iwgfr.CODES
    SOURCES[_name] = {
        "label": "IWGFR/75 DT1BOW, participants" if _code else "IWGFR/76, sensitivity studies",
        "T_ref_C": 20.0, "same_case": _c.get("same_case", True), "marker": "x" if _code else "+", "color": "0.45",
        "Y": dict(zip(iwgfr.STATIONS, _c[EX])), "Z": {"URP": iwgfr.AXIAL_URP[_name]} if _name in iwgfr.AXIAL_URP else {}, "dims": {},
        "dH": {"hot, above core top": iwgfr.LATERAL[EX][_name]} if _name in iwgfr.LATERAL[EX] else {},
    }

# Interpretation of a benchmark value, shown next to it
NOTES = {
    ("S1", "Y", "Core top"): "misprint: Table I and DT1BOW give 1.000 (0.125 = w at z = 2000 mm)",
    ("S1", "Z", None): "S1 Z matches T_ref = 20 C up to a constant +0.11 mm",
    ("Table II, non adiabatic", None, None): "other T field (sensitivity study), not the z3st case",
    ("S1", "Y", "LRP"): "misprint: DT1BOW gives 3.250",
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
    """benchmark.compare of SOURCES, and the bow at the digitized points of Fig. 2."""
    rows = bc.compare(SOURCES, NOTES, p, tol_mm, tol_rel)
    for z_mm, dt in FIG2:
        w = float(np.interp(z_mm / 1e3, p["z"], p["w"]) * 1e3)
        rows.append(bc.row("S1", "DT (Fig. 2)", "", z_mm / 1e3, w, dt, note("S1", "DT", None), tol_mm, tol_rel))
    return rows


def figure(p, rows, tol_mm, tol_rel):
    n = max(sum(r["flagged"] for r in rows), 1)
    fig = plt.figure(figsize=(16, 9 + 0.3 * n))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.2, 0.5 + 0.3 * n], width_ratios=[1.25, 1, 1])
    zmm = p["z"] * 1e3

    # bow
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(p["w"] * 1e3, zmm, "-", color="tab:blue", lw=2, label="z3st")
    ax.plot(p["w_cf"] * 1e3, zmm, "--", color="k", lw=1, label=r"closed form $\kappa=\alpha\,\Delta T/D$")
    ax.errorbar(FIG2[:, 1], FIG2[:, 0], xerr=FIG2_ERR, fmt="s", ms=4, color="tab:red",
                capsize=2, label="S1, Fig. 2 (digitized)")
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
    fig2 = [r for r in rows if r["quantity"] == "DT (Fig. 2)"]
    ax.axvspan(-FIG2_ERR, FIG2_ERR, color="tab:red", alpha=0.1, label="Fig. 2 reading band")
    ax.plot([r["diff_mm"] for r in fig2], [r["z_mm"] for r in fig2], "s-", ms=4, color="tab:red",
            label="z3st - Fig. 2")
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
    return hp.finish(fig, TITLE, "bottom clamped at z = 0, core 1500-2500 mm, corners 500/550 C above "
                      "the core; Z rebased to the T_ref of each source (S1 20 C, Table I/II 0 C)", top=0.91)


def figure_tables(rows, tol_mm, tol_rel):
    """Every Table I / II value (and S1 Table 1) against z3st, flagged rows highlighted."""
    tab = [r for r in rows if r["quantity"] != "DT (Fig. 2)"]
    fig = plt.figure(figsize=(17, 1.6 + 0.32 * (len(tab) + 1)))
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.84])
    bc.put_table(ax, tab, highlight=[r["flagged"] for r in tab])
    return hp.finish(fig, TITLE, f"S1 Table 1, Table I / II and IWGFR/75 (../iwgfr.py) against z3st; highlighted: |diff| > {tol_mm} mm "
                      f"and > {100 * tol_rel:.0f} %. Dimensions at mid-wall, z = {Z_BELOW * 1e3:.0f} / "
                      f"{Z_ABOVE * 1e3:.0f} mm, rebased to T_ref = 0 C; cold = undeformed mesh", top=0.88)


def figure_fig2(p):
    """Bow along z in the layout of the benchmark Fig. 2: DT (mm) on x, height
    above the support (COTE, mm) on y, 0-4000 mm every 250 mm."""
    fig, ax = plt.subplots(figsize=(7, 8.5))
    zmm = p["z"] * 1e3
    ax.plot(FIG2[:, 1], FIG2[:, 0], "-s", color="k", ms=4, lw=1.2, label="benchmark AIEA N1, Fig. 2 (S1, digitized)")
    ax.plot(p["w"] * 1e3, zmm, "-", color="tab:blue", lw=2, alpha=0.85, label="z3st")
    ax.plot(p["w_cf"] * 1e3, zmm, "--", color="tab:red", lw=1, label=r"closed form $\kappa=\alpha\,\Delta T/D$")
    bc.source_markers(ax, SOURCES, "Y")
    for zc in (1500, 2500):
        ax.axhline(zc, color="0.7", lw=0.8, ls=":")
    ax.text(12.4, 1500 + 40, "core bottom", fontsize=8, color="0.4", ha="right")
    ax.text(12.4, 2500 + 40, "core top", fontsize=8, color="0.4", ha="right")
    ax.set_xlim(-0.3, 13)
    ax.set_ylim(-50, 4100)
    ax.set_yticks(np.arange(0, 4001, 250))
    ax.set_xticks([0, 5, 10])
    ax.set_xticks(np.arange(0, 13), minor=True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_position(("data", 0))
    ax.set_xlabel("DT (MM)", loc="right")
    ax.set_ylabel("COTE (MM)", rotation=0)
    ax.yaxis.set_label_coords(0.12, 1.0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="center right")
    ax.set_title("BENCHMARK AIEA N1: thermal bow, z3st overlaid", fontstyle="italic")
    fig.tight_layout()
    return fig


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    p = z3st_profiles()
    rows = compare(p)
    bc.write(PREFIX, p, rows, {"": figure(p, rows, bc.TOL_MM, bc.TOL_REL), "_fig2": figure_fig2(p),
                               "_tables": figure_tables(rows, bc.TOL_MM, bc.TOL_REL)}, bc.TOL_MM, bc.TOL_REL)


if __name__ == "__main__":
    main()
