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

    python vp1_benchmark.py [--xdmf output/<run>.xdmf] [--tol-mm 0.1] [--tol-rel 0.05]

Outputs in ``U_hex_shell/results/``, named ``IAEA_VP1_vs_benchmark*``: the
overview figure (.png), the bow along z in the layout of Fig. 2 (_fig2.png),
the full comparison with all the sources (_tables.png) and the summary (.json).
A point is flagged when |z3st - benchmark| exceeds both tol-mm and tol-rel of
the benchmark value.
"""

import argparse
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 4))
sys.path.insert(0, U_HEX)
sys.path.insert(0, os.path.dirname(HERE))
import hex_post as hp  # noqa: E402
import iwgfr  # noqa: E402

PREFIX, TITLE = "IAEA_VP1_vs_benchmark", "IAEA VP1: z3st vs IAEA benchmark N1, example 1"
C2K = 273.15

# Fig. 2 (S1), digitized: (z above support (mm), DT (mm))
FIG2 = np.array([
    (1500, 0.00), (1600, 0.00), (1700, 0.00), (1800, 0.02), (1900, 0.05),
    (2000, 0.10), (2100, 0.18), (2200, 0.30), (2300, 0.50), (2400, 0.70),
    (2500, 0.95), (2700, 1.65), (2900, 2.60), (3100, 3.85), (3300, 5.35),
    (3500, 7.10), (3700, 8.95), (3900, 11.20), (4000, 12.35),
])
FIG2_ERR = 0.15   # (mm) reading uncertainty of the scan

STATIONS = {"Core top": 2.5, "LRP": 3.0, "URP": 4.0}   # (m) above the support
Z_BELOW, Z_ABOVE = 1.0, 3.5   # (m) levels of the dimensions: below the core bottom, above the core top
LEVELS = {"hot, below core bottom": Z_BELOW, "hot, above core top": Z_ABOVE, "cold": None}

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


def _note(src, qty, station):
    for key in ((src, qty, station), (src, qty, None), (src, None, None)):
        if key in NOTES:
            return NOTES[key]
    return ""


def section_dims(d, z):
    """Mid-wall distance across the corners and across the flats at height z
    (m): undeformed and deformed. The corner axis is the in-plane axis
    holding the outermost node (y for 'vertical', x for 'flat')."""
    X, u = d["X"], d["u"]
    on = np.abs(X[:, 2] - z) < 1e-6
    if not on.any():
        z = X[np.argmin(np.abs(X[:, 2] - z)), 2]
        on = np.abs(X[:, 2] - z) < 1e-6
    Xs, us = X[on], u[on]
    ic = int(np.argmax(np.abs(Xs[:, :2]).max(axis=0)))   # corner axis
    out = {}
    for name, axis in (("corners", ic), ("flats", 1 - ic)):
        other = 1 - axis
        line = np.abs(Xs[:, other]) < 1e-9 * max(1.0, np.abs(Xs[:, axis]).max())
        ends = []
        for sign in (1, -1):
            sel = line & (sign * Xs[:, axis] > 0)
            if sel.sum() < 2:
                raise RuntimeError(f"no inner/outer node pair on the {name} axis at z = {z}")
            ends.append(((Xs[sel, axis]).mean(), (Xs[sel, axis] + us[sel, axis]).mean()))
        out[name] = (ends[0][0] - ends[1][0], ends[0][1] - ends[1][1])
    return z, out


def z3st_profiles(xdmf):
    """Section-mean bow w (towards the cold corner), u_z and the section
    dimensions from the z3st run, with the closed forms."""
    with open(os.path.join(HERE, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    mat = hp.load_material(os.path.join(HERE, cfg["material"]))
    d = hp.load_z3st(xdmf, temperature=cfg["temperature"], case_dir=HERE)
    cf = hp.closed_form(cfg, mat, HERE, hp.frame(d))
    res = hp.analyse(d, mat, cfg.get("support", "clamped"), cf)
    z = res["z_nodes"]
    dims = {lev: section_dims(d, zl if zl is not None else Z_BELOW) for lev, zl in LEVELS.items()}
    return {"z": z, "w": res["w"], "u_z": res["u_bar"][:, 2],
            "w_cf": cf.w(z), "u_z_cf": cf.u_z(z), "mat": mat, "dims": dims, "dH": iwgfr.outer_flats(d, Z_ABOVE)}


def shift_per_m(p, T_ref_C):
    """Free dilatation (1/1) taking the case T_ref to ``T_ref_C``."""
    return p["mat"]["alpha"] * (p["mat"]["T_ref"] - (T_ref_C + C2K))


def compare(p, tol_mm, tol_rel):
    """One row per benchmark value: source, quantity, station, z, z3st,
    benchmark, diff, rel, flagged, note."""
    rows = []

    def add(src, qty, station, z_m, fem, ref, note=""):
        diff = fem - ref
        rel = diff / ref if abs(ref) > 1e-12 else np.inf
        flag = abs(diff) > tol_mm and abs(rel) > tol_rel
        rows.append({"source": src, "quantity": qty, "station": station, "z_mm": z_m * 1e3,
                     "z3st_mm": fem, "benchmark_mm": ref, "diff_mm": diff, "rel": rel,
                     "flagged": bool(flag), "note": note})

    at = lambda f, z_m: float(np.interp(z_m, p["z"], f) * 1e3)
    for src, data in SOURCES.items():
        sh = shift_per_m(p, data["T_ref_C"])
        for st, z_m in STATIONS.items():
            if st in data["Y"]:
                add(src, "Y", st, z_m, at(p["w"], z_m), data["Y"][st], _note(src, "Y", st))
            if st in data["Z"]:
                add(src, f"Z (T_ref {data['T_ref_C']:.0f} C)", st, z_m,
                    at(p["u_z"], z_m) + sh * z_m * 1e3, data["Z"][st], _note(src, "Z", st))
        for lev, ref in data.get("dH", {}).items():
            z_m, L0, L = p["dH"]
            add(src, f"outer flats increase (T_ref {data['T_ref_C']:.0f} C)", lev, z_m,
                ((L - L0) + sh * L0) * 1e3, ref, _note(src, "dH", lev))
        for lev, (oy, ox) in data["dims"].items():
            z_m, dm = p["dims"][lev]
            hot = LEVELS[lev] is not None
            for name, ref, axis in (("across corners OY", oy, "corners"), ("across flats OX", ox, "flats")):
                L0, L = dm[axis]
                fem = (L + sh * L0) * 1e3 if hot else L0 * 1e3
                add(src, name, lev, z_m if hot else np.nan, fem, ref, _note(src, name, lev))
    for z_mm, dt in FIG2:
        add("S1", "DT (Fig. 2)", "", z_mm / 1e3, at(p["w"], z_mm / 1e3), dt, _note("S1", "DT", None))
    return rows


def _fmt_rel(r):
    return "-" if not np.isfinite(r["rel"]) else f"{100 * r['rel']:+.2f} %"


def _fmt_z(r):
    return "-" if not np.isfinite(r["z_mm"]) else f"{r['z_mm']:.0f}"


def _table_rows(rows):
    return [[r["source"], r["quantity"], r["station"] or "-", _fmt_z(r), f"{r['z3st_mm']:.3f}",
             f"{r['benchmark_mm']:.3f}", f"{r['diff_mm']:+.3f}", _fmt_rel(r), r["note"] or "-"]
            for r in rows]


_HEADER = ["source", "quantity", "station / level", "z (mm)", "z3st (mm)", "benchmark (mm)",
           "diff (mm)", "diff (%)", "note"]
_WIDTHS = [0.13, 0.12, 0.12, 0.05, 0.065, 0.075, 0.065, 0.065, 0.31]


def _put_table(ax, rows, highlight=None):
    tb = hp._table(ax, _HEADER, _table_rows(rows), fill=True)
    for (i, c), cell in tb.get_celld().items():
        cell.set_width(_WIDTHS[c])
        if c == len(_HEADER) - 1 and i > 0:
            cell.set_text_props(ha="left")
            cell.PAD = 0.02
        if highlight is not None and i > 0 and highlight[i - 1]:
            cell.set_facecolor("#ffe0b2")
    tb.set_fontsize(8)
    return tb


def _source_markers(ax, key, p=None, residual=False):
    """Station values of every source as markers: ``key`` 'Y' or 'Z'; one legend entry per label."""
    seen = set()
    for src, data in SOURCES.items():
        st = [k for k in STATIONS if k in data[key]]
        z = np.array([STATIONS[k] for k in st]) * 1e3
        ref = np.array([data[key][k] for k in st])
        if residual:
            fem = np.interp(z / 1e3, p["z"], p["w"]) * 1e3
            x = fem - ref
        else:
            x = ref
        label = f"{data['label']}" + ("" if data["same_case"] else " (other T field)")
        ax.plot(x, z, data["marker"], ms=8, mfc="none", mec=data["color"], mew=1.6,
                label=None if label in seen else label)
        seen.add(label)


def figure(p, rows, tol_mm, tol_rel):
    flagged = [r for r in rows if r["flagged"]]
    n = max(len(flagged), 1)
    fig = plt.figure(figsize=(16, 9 + 0.3 * n))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.2, 0.5 + 0.3 * n], width_ratios=[1.25, 1, 1])
    zmm = p["z"] * 1e3

    # bow
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(p["w"] * 1e3, zmm, "-", color="tab:blue", lw=2, label="z3st")
    ax.plot(p["w_cf"] * 1e3, zmm, "--", color="k", lw=1, label=r"closed form $\kappa=\alpha\,\Delta T/D$")
    ax.errorbar(FIG2[:, 1], FIG2[:, 0], xerr=FIG2_ERR, fmt="s", ms=4, color="tab:red",
                capsize=2, label="S1, Fig. 2 (digitized)")
    _source_markers(ax, "Y")
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
    _source_markers(ax, "Y", p, residual=True)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("z3st - benchmark (mm)")
    ax.set_title("Bow: z3st - benchmark")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower right")

    # axial
    ax = fig.add_subplot(gs[0, 2], sharey=fig.axes[0])
    ax.plot(p["u_z"] * 1e3, zmm, "-", color="tab:blue", lw=2, label="z3st, T_ref 400 C")
    ax.plot(p["u_z_cf"] * 1e3, zmm, "--", color="k", lw=1, label="closed form, T_ref 400 C")
    for T_C, col in ((20.0, "tab:green"), (0.0, "tab:purple")):
        ax.plot(p["u_z"] * 1e3 + shift_per_m(p, T_C) * zmm, zmm, "-", color=col, lw=1.5, alpha=0.7,
                label=f"z3st rebased, T_ref {T_C:.0f} C")
    _source_markers(ax, "Z")
    ax.set_xlabel("axial displacement u_z (mm)")
    ax.set_title("Axial growth (Z)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower right")

    # flagged points
    ax = fig.add_subplot(gs[1, :])
    if flagged:
        _put_table(ax, flagged)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no point beyond the tolerances", ha="center", va="center")
    ax.set_title(f"Points with |diff| > {tol_mm} mm and > {100 * tol_rel:.0f} % of the benchmark",
                 fontsize=10)
    return hp._finish(fig, TITLE, "bottom clamped at z = 0, core 1500-2500 mm, corners 500/550 C above "
                      "the core; Z rebased to the T_ref of each source (S1 20 C, Table I/II 0 C)", top=0.91)


def figure_tables(rows, tol_mm, tol_rel):
    """Every Table I / II value (and S1 Table 1) against z3st, flagged rows highlighted."""
    tab = [r for r in rows if r["quantity"] != "DT (Fig. 2)"]
    fig = plt.figure(figsize=(17, 1.6 + 0.32 * (len(tab) + 1)))
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.84])
    _put_table(ax, tab, highlight=[r["flagged"] for r in tab])
    return hp._finish(fig, TITLE, f"S1 Table 1, Table I / II and IWGFR/75 (../iwgfr.py) against z3st; highlighted: |diff| > {tol_mm} mm "
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
    _source_markers(ax, "Y")
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
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    with open(os.path.join(HERE, "post.yaml")) as fh:
        default = os.path.join(HERE, yaml.safe_load(fh)["fields"])
    ap.add_argument("--xdmf", default=default)
    ap.add_argument("--tol-mm", type=float, default=0.1)
    ap.add_argument("--tol-rel", type=float, default=0.05)
    ap.add_argument("--out", default=hp.RESULTS)
    args = ap.parse_args()

    p = z3st_profiles(args.xdmf)
    rows = compare(p, args.tol_mm, args.tol_rel)
    os.makedirs(args.out, exist_ok=True)
    for suffix, fig in (("", figure(p, rows, args.tol_mm, args.tol_rel)),
                        ("_fig2", figure_fig2(p)),
                        ("_tables", figure_tables(rows, args.tol_mm, args.tol_rel))):
        path = os.path.join(args.out, f"{PREFIX}{suffix}.png")
        fig.savefig(path, dpi=130)
        plt.close(fig)
        print(f"[INFO] {path}")

    print(f"{'source':<34}{'quantity':<36}{'station':<24}{'z':>6}{'z3st':>10}{'bench':>10}{'diff':>9}{'%':>10}")
    for r in rows:
        mark = "  <--" if r["flagged"] else ""
        print(f"{r['source']:<34}{r['quantity']:<36}{r['station']:<24}{_fmt_z(r):>6}{r['z3st_mm']:>10.3f}"
              f"{r['benchmark_mm']:>10.3f}{r['diff_mm']:>+9.3f}{_fmt_rel(r):>10}{mark}")
    with open(os.path.join(args.out, f"{PREFIX}_summary.json"), "w") as fh:
        json.dump({"source": args.xdmf, "tol_mm": args.tol_mm, "tol_rel": args.tol_rel,
                   "z_levels_m": {k: v for k, v in LEVELS.items()}, "rows": rows},
                  fh, indent=2, default=lambda v: None if not np.isfinite(v) else float(v))


if __name__ == "__main__":
    main()
