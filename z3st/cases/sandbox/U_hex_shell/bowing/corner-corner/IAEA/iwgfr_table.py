#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 / VP2: z3st against every reference of the CRP (iwgfr.py and the
SOURCES of vp#_benchmark.py), with the differences in % and the IAEA
declared values (theory, AVERAGE(10), spread, range of the accepted codes).

Table 1, bow Y at core top / LRP / URP: one row per reference, the value
and z3st - reference (%). The IAEA rows: AVERAGE(10) (declared at the URP,
computed here from DT1BOW at core top and LRP), the range min - max of the
10 accepted codes (z3st inside or not) and the accuracy range of those codes
against AVERAGE(10), next to the z3st accuracy (IWGFR/75 paras 17-19).

Table 2: axial growth Z (T_ref 20 C), outer across-flats increase above the
core, the S1 Fig. 2 bow profile (VP1, digitized, +-0.15 mm reading
uncertainty) and the second-order solutions, not comparable with a
small-strain run.

Usage (z3st conda env)::

    python iwgfr_table.py [VP1] [VP2]

Outputs in ``U_hex_shell/results/``: ``IAEA_<case>_vs_codes.png``, ``.csv``
(results_organize.py moves them to the folder of the case).
"""

import csv
import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 3))
sys.path.insert(0, U_HEX)
sys.path.insert(0, HERE)
import hex_post as hp  # noqa: E402
import iwgfr  # noqa: E402

CASES = {"VP1": "ex1", "VP2": "ex2"}
HEAD_Y = ["reference", "kind", "core top (mm)", "diff (%)", "LRP (mm)", "diff (%)", "URP (mm)", "diff (%)", "note"]
HEAD_O = ["quantity", "reference", "station", "reference (mm)", "z3st (mm)", "diff (%)", "note"]
FILL = {"IAEA": "#dbe9f6", "out": "#eeeeee"}


def pct(fem, ref):
    return f"{100 * (fem - ref) / ref:+.2f}"


def z3st_values(case_dir):
    bm = hp._import(f"{os.path.basename(case_dir).lower()}_benchmark", case_dir)
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        xdmf = os.path.join(case_dir, yaml.safe_load(fh)["fields"])
    p = bm.z3st_profiles(xdmf)
    zs = np.array([bm.STATIONS[s] for s in iwgfr.STATIONS])
    Y = np.interp(zs, p["z"], p["w"]) * 1e3
    Z20 = (np.interp(zs, p["z"], p["u_z"]) + bm.shift_per_m(p, 20.0) * zs) * 1e3
    _, L0, L = p["dH"]
    dH20 = ((L - L0) + bm.shift_per_m(p, 20.0) * L0) * 1e3
    return bm, p, Y, Z20, dH20, L0


def table_y(ex, bm, Y):
    rows, fill = [["z3st", "this work", *sum(([f"{y:.3f}", "-"] for y in Y), []), "section-mean bow"]], [None]

    def add(name, kind, ref, note="", shade=None):
        cells = []
        for i, st in enumerate(iwgfr.STATIONS):
            r = ref.get(st) if isinstance(ref, dict) else ref[i]
            cells += ["-", "-"] if r is None else [f"{r:.3f}", pct(Y[i], r)]
        rows.append([name, kind, *cells, note])
        fill.append(shade)

    add("linear theory", "IAEA (DT1BOW)", iwgfr.THEORY[ex], "small strain, 20 C geometry", "IAEA")
    b = iwgfr.band(ex)
    add("AVERAGE(10)", "IAEA (AT1BOW)", [b[s]["average"] for s in iwgfr.STATIONS],
        "URP declared, spread " + ", ".join(f"{s} {b[s]['spread']:.2f} %" for s in iwgfr.STATIONS)
        + " (core top, LRP computed here)", "IAEA")
    for name, (avg, spread) in iwgfr.DECLARED_URP[ex].items():
        if name != "AVERAGE(10)":
            add(name, "IAEA (AT1BOW)", {"URP": avg}, f"spread {spread} % (RAINBOW included)", "IAEA")
    cells, acc = [], []
    for i, st in enumerate(iwgfr.STATIONS):
        inside = b[st]["min"] <= Y[i] <= b[st]["max"]
        cells += [f"{b[st]['min']:.3f}-{b[st]['max']:.3f}", "inside" if inside else "OUTSIDE"]
        a = 100 * (iwgfr.accepted(ex)[:, i] - b[st]["average"]) / b[st]["average"]
        acc += [f"{a.min():+.2f} / {a.max():+.2f}", pct(Y[i], b[st]["average"])]
    rows += [["range of the 10 codes", "IAEA (DT1BOW)", *cells, "min - max of the codes in AVERAGE(10)"],
             ["accuracy of the 10 codes", "IAEA (paras 17-19)", *acc,
              "codes (min / max %) vs AVERAGE(10); diff = z3st accuracy"]]
    fill += ["IAEA", "IAEA"]
    for name, c in iwgfr.CODES.items():
        out = name == iwgfr.REJECTED
        add(name, "DT1BOW", c[ex], f"IWGFR/{c['page']}" + ("; out of AVERAGE(10)" if out else ""),
            "out" if out else None)
    for name, c in iwgfr.SENSITIVITY.items():
        add(name, "sensitivity study", c[ex], f"IWGFR/{c['page']}; {c['note']}", None if c["same_case"] else "out")
    for name, s in bm.SOURCES.items():   # the Part II tables already in the benchmark (not from iwgfr.py)
        if name in iwgfr.CODES or name in iwgfr.SENSITIVITY or name.startswith("IWGFR/75") or not s["Y"]:
            continue
        notes = {st: bm._note(name, "Y", st) for st in iwgfr.STATIONS}
        note = (next(iter(notes.values())) if len(set(notes.values())) == 1
                else "; ".join(f"{st}: {n}" for st, n in notes.items() if n))
        add(s["label"], "IWGFR/76 table", s["Y"], note,
            None if s["same_case"] else "out")
    return rows, fill


def table_other(ex, bm, p, Y, Z20, dH20, L0):
    rows, fill = [], []

    def add(qty, ref_name, station, ref, fem, note="", shade=None):
        rows.append([qty, ref_name, station, f"{ref:.3f}", f"{fem:.3f}", pct(fem, ref), note])
        fill.append(shade)

    for i, st in enumerate(iwgfr.STATIONS):
        add("Z (T_ref 20 C)", "linear theory (DT1LEN)", st, iwgfr.AXIAL_THEORY[i], Z20[i],
            "average temperature inlet -> station", "IAEA")
    if ex == "ex1":   # AT1LEN gives the participants of example 1
        for name, (avg, spread) in iwgfr.AXIAL_DECLARED.items():
            add("Z (T_ref 20 C)", f"{name} (AT1LEN)", "URP", avg, Z20[2], f"spread {spread} %, STRAW (T_ref 0 C) out",
                "IAEA")
        for name, z in iwgfr.AXIAL_URP.items():
            add("Z (T_ref 20 C)", name, "URP", z, Z20[2], "AT1LEN")
    alpha = p["mat"]["alpha"]
    add("outer flats increase (T_ref 20 C)", "linear theory (computed here)", "z = 3500 mm",
        alpha * (525.0 - 20.0) * L0 * 1e3, dH20, "alpha (525 - 20) * outer across flats", "IAEA")
    for name, v in iwgfr.LATERAL[ex].items():
        add("outer flats increase (T_ref 20 C)", name, "z = 3500 mm", v, dH20, f"IWGFR/{iwgfr.CODES[name]['page']}")
    if hasattr(bm, "FIG2"):
        fem = np.interp(bm.FIG2[:, 0] / 1e3, p["z"], p["w"]) * 1e3
        dev = fem - bm.FIG2[:, 1]
        k = np.argmax(np.abs(dev))
        rows.append(["DT profile (S1 Fig. 2)", "S1, digitized", f"{len(dev)} points 0-4000 mm",
                     f"+-{bm.FIG2_ERR:.2f}", f"max |diff| {abs(dev[k]):.3f}", "-",
                     f"{np.sum(np.abs(dev) <= bm.FIG2_ERR)}/{len(dev)} within the reading uncertainty; "
                     f"max at z = {bm.FIG2[k, 0]:.0f} mm"])
        fill.append(None)
    for name, v in iwgfr.SECOND_ORDER_URP[ex].items():
        add("Y", f"linear theory at temperature, {name}", "URP", v, Y[2],
            "second order (DT1BOW): not comparable with a small-strain run", "out")
    return rows, fill


def _draw(ax, header, rows, fill, widths):
    tb = hp._table(ax, header, rows, fill=True)
    for (i, c), cell in tb.get_celld().items():
        cell.set_width(widths[c])
        if c == len(header) - 1 and i > 0:
            cell.set_text_props(ha="left")
            cell.PAD = 0.02
        if i > 0 and fill[i - 1]:
            cell.set_facecolor(FILL[fill[i - 1]])
        if i == 1 and header is HEAD_Y:
            cell.set_text_props(fontweight="bold")
    tb.set_fontsize(8)


def run(case):
    case_dir = os.path.join(HERE, case)
    ex = CASES[case]
    bm, p, Y, Z20, dH20, L0 = z3st_values(case_dir)
    ty, fy = table_y(ex, bm, Y)
    to, fo = table_other(ex, bm, p, Y, Z20, dH20, L0)
    ny, no = len(ty) + 1, len(to) + 1
    H = 1.6 + 0.3 * (ny + no)
    fig = plt.figure(figsize=(20, H))
    gs = fig.add_gridspec(2, 1, height_ratios=[ny, no], hspace=0.08, top=1 - 1.1 / H, bottom=0.01,
                          left=0.01, right=0.99)
    _draw(fig.add_subplot(gs[0]), HEAD_Y, ty, fy, [0.15, 0.10, 0.065, 0.05, 0.065, 0.05, 0.065, 0.05, 0.405])
    _draw(fig.add_subplot(gs[1]), HEAD_O, to, fo, [0.13, 0.20, 0.10, 0.07, 0.09, 0.05, 0.36])
    sub = (f"example {ex[-1]}: diff = (z3st - reference) / reference; blue = IAEA theory / declared values, "
           "grey = out of the IAEA average or not comparable (other T field, second order)")
    fig.suptitle(f"IAEA {case}: z3st vs the CRP codes and the IAEA values", fontsize=15, fontweight="bold",
                 y=1 - 0.15 / H)
    fig.text(0.5, 1 - 0.55 / H, sub, ha="center", va="top", fontsize=10)
    out = os.path.join(hp.RESULTS, f"IAEA_{case}_vs_codes")
    fig.savefig(out + ".png", dpi=130)
    plt.close(fig)
    with open(out + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["table"] + HEAD_Y)
        w.writerows(["Y"] + r for r in ty)
        w.writerow(["table"] + HEAD_O)
        w.writerows(["other"] + r for r in to)
    print(f"[INFO] {out}.png, .csv")
    for r in ty[:6]:
        print("   ", " | ".join(r[:8]))


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    os.makedirs(hp.RESULTS, exist_ok=True)
    for case in sys.argv[1:] or list(CASES):
        run(case)


if __name__ == "__main__":
    main()
