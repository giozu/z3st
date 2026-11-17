# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Comparison of an IAEA VP1 / VP2 run with the published benchmark values,
shared by VP1/vp1_benchmark.py and VP2/vp2_benchmark.py, which hold the data
(SOURCES, NOTES) and the case figures.

A source is a dict with its T_ref (C), Y and Z per station (mm), the mid-wall
dimensions per level (across corners OY, across flats OX, mm), the outer
across-flats increase dH (mm), a plot style and whether it is the same
problem as the z3st case. The z3st Z and dimensions are rebased to the T_ref
of each source by superposing the free dilatation alpha (T_ref,case -
T_ref,source), valid away from the clamp: u_z + alpha dT z and L + alpha dT L0.
A point is flagged when |z3st - benchmark| exceeds both TOL_MM and TOL_REL of
the benchmark value.
"""

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, *[os.pardir] * 3)))
sys.path.insert(0, HERE)
import hex_post as hp  # noqa: E402
import iwgfr  # noqa: E402

C2K = 273.15
TOL_MM, TOL_REL = 0.1, 0.05   # (mm), (-) flag thresholds
STATIONS = {"Core top": 2.5, "LRP": 3.0, "URP": 4.0}   # (m) above the support
Z_BELOW, Z_ABOVE = 1.0, 3.5   # (m) levels of the dimensions: below the core bottom, above the core top
LEVELS = {"hot, below core bottom": Z_BELOW, "hot, above core top": Z_ABOVE, "cold": None}
HEADER = ["source", "quantity", "station / level", "z (mm)", "z3st (mm)", "benchmark (mm)",
          "diff (mm)", "diff (%)", "note"]
WIDTHS = [0.13, 0.12, 0.12, 0.05, 0.065, 0.075, 0.065, 0.065, 0.31]


def note(notes, src, qty, station):
    """Most specific note of ``notes`` for (source, quantity, station)."""
    for key in ((src, qty, station), (src, qty, None), (src, None, None)):
        if key in notes:
            return notes[key]
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
        line = np.abs(Xs[:, 1 - axis]) < 1e-9 * max(1.0, np.abs(Xs[:, axis]).max())
        ends = []
        for sign in (1, -1):
            sel = line & (sign * Xs[:, axis] > 0)
            if sel.sum() < 2:
                raise RuntimeError(f"no inner/outer node pair on the {name} axis at z = {z}")
            ends.append(((Xs[sel, axis]).mean(), (Xs[sel, axis] + us[sel, axis]).mean()))
        out[name] = (ends[0][0] - ends[1][0], ends[0][1] - ends[1][1])
    return z, out


def profiles(case_dir, xdmf=None):
    """Section-mean bow w (towards the cold corner), u_z and the section
    dimensions of the run of ``case_dir`` (post.yaml), with the closed forms;
    ``C`` is the section factor of the closed form (1 for VP1)."""
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    xdmf = xdmf or os.path.join(case_dir, cfg["fields"])
    mat = hp.load_material(os.path.join(case_dir, cfg["material"]))
    d = hp.load_z3st(xdmf, temperature=cfg["temperature"], case_dir=case_dir)
    cf = hp.closed_form(cfg, mat, case_dir, hp.frame(d))
    res = hp.analyse(d, cfg.get("support", "clamped"), cf)
    z = res["z_nodes"]
    dims = {lev: section_dims(d, zl if zl is not None else Z_BELOW) for lev, zl in LEVELS.items()}
    return {"xdmf": xdmf, "z": z, "w": res["w"], "u_z": res["u_bar"][:, 2],
            "w_cf": cf.w(z), "u_z_cf": cf.u_z(z), "C": getattr(cf, "C", 1.0), "mat": mat, "dims": dims,
            "dH": iwgfr.outer_flats(d, Z_ABOVE)}


def shift_per_m(p, T_ref_C):
    """Free dilatation (1/1) taking the case T_ref to ``T_ref_C`` (C)."""
    return p["mat"]["alpha"] * (p["mat"]["T_ref"] - (T_ref_C + C2K))


def row(src, qty, station, z_m, fem, ref, note_, tol_mm, tol_rel):
    """One comparison row; lengths in mm, z in m."""
    diff = fem - ref
    rel = diff / ref if abs(ref) > 1e-12 else np.inf
    return {"source": src, "quantity": qty, "station": station, "z_mm": z_m * 1e3,
            "z3st_mm": fem, "benchmark_mm": ref, "diff_mm": diff, "rel": rel,
            "flagged": bool(abs(diff) > tol_mm and abs(rel) > tol_rel), "note": note_}


def compare(sources, notes, p, tol_mm, tol_rel):
    """One row per benchmark value of ``sources``: Y and Z at the stations,
    the outer across-flats increase and the mid-wall dimensions."""
    rows = []
    at = lambda f, z_m: float(np.interp(z_m, p["z"], f) * 1e3)
    add = lambda *a: rows.append(row(*a, tol_mm, tol_rel))
    for src, data in sources.items():
        sh = shift_per_m(p, data["T_ref_C"])
        for st, z_m in STATIONS.items():
            if st in data["Y"]:
                add(src, "Y", st, z_m, at(p["w"], z_m), data["Y"][st], note(notes, src, "Y", st))
            if st in data["Z"]:
                add(src, f"Z (T_ref {data['T_ref_C']:.0f} C)", st, z_m,
                    at(p["u_z"], z_m) + sh * z_m * 1e3, data["Z"][st], note(notes, src, "Z", st))
        for lev, ref in data.get("dH", {}).items():
            z_m, L0, L = p["dH"]
            add(src, f"outer flats increase (T_ref {data['T_ref_C']:.0f} C)", lev, z_m,
                ((L - L0) + sh * L0) * 1e3, ref, note(notes, src, "dH", lev))
        for lev, (oy, ox) in data["dims"].items():
            z_m, dm = p["dims"][lev]
            hot = LEVELS[lev] is not None
            for name, ref, axis in (("across corners OY", oy, "corners"), ("across flats OX", ox, "flats")):
                L0, L = dm[axis]
                fem = (L + sh * L0) * 1e3 if hot else L0 * 1e3
                add(src, name, lev, z_m if hot else np.nan, fem, ref, note(notes, src, name, lev))
    return rows


def fmt_rel(r):
    return "-" if not np.isfinite(r["rel"]) else f"{100 * r['rel']:+.2f} %"


def fmt_z(r):
    return "-" if not np.isfinite(r["z_mm"]) else f"{r['z_mm']:.0f}"


def put_table(ax, rows, highlight=None):
    """Comparison rows as a table on ``ax``; ``highlight`` marks rows (flagged)."""
    body = [[r["source"], r["quantity"], r["station"] or "-", fmt_z(r), f"{r['z3st_mm']:.3f}",
             f"{r['benchmark_mm']:.3f}", f"{r['diff_mm']:+.3f}", fmt_rel(r), r["note"] or "-"] for r in rows]
    tb = hp.table(ax, HEADER, body, fill=True)
    for (i, c), cell in tb.get_celld().items():
        cell.set_width(WIDTHS[c])
        if c == len(HEADER) - 1 and i > 0:
            cell.set_text_props(ha="left")
            cell.PAD = 0.02
        if highlight is not None and i > 0 and highlight[i - 1]:
            cell.set_facecolor("#ffe0b2")
    tb.set_fontsize(8)
    return tb


def source_markers(ax, sources, key, p=None, residual=False):
    """Station values of every source as markers: ``key`` 'Y' or 'Z'; with
    ``residual`` z3st - source; one legend entry per label."""
    seen = set()
    for data in sources.values():
        st = [k for k in STATIONS if k in data[key]]
        if not st:
            continue
        z = np.array([STATIONS[k] for k in st]) * 1e3
        ref = np.array([data[key][k] for k in st])
        x = np.interp(z / 1e3, p["z"], p["w"]) * 1e3 - ref if residual else ref
        label = data["label"] + ("" if data["same_case"] else " (other T field)")
        ax.plot(x, z, data["marker"], ms=8, mfc="none", mec=data["color"], mew=1.6,
                label=None if label in seen else label)
        seen.add(label)


def flagged_table(ax, rows, tol_mm, tol_rel):
    """Flagged rows on ``ax``, or a note when there is none."""
    flagged = [r for r in rows if r["flagged"]]
    if flagged:
        put_table(ax, flagged)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "no point beyond the tolerances", ha="center", va="center")
    ax.set_title(f"Points with |diff| > {tol_mm} mm and > {100 * tol_rel:.0f} % of the benchmark", fontsize=10)


def axial_panel(ax, p, sources):
    """u_z of z3st and of the closed form at the case T_ref, z3st rebased to 20 C and 0 C, sources' Z."""
    zmm = p["z"] * 1e3
    T_case = p["mat"]["T_ref"] - C2K
    ax.plot(p["u_z"] * 1e3, zmm, "-", color="tab:blue", lw=2, label=f"z3st, T_ref {T_case:.0f} C")
    ax.plot(p["u_z_cf"] * 1e3, zmm, "--", color="k", lw=1, label=f"closed form, T_ref {T_case:.0f} C")
    for T_C, col in ((20.0, "tab:green"), (0.0, "tab:purple")):
        ax.plot(p["u_z"] * 1e3 + shift_per_m(p, T_C) * zmm, zmm, "-", color=col, lw=1.5, alpha=0.7,
                label=f"z3st rebased, T_ref {T_C:.0f} C")
    source_markers(ax, sources, "Z")
    ax.set_xlabel("axial displacement u_z (mm)")
    ax.set_title("Axial growth (Z)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower right")


def write(prefix, p, rows, figures, tol_mm, tol_rel):
    """Figures {suffix: fig} as results/<prefix><suffix>.png, the rows on
    stdout and results/<prefix>_summary.json."""
    os.makedirs(hp.RESULTS, exist_ok=True)
    for suffix, fig in figures.items():
        path = os.path.join(hp.RESULTS, f"{prefix}{suffix}.png")
        fig.savefig(path, dpi=130)
        plt.close(fig)
        print(f"[INFO] {path}")
    print(f"{'source':<34}{'quantity':<36}{'station':<24}{'z':>6}{'z3st':>10}{'bench':>10}{'diff':>9}{'%':>10}")
    for r in rows:
        mark = "  <--" if r["flagged"] else ""
        print(f"{r['source']:<34}{r['quantity']:<36}{r['station']:<24}{fmt_z(r):>6}{r['z3st_mm']:>10.3f}"
              f"{r['benchmark_mm']:>10.3f}{r['diff_mm']:>+9.3f}{fmt_rel(r):>10}{mark}")
    summary = {"source": p["xdmf"], "tol_mm": tol_mm, "tol_rel": tol_rel}
    if p["C"] != 1.0:
        summary["C"] = p["C"]
    summary.update({"z_levels_m": dict(LEVELS), "rows": rows})
    with open(os.path.join(hp.RESULTS, f"{prefix}_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=lambda v: None if not np.isfinite(v) else float(v))
