#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 (z-gradient): error of the constant 15-15Ti card (../) against the
T-dependent card, Sobolev set (card of this folder, fields in Sobolev/output)
and Di Marcello set (DiMarcello/).
Same mesh, BCs and T in the three runs (checked on the nodes, T_TOL); every
figure has the three runs, never two laws alone.

Error of the constant card against law L, on a profile q (along z, or along s
at one height):

    Δ_L   = q_constant - q_L                              (unit of q)
    max_L = max |Δ_L|,  rel_L = max_L / max |q_L|  (%)

No closed form, benchmark or tolerance on the runs: eta(T) and E(T) make the
problem non-linear in T, and the differences are the result.

Usage (z3st conda env)::

    python compare_const.py

Outputs in ``U_hex_shell/results/bowing_corner-corner_IAEA_VP1/1515Ti_Tdep/vs_const/``,
named ``IAEA_VP1_1515Ti_Tdep_vs_const_<figure>``: bow, curvature, u_z and
largest section stresses along z (runs and Δ_L), the hex_post along-s and
max-stress figures with max_L, rel_L in the band, the table ``errors`` and
the summary (.json).
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

TDEP = {"Sobolev": HERE, "DiMarcello": os.path.join(HERE, "DiMarcello")}   # expansion_law: case folder
OUTPUT = {"Sobolev": "Sobolev"}   # output/ of the set, when output is the link of Allrun
OUT = os.path.join(hp.RESULTS, "bowing_corner-corner_IAEA_VP1", "1515Ti_Tdep", "vs_const")
PREFIX = "IAEA_VP1_1515Ti_Tdep_vs_const_"
TITLE = "IAEA_VP1 (z-gradient): 1515Ti_Tdep Sobolev and Di Marcello vs constant card"
T_TOL = 1e-6   # (K) max |T - T constant| on the nodes: above it the runs are not the same case
ABOVE = (2.6, 3.9)   # (m) above the core, away from the kink and the top
STYLE = {lab: dict(color=c) for lab, c in zip(("Sobolev", "DiMarcello", "constant"), hp.RUN_COLORS)}   # as hex_post


def run(case_dir, cfg):
    """hex_post analysis of one run, without closed form."""
    mat = hp.load_material(os.path.join(case_dir, cfg["material"]))
    d = hp.load_z3st(os.path.join(case_dir, cfg["fields"]), temperature=cfg["temperature"], case_dir=case_dir)
    return mat, hp.analyse(d, mat, cfg["support"])


def load_runs():
    """[(label, res)]: Sobolev, DiMarcello, constant (last); cfg of this folder.
    None when a run has no output."""
    with open(os.path.join(HERE, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    with open(os.path.join(CONST, "input.yaml")) as fh:
        inp = yaml.safe_load(fh)
    cfgs = {}
    for law, d in TDEP.items():
        with open(os.path.join(d, "post.yaml")) as fh:
            c = yaml.safe_load(fh)
        cfgs[law] = (d, dict(c, fields=os.path.join(OUTPUT.get(law, ""), c["fields"])))
    cfgs["constant"] = (CONST, dict(cfg, fields=f"output/{inp['output']['filename']}.xdmf",
                                    material=inp["materials"]["steel"]))
    for lab, (d, c) in cfgs.items():
        h5 = os.path.splitext(os.path.join(d, c["fields"]))[0] + ".h5"
        if not os.path.isfile(h5):
            print(f"[SKIP] {lab}: no {h5} (run the case first)")
            return None
    runs = []
    for law in TDEP:
        mat, res = run(*cfgs[law])
        if mat["expansion_law"] != law:
            sys.exit(f"{cfgs[law][0]}: card of the {mat['expansion_law']} set, expected {law}")
        runs.append((law, res))
    runs.append(("constant", run(*cfgs["constant"])[1]))
    rc = runs[-1][1]
    for lab, r in runs[:-1]:
        dT = np.abs(r["data"]["T"] - rc["data"]["T"]).max()
        if dT > T_TOL:
            sys.exit(f"T of {lab} and constant differ by up to {dT:.2e} K (> {T_TOL} K): not the same case")
        if not (np.array_equal(r["data"]["cells"], rc["data"]["cells"]) and np.allclose(r["data"]["X"], rc["data"]["X"])):
            sys.exit(f"{lab}: mesh differs from the constant run")
        if not np.allclose(r["e"], rc["e"], atol=1e-6):   # w of every run along the same direction
            sys.exit(f"{lab}: bowing direction {r['e']} differs from the constant run {rc['e']}")
    return runs, cfg, {lab: os.path.join(d, os.path.dirname(c["fields"])) for lab, (d, c) in cfgs.items()}


def error(q_c, q_l):
    """max |q_c - q_l| and the same in % of max |q_l|."""
    m = np.abs(q_c - q_l).max()
    return m, 100.0 * m / np.abs(q_l).max()


# along z: (name, unit, scale, z key, q(res)); stresses as the section maximum
def _sec_max(r, v):
    return hp._max_per_z(r["frame"]["cells"]["z"], v)[1]


PROFILES = {
    "vtb_bowing": (("bowing w", "mm", 1e3, "z_nodes", lambda r: r["w"]),),
    "curvature": (("curvature κ", "1/m", 1.0, "z_sec", lambda r: r["kappa"]),),
    "displacement_z": (("u_z", "mm", 1e3, "z_nodes", lambda r: r["u_bar"][:, 2]),),
    "stress_z": (("max σ_vm in the section", "MPa", 1e-6, None, lambda r: _sec_max(r, r["vm"])),
                 ("max |σ_ss| in the section", "MPa", 1e-6, None, lambda r: _sec_max(r, np.abs(r["ss"])))),
}


def _z(r, key):
    return r[key] if key else hp._max_per_z(r["frame"]["cells"]["z"], r["vm"])[0]


def fig_profiles(runs, items):
    """Per item: the three runs along z (left) and constant - law (right), max_L and rel_L in the legend."""
    rc = runs[-1][1]
    fig, axes = plt.subplots(len(items), 2, figsize=(14, 7.5 * len(items)), sharey=True, squeeze=False)
    errs = {}
    for (ax, dax), (name, unit, scale, zk, get) in zip(axes, items):
        for lab, r in runs:
            ax.plot(get(r) * scale, _z(r, zk) * 1e3, lw=1.6, label=f"{lab}", **STYLE[lab])
        for lab, r in runs[:-1]:
            m, rel = error(get(rc) * scale, get(r) * scale)
            errs[f"{name} ({unit})"] = errs.get(f"{name} ({unit})", {}) | {lab: (m, rel)}
            dax.plot((get(rc) - get(r)) * scale, _z(r, zk) * 1e3, lw=1.6, color=STYLE[lab]["color"],
                     label=f"constant − {lab}: max |Δ| = {m:.4g} {unit} ({rel:.2f} % of max |{lab}|)")
        dax.axvline(0, color="0.4", lw=0.8)
        ax.set(xlabel=f"{name} ({unit})", ylabel="Axial location (mm)")
        dax.set(xlabel=f"constant − T-dependent ({unit})", title=f"Error of the constant card: {name}")
        for a in (ax, dax):
            for zc in (1500, 2500):
                a.axhline(zc, color="0.6", lw=0.8, ls=":")
            a.grid(alpha=0.4)
            a.legend(loc="upper left", fontsize=9)
    sub = ("Same mesh, BCs and T; w and u_z: section mean of the nodes; stresses: cell values (DG0); "
           "dotted lines: core bottom and top\nΔ = constant − T-dependent; % of the largest |value| of that "
           "T-dependent run along z")
    return hp._finish(fig, TITLE, sub, top=0.95 if len(items) > 1 else 0.90), errs


def along_s_errors(runs, levels, keys):
    """max_L and rel_L of each along-s panel key at each height, on the curves
    of hex_post.fig_along_s_levels (_ring, every wall layer drawn):
    [(level, key, law, max, rel, unit)]."""
    rc = runs[-1][1]
    lay, single = hp._layers(rc["frame"])
    out = []
    for tag, z in levels.items():
        for key in keys:
            scale, label = hp.PANEL[key]
            unit = label[label.rfind("(") + 1:-1]
            f = rc["frame"]["nodes"] if key == "u_n" else rc["frame"]["cells"]
            xis = [lay[s][1] if key == "u_n" else lay[s][0] for s in hp._sides(key, single)]
            ring = lambda r: np.concatenate([hp._ring(f, r[key], z, xi)[1] for xi in xis]) * scale
            for lab, r in runs[:-1]:
                m, rel = error(ring(rc), ring(r))
                out.append((tag, key, lab, m, rel, unit))
    return out


def station_rows(runs, zs):
    """Values of the three runs, constant - law and its % of the law, at the stations and maxima."""
    rc = runs[-1][1]
    zc = rc["frame"]["cells"]["z"]
    above = (zc > ABOVE[0]) & (zc < ABOVE[1])
    q = [(f"w at z = {z * 1e3:.0f} mm (mm)", lambda r, z=z: np.interp(z, r["z_nodes"], r["w"]) * 1e3) for z in zs]
    q += [("κ at z = 3000 mm (1/m)", lambda r: np.interp(3.0, r["z_sec"], r["kappa"])),
          ("u_z at the top (mm)", lambda r: r["u_bar"][-1, 2] * 1e3),
          ("max σ_vm (MPa)", lambda r: r["vm"].max() * 1e-6),
          (f"max σ_vm, {ABOVE[0] * 1e3:.0f} < z < {ABOVE[1] * 1e3:.0f} mm (MPa)",
           lambda r: r["vm"][above].max() * 1e-6),
          ("max |σ_ss| (MPa)", lambda r: np.abs(r["ss"]).max() * 1e-6)]
    rows = []
    for name, get in q:
        v = {lab: get(r) for lab, r in runs}
        row = [name] + [f"{v[lab]:.4g}" for lab, _ in runs]
        for lab, _ in runs[:-1]:
            row += [f"{v['constant'] - v[lab]:+.4g}", f"{100 * (v['constant'] - v[lab]) / v[lab]:+.2f}"]
        rows.append(row)
    return rows


def fig_errors(runs, rows, prof):
    laws = [lab for lab, _ in runs[:-1]]
    head = ["quantity"] + [lab for lab, _ in runs] + [h for lab in laws for h in (f"const − {lab}", "(%)")]
    body = rows + [[f"max over z, {name}"] + [""] * len(runs)
                   + [x for lab in laws for x in (f"{e[lab][0]:.4g}", f"{e[lab][1]:.2f}")]
                   for name, e in prof.items()]
    fig, ax = plt.subplots(figsize=(17, 1.5 + 0.36 * (len(body) + 1)))
    tb = hp._table(ax, head, body, fill=True)
    for (_, c), cell in tb.get_celld().items():   # first column wide enough for the quantity names
        cell.set_width(0.25 if c == 0 else 0.75 / (len(head) - 1))
    return hp._finish(fig, TITLE, "Error of the constant card, Δ = constant − T-dependent; (%) of the "
                      "T-dependent value (stations), of max |T-dependent| along z (max over z)", top=0.88)


def save(fig, name):
    path = os.path.join(OUT, f"{PREFIX}{name}.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")
    return path


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    loaded = load_runs()
    if loaded is None:
        return
    runs, cfg, sources = loaded
    os.makedirs(OUT, exist_ok=True)
    prof = {}
    for name, items in PROFILES.items():
        fig, errs = fig_profiles(runs, items)
        save(fig, name)
        prof.update(errs)

    levels = cfg["along_s_z"]
    res = runs[-1][1]
    notes = "\n".join(ro.hot_cold_text(res, z, ro.thermal_sides(res, z), tag) for tag, z in levels.items())
    s_err = []
    for group, keys in hp.ALONG_S_GROUPS:
        path = save(hp.fig_along_s_levels(runs, levels, keys, TITLE), f"along_s_{group}")
        e = along_s_errors(runs, levels, keys)
        s_err += e
        lines = [f"{tag}, {hp.PANEL[key][1].split(' (')[0]}: max |constant − {lab}| = {m:.4g} {unit} ({rel:.2f} %)"
                 for tag, key, lab, m, rel, unit in e]
        ro.annotate(path, notes + "\nΔ over the cell layer (nodes for u_n), % of max |T-dependent| there:\n"
                    + "\n".join(lines),
                    [f"{tag}: {x}" for tag, z in levels.items() for x in ro.along_s_flags(res, z, keys)])
    save(hp.fig_max_table(runs, TITLE), "max_stress")

    zs = [2.5, 3.0, 4.0]   # core top, ACLP, duct top
    rows = station_rows(runs, zs)
    save(fig_errors(runs, rows, prof), "errors")
    summary = {"max_stresses": {lab: hp.max_stresses(r) for lab, r in runs},
               "stations": rows,
               "max_error_along_z": {name: {lab: {"max_abs": m, "rel_pct": rel} for lab, (m, rel) in e.items()}
                                     for name, e in prof.items()},
               "max_error_along_s": [{"level": t, "key": k, "law": lab, "max_abs": m, "rel_pct": rel, "unit": u}
                                     for t, k, lab, m, rel, u in s_err],
               "sources": sources}
    with open(os.path.join(OUT, f"{PREFIX}summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
