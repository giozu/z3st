#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 step 1: z3st against Ansys Mechanical, with the outputs of
``U_hex_shell/hex_post.py`` overlaid (curvature, displacement and deflection along z,
von Mises and hoop stress and strains along s, max stresses and their
position).

Ansys side: ``vp1_ansys_nodes.txt`` and ``vp1_ansys_elems.txt`` written by
``export_ansys_vp1.mac`` (FEM/ansys_vp1). Both codes share the mesh; nodes
and elements are paired by coordinates. Ansys gives stress per element (ETABLE
S, average of the element nodal values) and the thermal strain per node: the
total strain is rebuilt as  C(T)^-1 sigma + eps_th(T) I,  with T read back
from the thermal strain; eps_th and E are those of the card, alpha (T - T_ref)
and E constant or the hooks of a T-dependent card (1515Ti_Tdep). z3st stress
and strain are DG0 (cell midpoint).

Usage (z3st conda env)::

    python compare_ansys.py                                  # both default Ansys exports
    python compare_ansys.py --ansys <dir> [--ansys <dir> ...]
    python compare_ansys.py --case <z3st case folder>          # default z_uniform/: same option as Ansys
    python compare_ansys.py --case . --ansys <dir>             # with the axial gradient (option z-gradient)
    python compare_ansys.py --case 1515Ti_Tdependent --ansys <FEM_1515Ti>/Sobolev               # T-dependent card
    python compare_ansys.py --case 1515Ti_Tdependent/DiMarcello --ansys <FEM_1515Ti>/DiMarcello
    (<FEM_1515Ti> = /mnt/c/Users/Mat/Desktop/Newcleo/FEM/ansys_1515Ti)

The z3st case and the Ansys export must be the same option (post.yaml ``option``):
the T rebuilt from the Ansys thermal strain is checked against the z3st T and the
script stops when they differ by more than T_TOL.

Outputs in ``U_hex_shell/results/`` (with the other subfolders), named
``IAEA_VP1_<option>_vs_Ansys_<figure>``: the hex_post figures with one curve set per
code and the closed forms of ``post.yaml``, a table of the largest
differences, and the summary (.json). With a T-dependent card, one set per
expansion law: ``results/steel_Tdependent/IAEA_VP1_<law>_vs_Ansys_<figure>``.
"""

import argparse
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import yaml
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
U_HEX = os.path.abspath(os.path.join(HERE, *[os.pardir] * 4))   # U_hex_shell, home of hex_post.py
sys.path.insert(0, U_HEX)
import hex_post as hp  # noqa: E402

FEM = "/mnt/c/Users/Mat/Desktop/Newcleo/FEM/ansys_vp1"
# Default Ansys exports per option (post.yaml ``option``): (label, nodes file, elements file).
# z-uniform: 1step/, where vp1_ansys_{nodes,elems}.txt are a copy of keyopt2_0/ and the first
# export (results/) is in the *0.txt files. Other options: give the folder with --ansys.
DEFAULT_ANSYS = {"z-uniform": [
    ("Ansys", os.path.join(FEM, "1step", "vp1_ansys_nodes0.txt"), os.path.join(FEM, "1step", "vp1_ansys_elems0.txt")),
    ("Ansys keyopt2_0", os.path.join(FEM, "keyopt2_0", "vp1_ansys_nodes.txt"),
     os.path.join(FEM, "keyopt2_0", "vp1_ansys_elems.txt")),
]}
# Default case: z_uniform/, the option of the default exports
CASE = os.path.join(HERE, "z_uniform")
T_TOL = 1.0   # (K) max |T Ansys - T z3st| on the nodes: above it the two runs are different options
VOIGT = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))   # Ansys XX YY ZZ XY YZ XZ


def pair(ref, other, what):
    """Row of ``other`` matching each point of ``ref`` (same mesh, any numbering)."""
    dist, idx = cKDTree(other).query(ref)
    tol = 1e-6 * np.ptp(ref, axis=0).max()
    if dist.max() > tol or np.unique(idx).size != idx.size:
        sys.exit(f"{what}: Ansys and z3st do not share the mesh "
                 f"(max pairing distance {dist.max():.3e} m, tolerance {tol:.1e} m)")
    return idx


def card_thermal(mat):
    """eps_th(T) and E(T) of the card: alpha and E constant, or the hooks of a
    T-dependent card (E(T), thermal_strain of the eigenstrain module)."""
    if "alpha" in mat:
        return (lambda T: mat["alpha"] * (T - mat["T_ref"])), (lambda T: np.full_like(T, mat["E"]))
    eps = hp._hook(mat["eigenstrain"].rsplit(".", 1)[0] + ".thermal_strain")
    law = mat["expansion_law"] or "Sobolev"
    return (lambda T: eps(T, mat["T_ref"], law)), hp._hook(mat["E"])


def load_ansys(nodes_txt, elems_txt, z, mat):
    """Ansys solution on the z3st mesh, as a hex_post dataset."""
    nodes = np.loadtxt(nodes_txt, ndmin=2)
    elems = np.loadtxt(elems_txt, ndmin=2)
    C = z["X"][z["cells"]].mean(axis=1)
    nodes = nodes[pair(z["X"], nodes[:, 1:4], f"{nodes_txt}: nodes")]
    elems = elems[pair(C, elems[:, 1:4], f"{elems_txt}: elements")]

    S = np.zeros((len(elems), 3, 3))
    for k, (i, j) in enumerate(VOIGT):
        S[:, i, j] = S[:, j, i] = elems[:, 4 + k]
    eps_th, E_of = card_thermal(mat)
    Tg = mat["T_ref"] + np.linspace(-500.0, 1000.0, 150001)   # eps_th increasing: T by inverse interpolation
    T = np.interp(nodes[:, 7], eps_th(Tg), Tg)
    Tc = T[z["cells"]].mean(axis=1)
    tr = np.trace(S, axis1=1, axis2=2)[:, None, None]
    E = ((1 + mat["nu"]) * S - mat["nu"] * tr * np.eye(3)) / E_of(Tc)[:, None, None]
    E += eps_th(Tc)[:, None, None] * np.eye(3)
    return {"X": z["X"], "cells": z["cells"], "u": nodes[:, 4:7], "T": T, "S": S, "E": E}


def diff_rows(runs):
    """Max |Ansys - z3st| of each Ansys run, on the paired nodes and cells."""
    ref = runs[0][1]
    far = ref["frame"]["cells"]["z"] > 0.5
    rows = []

    def row(name, get, scale, unit, fmt="{:.4g}"):
        a = get(ref)
        rows.append([f"{name} ({unit})", fmt.format(np.abs(a).max() * scale)]
                    + [fmt.format(np.abs(get(r) - a).max() * scale) for _, r in runs[1:]])

    for i, c in enumerate("xyz"):
        row(f"u_{c}", lambda r, i=i: r["data"]["u"][:, i], 1e3, "mm")
    row("bowing w", lambda r: r["w"], 1e3, "mm")
    row("curvature κ", lambda r: r["kappa"], 1.0, "1/m")
    row("σ_vm", lambda r: r["vm"], 1e-6, "MPa")
    row("σ_vm, z > 500 mm", lambda r: r["vm"][far], 1e-6, "MPa")
    row("hoop σ_ss", lambda r: r["ss"], 1e-6, "MPa")
    row("ε_ss", lambda r: r["eps_ss"], 1.0, "-")
    row("ε_zz", lambda r: r["eps_zz"], 1.0, "-")
    return rows


def fig_differences(runs, title):
    rows = diff_rows(runs)
    fig, ax = plt.subplots(figsize=(13, 1.3 + 0.36 * (len(rows) + 1)))
    hp._table(ax, ["quantity", "max |z3st|"] + [f"max |{lab} - z3st|" for lab, _ in runs[1:]], rows, fill=True)
    return hp._finish(fig, title, "Largest differences over the paired nodes (u) and cells (σ, ε); "
                      "z3st: DG0 at the cell midpoint, Ansys: element average of nodal values", top=0.86)


def main():
    matplotlib.use("Agg")   # headless when run as a script; imported, the caller keeps its backend
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ansys", action="append",
                    help="folder with vp1_ansys_nodes.txt and vp1_ansys_elems.txt (repeatable)")
    ap.add_argument("--out", help="default: results/, or results/steel_Tdependent/ for a T-dependent card")
    ap.add_argument("--case", default=CASE,
                    help="z3st case folder (post.yaml + output) run with the same option as the Ansys exports")
    args = ap.parse_args()
    case = os.path.abspath(args.case)
    with open(os.path.join(case, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    option = cfg["option"]
    if args.ansys:
        sets = [(f"Ansys {os.path.basename(os.path.normpath(f))}",
                 os.path.join(f, "vp1_ansys_nodes.txt"), os.path.join(f, "vp1_ansys_elems.txt"))
                for f in args.ansys]
    elif option in DEFAULT_ANSYS:
        sets = DEFAULT_ANSYS[option]
    else:
        sys.exit(f"no default Ansys export for option '{option}' ({case}): give --ansys <dir>")
    if not args.ansys and not all(os.path.isfile(f) for _, *files in sets for f in files):
        print(f"SKIP: default Ansys exports not found under {FEM}: give --ansys <dir>")
        return
    xdmf = os.path.join(case, cfg["fields"])
    mat = hp.load_material(os.path.join(case, cfg["material"]))
    if "alpha" in mat:
        prefix, title, out = f"IAEA_VP1_{option}_vs_Ansys", f"IAEA_VP1 ({option}): z3st vs Ansys", hp.RESULTS
    else:   # T-dependent card: results/steel_Tdependent/, one prefix per expansion law
        law = mat["expansion_law"] or "Sobolev"
        prefix, title = f"IAEA_VP1_{law}_vs_Ansys", f"IAEA_VP1 ({option}), 1515Ti_Tdep {law}: z3st vs Ansys"
        out = os.path.join(hp.RESULTS, "steel_Tdependent")
    args.out = args.out or out
    os.makedirs(args.out, exist_ok=True)
    z = hp.load_z3st(xdmf, temperature=cfg["temperature"], case_dir=case)
    cf = hp.closed_form(cfg, mat, case, hp.frame(z))
    runs = [("z3st", hp.analyse(z, mat, "clamped", cf))]
    for label, nodes_txt, elems_txt in sets:
        a = load_ansys(nodes_txt, elems_txt, z, mat)
        dT = np.abs(a["T"] - z["T"]).max()
        if dT > T_TOL:   # never compare different options (z gradient, reference temperature, ...)
            sys.exit(f"{label}: T rebuilt from the Ansys thermal strain differs from the z3st T of "
                     f"{case} by up to {dT:.1f} K (> {T_TOL} K): not the same option "
                     f"('{option}', T_ref {mat['T_ref'] - 273.15:.1f} °C); nothing written")
        runs.append((label, hp.analyse(a, mat, "clamped", cf)))

    summary = hp.save_all(runs, args.out, f"{prefix}_", title, cfg.get("stations_z"), vtb_3d=False,
                          levels=cfg.get("along_s_z"))
    path = os.path.join(args.out, f"{prefix}_differences.png")
    fig = fig_differences(runs, title)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")
    summary["sources"] = {"z3st": xdmf, **{lab: [n, e] for lab, n, e in sets}}
    summary["max_abs_difference"] = diff_rows(runs)
    summary["option"] = option
    with open(os.path.join(args.out, f"{prefix}_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
