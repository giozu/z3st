#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Catalogue of the stress peaks of the existing 3D runs, to tell physical peaks
from numerical errors. Reads the fields only (no solve).

For every run, height and position of the section it stores the signed
extreme (largest |value|) of each stress component, cell values (DG0):

  local frame of the wall (hex_post.frame; tube: radial n, hoop t)
      nn, ss, zz, sz, ns, nz      sigma_ab = a . sigma . b, a, b in (n, t = s, e_z)
  global frame                    xx, yy, xy
  von Mises                       vm

Positions: 'vertex' (cells within one cell width of a corner, along the
perimeter), 'mid-flat' (within one cell width of the mid-flat), 'section'
(all cells); for the tube the same polar angles as the hexagon of VP1.
Wall layers: inner, mid (only for an odd number of cells), outer; one cell
through the wall gives the single layer 'mid-wall'.

The same peak on different meshes or orders (Q1, Q2, fine), the tube at the
same polar angle and the closed forms classify it (thesis log, sec:peaks):
physical regular, physical singular, numerical, not converged.

Usage (z3st conda env)::

    python peaks_catalogue.py

Writes ``results/peaks/summary.json`` and ``results/peaks/peaks_along_angle.png``.
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hex_post as hp  # noqa: E402

IAEA = os.path.join(HERE, "bowing", "corner-corner", "IAEA")
VP1, VP2 = os.path.join(IAEA, "VP1"), os.path.join(IAEA, "VP2")
OUT = os.path.join(hp.RESULTS, "peaks")
TITLE = "U_hex_shell: stress peaks, physical or numerical"

VP1_LEVELS = {"core entry": 1.5, "mid core": 2.0, "core exit": 2.5, "above core": 3.0, "tip": 4.0}
# label: (case dir, fields, tube?, levels)
RUNS = {
    "VP1 Q1": (VP1, "output/vp1_const_fields_3d.xdmf", False, VP1_LEVELS),
    "VP1 Q2": (os.path.join(VP1, "order2"), "output/vp1_const_Q2_fields_3d.xdmf", False, VP1_LEVELS),
    "VP1 fine (nt 3)": (os.path.join(VP1, "cylinder", "fine_hex"), "output/vp1_hex_fine_fields_3d.xdmf", False, VP1_LEVELS),
    "VP1 tube": (os.path.join(VP1, "cylinder"), "output/vp1_tube_fields_3d.xdmf", True, VP1_LEVELS),
    "VP1 tube fine": (os.path.join(VP1, "cylinder", "fine"), "output/vp1_tube_fine_fields_3d.xdmf", True, VP1_LEVELS),
    "VP1 z-uniform Q1": (os.path.join(VP1, "z_uniform"), "output/vp1_step1_fields_3d.xdmf", False,
                         {"clamp": 0.0, "mid height": 2.0}),
    "VP2 Q1": (VP2, "output/vp2_fields_3d.xdmf", False, {"core entry": 1.5, "core exit": 2.5, "above core": 3.0}),
    "VP2 Q2": (os.path.join(VP2, "order2"), "output/vp2_Q2_fields_3d.xdmf", False,
               {"core entry": 1.5, "core exit": 2.5, "above core": 3.0}),
    "inner_q_outer_T": (os.path.join(HERE, "inner_q_outer_T"), "output/fields.xdmf", False, None),
}
COMPS = ("nn", "ss", "zz", "sz", "ns", "nz", "xx", "yy", "xy", "vm")


def geometry(d, tube, ref=None):
    """Per cell: polar angle, z, local frame (n, t), wall position xi, and the
    perimeter classes; ref is the hexagon frame used for the tube's angles."""
    C = d["X"][d["cells"]].mean(axis=1)
    th = np.arctan2(C[:, 1], C[:, 0])
    if not tube:
        fr = hp.frame(d)
        c = fr["cells"]
        L = fr["L_mid"]
        q = c["p"] % L
        dv, dm = np.minimum(q, L - q), np.abs(q - L / 2)
        return {"theta": th, "z": C[:, 2], "n": c["n"], "t": c["t"], "xi": c["xi"], "dv": dv, "dm": dm, "frame": fr,
                "perimeter": 6 * L}
    r = np.hypot(C[:, 0], C[:, 1])
    n = np.stack([np.cos(th), np.sin(th), 0 * th], axis=1)
    t = np.stack([-np.sin(th), np.cos(th), 0 * th], axis=1)
    xi = r - 0.5 * (r.min() + r.max())
    # arc distance at mid-wall to the vertex and mid-flat angles of the reference hexagon
    a0, R = ref["a0"], 0.5 * (r.min() + r.max())
    def arc(angles):
        dth = np.abs((th[:, None] - angles[None, :] + np.pi) % (2 * np.pi) - np.pi)
        return R * dth.min(axis=1)
    k = np.arange(6)
    return {"theta": th, "z": C[:, 2], "n": n, "t": t, "xi": xi,
            "dv": arc(a0 + k * np.pi / 3), "dm": arc(a0 + k * np.pi / 3 + np.pi / 6), "frame": None,
            "perimeter": 2 * np.pi * R}


def components(S, g):
    ez = np.array([0.0, 0.0, 1.0])
    proj = lambda a, b: np.einsum("ci,cij,cj->c", a, S, b)   # noqa: E731
    Z = np.broadcast_to(ez, g["n"].shape)
    return {"nn": proj(g["n"], g["n"]), "ss": proj(g["t"], g["t"]), "zz": S[:, 2, 2],
            "sz": proj(g["t"], Z), "ns": proj(g["n"], g["t"]), "nz": proj(g["n"], Z),
            "xx": S[:, 0, 0], "yy": S[:, 1, 1], "xy": S[:, 0, 1], "vm": hp.von_mises(S)}


def layers(xi):
    u = np.unique(np.round(xi, 9))
    if len(u) == 1:
        return {"mid-wall": u[0]}
    out = {"inner": u[0], "outer": u[-1]}
    if len(u) % 2 == 1:
        out["mid"] = u[len(u) // 2]
    return out


def signed_extreme(v):
    return float(v[np.argmax(np.abs(v))]) * 1e-6 if len(v) else None


def catalogue(g, comp, levels):
    """{level: {position: {layer: {comp: MPa}}}}, with the cell width along the perimeter."""
    lay = layers(g["xi"])
    out = {}
    for tag, z in levels.items():
        at_z = hp.nearest(g["z"], z)
        first = at_z & (np.abs(g["xi"] - next(iter(lay.values()))) < 1e-9)
        n_ring = len(np.unique(np.round(g["theta"][first], 9)))      # cells around one layer of the ring
        h = g["perimeter"] / max(n_ring, 1)                           # cell width along the perimeter
        pos = {"vertex": g["dv"] < h, "mid-flat": g["dm"] < h, "section": np.ones_like(at_z)}
        out[tag] = {"z_m": float(np.mean(g["z"][at_z])), "cell_width_mm": float(h * 1e3), "positions": {}}
        for pn, pm in pos.items():
            out[tag]["positions"][pn] = {}
            for ln, xv in lay.items():
                m = at_z & pm & (np.abs(g["xi"] - xv) < 1e-9)
                out[tag]["positions"][pn][ln] = {k: signed_extreme(comp[k][m]) for k in COMPS}
    return out


def ring(g, comp, key, z, layer):
    """Values of one component on a ring (height z, wall layer), sorted by polar angle (deg)."""
    xv = layers(g["xi"])[layer]
    m = hp.nearest(g["z"], z) & (np.abs(g["xi"] - xv) < 1e-9)
    o = np.argsort(g["theta"][m])
    return np.degrees(g["theta"][m][o]), comp[key][m][o] * 1e-6


def vp2_closed_form(fr, z):
    """Plane-section sigma_zz of VP2 on the mid-wall line: -E alpha (T - T_plane), T_plane the ds-weighted plane fit."""
    T_vp2 = hp.import_case("vp2_temperature", VP2).T_vp2
    mat = hp.load_material(os.path.join(VP2, "15_15Ti_local.yaml"))
    p = (np.arange(6 * 400) + 0.5) * fr["L_mid"] / 400
    x, y, _ = hp.perimeter_point(fr, p, 0.0)
    T = T_vp2(np.array([x, y, np.full_like(x, z)]))
    A = np.stack([np.ones_like(x), x, y], axis=1)
    c = np.linalg.lstsq(A, T, rcond=None)[0]
    o = np.argsort(np.arctan2(y, x))
    return np.degrees(np.arctan2(y, x))[o], (-mat["E"] * mat["alpha"] * (T - A @ c) * 1e-6)[o]


def main():
    os.makedirs(OUT, exist_ok=True)
    data, summary = {}, {"components": COMPS, "runs": {}}
    ref = None
    for lab, (d, fields, tube, levels) in RUNS.items():
        path = os.path.join(d, fields)
        if not os.path.isfile(os.path.splitext(path)[0] + ".h5"):
            print(f"[SKIP] {lab}: no {path}")
            continue
        if tube and ref is None:
            print(f"[SKIP] {lab}: the tube takes the vertex angles of a hexagon run, none loaded before it")
            continue
        raw = hp.load_z3st(path, temperature="field", case_dir=d)
        g = geometry(raw, tube, ref)
        if ref is None and not tube:
            ref = g["frame"]
        if levels is None:   # mid height of the case
            levels = {"mid height": float(g["z"].min() + 0.5 * np.ptp(g["z"]))}
        comp = components(raw["S"], g)
        data[lab] = (g, comp)
        summary["runs"][lab] = {"source": path, "cells": int(len(raw["cells"])),
                                "layers": {k: float(v * 1e3) for k, v in layers(g["xi"]).items()},
                                "levels": catalogue(g, comp, levels)}
        print(f"[INFO] {lab}: {len(raw['cells'])} cells, layers {list(layers(g['xi']))}")
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"[INFO] {os.path.join(OUT, 'summary.json')}")

    fig, ax = plt.subplots(2, 2, figsize=(16, 11))
    panels = (
        (ax[0, 0], ("VP1 Q1", "VP1 Q2", "VP1 fine (nt 3)", "VP1 tube", "VP1 tube fine"), "zz", 2.5,
         "VP1 core exit (z = 2.5 m): σ_zz, mid-wall"),
        (ax[0, 1], ("VP1 Q1", "VP1 Q2", "VP1 fine (nt 3)", "VP1 tube", "VP1 tube fine"), "ss", 2.5,
         "VP1 core exit (z = 2.5 m): σ_ss, mid-wall"),
        (ax[1, 0], ("VP2 Q1", "VP2 Q2"), "zz", 3.0, "VP2 above the core (z = 3 m): σ_zz, mid-wall"),
        (ax[1, 1], ("VP1 Q1", "VP1 Q2", "VP1 fine (nt 3)"), "vm", 3.0,
         "VP1 above the core (z = 3 m): σ_vm, stress free in closed form"),
    )
    for a, labs, key, z, title in panels:
        for lab in labs:
            if lab not in data:
                continue
            g, comp = data[lab]
            lay = "mid-wall" if "mid-wall" in layers(g["xi"]) else "mid"
            th, v = ring(g, comp, key, z, lay)
            a.plot(th, v, marker="none", lw=1.4, ls="--" if "tube" in lab else "-", label=lab)
        if key == "zz" and labs[0].startswith("VP2") and "VP2 Q1" in data:
            th, v = vp2_closed_form(data["VP2 Q1"][0]["frame"], z)
            a.plot(th, v, color="k", lw=3, alpha=0.3, marker="none", label="plane sections: −Eα (T − T_plane)")
        if ref is not None:
            for k in range(6):
                a.axvline(np.degrees((ref["a0"] + k * np.pi / 3 + np.pi) % (2 * np.pi) - np.pi), color="0.6", lw=0.7, ls=":")
        a.set(title=title, xlabel="polar angle (deg); dotted: vertices", ylabel="MPa")
        if key == "vm":
            a.set_yscale("log")
        a.grid(alpha=0.4)
        a.legend(fontsize=8)
    hp.finish(fig, TITLE, "Cell values (DG0) on the ring nearest to z, middle wall layer; tube at the same polar angle",
               top=0.91)
    path = os.path.join(OUT, "peaks_along_angle.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[INFO] {path}")


if __name__ == "__main__":
    main()
