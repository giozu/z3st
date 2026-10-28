#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Groups the figures of ``results/`` by problem, one subfolder each.

Usage (after ./Allpost)::

    python results_organize.py

For every problem in PROBLEMS:

  - moves ``results/<prefix>*`` into ``results/<subfolder>/`` (names kept);
  - writes ``geometry_2d.png`` (section at the along-s heights: vertices
    V1..V6, flats, s and xi, temperature, hot/cold flat and vertex) and
    ``geometry_3d.png`` (duct, vertex edges, along-s planes);
  - appends to every along-s figure a band with the hot/cold flats and
    vertices at its heights (post.yaml along_s_z, else hex_post.s_heights)
    and the panels whose scale is at the statistical zero (flag:
    convergence/statistics); a new along-s series replaces the older one;
  - writes ``README.md`` with the figure list and the same information.

Statistical-zero flag, per along-s panel (plotted units: MPa, -, mm):
max |v| < ZERO, or (max v - min v) / max |v| < ZERO (constant along s up to
the solver tolerance). Bowing and curvature figures are flagged on max |w|
(mm) and max |kappa| (1/m). An along-s figure already moved is not annotated
twice: only the files found at the top level of ``results/`` are processed.
"""

import glob
import os
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hex_post as hp  # noqa: E402

RESULTS = hp.RESULTS
ZERO = 1e-4

# (subfolder of results/, case folder, {prefix in results/: prefix in the subfolder}, one-line description).
# The along-s heights are post.yaml along_s_z when given (along_s_<group>.png), else hex_post.s_heights.
PROBLEMS = (
    ("inner_q_flux", "inner_q_outer_T", {"inner_q_outer_T_": "inner_q_outer_T_"},
     "Radial through-wall gradient: inner heat flux q, outer T_o fixed. D6 symmetric: no bowing."),
    ("bowing_flat-flat_T_asymmetries", "bowing/flat-flat/T_asymmetries", {"T_asymmetries_": "T_asymmetries_"},
     "Flat-to-flat bowing: Dirichlet T on the nord (cold) and sud (hot) flats, free duct."),
    # IAEA VP1, one folder per option, never compared with each other. hex_post names the figures
    # after the case folder (VP1_, z_uniform_); the option goes in every file name.
    ("bowing_corner-corner_IAEA_VP1/z-gradient", "bowing/corner-corner/IAEA/VP1",
     {"VP1_": "IAEA_VP1_z-gradient_", "IAEA_VP1_vs_benchmark": "IAEA_VP1_z-gradient_vs_benchmark",
      "IAEA_VP1_corners": "IAEA_VP1_z-gradient_corners",
      "IAEA_VP1_along_y": "IAEA_VP1_z-gradient_along_y", "IAEA_VP1_vs_codes": "IAEA_VP1_z-gradient_vs_codes",
      "IAEA_VP1_z-gradient_vs_Ansys_": "IAEA_VP1_z-gradient_vs_Ansys_"},
     "Corner-to-corner bowing, IAEA VP1 WITH axial gradient: T = 400 °C below the core, ramp across "
     "the core (1.5-2.5 m), cold/hot corner 500/550 °C above; clamped bottom. Comparison with the "
     "published benchmark and with Ansys when exported (same option)."),
    ("bowing_corner-corner_IAEA_VP1/z-uniform", "bowing/corner-corner/IAEA/VP1/z_uniform",
     {"z_uniform_": "IAEA_VP1_z-uniform_", "IAEA_VP1_z-uniform_vs_Ansys_": "IAEA_VP1_z-uniform_vs_Ansys_"},
     "Corner-to-corner bowing, IAEA VP1 step 1 WITHOUT axial gradient (Z_UNIFORM): the core-top section "
     "(cold/hot corner 500/550 °C) over the whole height; clamped bottom. Comparison with Ansys "
     "(same option)."),
    ("bowing_corner-corner_IAEA_VP2/z-gradient", "bowing/corner-corner/IAEA/VP2",
     {"VP2_": "IAEA_VP2_z-gradient_", "IAEA_VP2_vs_benchmark": "IAEA_VP2_z-gradient_vs_benchmark",
      "IAEA_VP2_corners": "IAEA_VP2_z-gradient_corners", "IAEA_VP2_along_y": "IAEA_VP2_z-gradient_along_y",
      "IAEA_VP2_vs_codes": "IAEA_VP2_z-gradient_vs_codes"},
     "Corner-to-corner bowing, IAEA VP2 (IWGFR/75 example 2) WITH axial gradient: as VP1, section "
     "non-linear (corners 1, 2, 6 at 500 °C, 3, 4, 5 at 550 °C above the core, linear between corners); "
     "clamped bottom. Comparison with the published benchmark and the CRP codes."),
)
GROUPS = dict(hp.ALONG_S_GROUPS)
PANEL_NAMES = {"vm": "σ_vm", "ss": "hoop σ_ss", "eps_ss": "ε_ss", "eps_zz": "ε_zz",
               "eps_nn": "ε_nn", "u_n": "u_n"}


# --. data --..
def fields_h5(case_rel):
    case_dir = os.path.join(HERE, case_rel)
    with open(os.path.join(case_dir, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    return os.path.splitext(os.path.normpath(os.path.join(case_dir, cfg["fields"])))[0] + ".h5", cfg


def load(case_rel):
    case_dir = os.path.join(HERE, case_rel)
    h5, cfg = fields_h5(case_rel)
    xdmf = os.path.splitext(h5)[0] + ".xdmf"
    if not os.path.isfile(h5):
        return None
    d = hp.load_z3st(xdmf, temperature=cfg.get("temperature", "field"), case_dir=case_dir)
    return hp.analyse(d, None, cfg.get("support", "free"), None)   # the material is not used without closed form


def thermal_sides(res, z):
    """Mean T of each flat, T of each vertex, inner and outer face, at z (K)."""
    fr, T = res["frame"], res["data"]["T"]
    nd = fr["nodes"]
    at_z = hp._nearest(nd["z"], z)
    flats = np.array([T[at_z & (nd["sector"] == m)].mean() for m in range(6)])
    P = 6 * fr["L_mid"]
    verts = []
    for m in range(6):
        dist = np.abs((nd["p"] - m * fr["L_mid"] + P / 2) % P - P / 2)
        sel = at_z & (dist < dist[at_z].min() + 1e-3 * fr["L_mid"])
        verts.append(T[sel].mean())
    xi_in, xi_out = nd["xi"][at_z].min(), nd["xi"][at_z].max()
    inner = T[at_z & hp._nearest(np.where(at_z, nd["xi"], np.inf), xi_in)].mean()
    outer = T[at_z & hp._nearest(np.where(at_z, nd["xi"], np.inf), xi_out)].mean()
    return {"flats": flats, "verts": np.array(verts), "inner": inner, "outer": outer,
            "span": np.ptp(T[at_z])}


def same_around(th):
    """Flats within 5 % of the T range at that height: no hot or cold side
    (the nodal mean of a flat moves by ~2 % with the corner nodes)."""
    return np.ptp(th["flats"]) < 0.05 * th["span"]


def hot_cold_text(res, z, th, name=None):
    fr = res["frame"]
    names = fr["flats"]
    head = (f"{name}, " if name else "") + f"z = {z * 1e3:.1f} mm: "
    if th["span"] < 1e-3:
        return head + f"T uniform ({th['flats'].mean() - 273.15:.1f} °C): no hot or cold flat/vertex"
    f, v = th["flats"], th["verts"]
    lines = [head + f"hot flat {names[f.argmax()]} ({f.max() - 273.15:.1f} °C), "
             f"cold flat {names[f.argmin()]} ({f.min() - 273.15:.1f} °C); "
             f"hot vertex V{v.argmax() + 1} ({v.max() - 273.15:.1f} °C), "
             f"cold vertex V{v.argmin() + 1} ({v.min() - 273.15:.1f} °C)"]
    if same_around(th):
        lines[0] = (head + f"every flat at the same T (spread {np.ptp(f):.1e} K, D6 symmetric): "
                    "no hot or cold flat/vertex")
    dT = th["inner"] - th["outer"]
    hot_face = (f"none, ΔT = {abs(dT):.2g} K negligible against the {th['span']:.3g} K range"
                if abs(dT) < 0.01 * th["span"] + 1e-6 else ("inner" if dT > 0 else "outer"))
    lines.append(f"through the wall: inner face {th['inner'] - 273.15:.2f} °C, outer face "
                 f"{th['outer'] - 273.15:.2f} °C (hot face: {hot_face})")
    return "\n".join(lines)


def along_s_flags(res, z, keys=None):
    """Panels (``keys``, default all) along s at z whose scale is at the statistical zero."""
    fr = res["frame"]
    layers, single = hp._layers(fr)
    out = []
    for key, scale, _ in hp.PANELS:
        if keys is not None and key not in keys:
            continue
        f = fr["nodes"] if key == "u_n" else fr["cells"]
        sides = ("inner", "outer") if (key == "u_n" or not single) else ("inner",)
        v = np.concatenate([hp._ring(f, res[key], z, layers[s][1] if key == "u_n" else layers[s][0])[1]
                            for s in sides]) * scale
        amax, span = np.abs(v).max(), np.ptp(v)
        if amax < ZERO:
            out.append(f"{PANEL_NAMES[key]}: max |v| = {amax:.1e} ≈ 0")
        elif span / amax < ZERO:
            out.append(f"{PANEL_NAMES[key]}: constant along s, variation {span:.1e} "
                       f"({span / amax:.1e} relative) — y-axis zoomed on solver noise")
    return out


def beam_flags(res):
    if "w" not in res:
        return {}
    w, k = np.abs(res["w"]).max() * 1e3, np.abs(res["kappa"]).max()
    out = {}
    if w < ZERO:
        out["vtb_bowing"] = f"max |w| = {w:.1e} mm ≈ 0"
    if k < ZERO:
        out["curvature"] = f"max |kappa| = {k:.1e} 1/m ≈ 0"
    return out


# --. geometry figures --..
def _outline(fr, xi):
    p = np.arange(7) * fr["L_mid"]
    x, y, _ = hp.perimeter_point(fr, np.minimum(p, 6 * fr["L_mid"] - 1e-12), xi)
    x[-1], y[-1] = x[0], y[0]
    return x * 1e3, y * 1e3


def fig_geometry_2d(res, zs, title, path):
    fr, T = res["frame"], res["data"]["T"]
    nd = fr["nodes"]
    half_t = fr["t"] / 2
    ncol = min(len(zs), 2)
    nrow = -(-len(zs) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.5 * ncol, 7.8 * nrow), squeeze=False)
    for ax, (tag, z) in zip(axes.flat, zs.items()):
        at_z = hp._nearest(nd["z"], z)
        X = res["data"]["X"][at_z]
        sc = ax.scatter(X[:, 0] * 1e3, X[:, 1] * 1e3, c=T[at_z] - 273.15, cmap="coolwarm", s=14, zorder=3)
        fig.colorbar(sc, ax=ax, shrink=0.7, label="T (°C)")
        for xi in (-half_t, half_t):
            ax.plot(*_outline(fr, xi), color="k", lw=0.8)
        th = thermal_sides(res, z)
        uniform = th["span"] < 1e-3 or same_around(th)
        R = fr["a_mid"] / np.cos(np.pi / 6)
        for m in range(6):
            ang = fr["a0"] + m * np.pi / 3
            vx, vy = R * np.cos(ang) * 1e3, R * np.sin(ang) * 1e3
            mark = ""
            if not uniform:
                mark = " HOT" if m == th["verts"].argmax() else (" COLD" if m == th["verts"].argmin() else "")
            ax.annotate(f"V{m + 1}{mark}\n{th['verts'][m] - 273.15:.1f} °C", (vx, vy),
                        xytext=(1.17 * vx, 1.17 * vy), ha="center", va="center", fontsize=8.5,
                        color="crimson" if "HOT" in mark else ("navy" if "COLD" in mark else "k"),
                        fontweight="bold" if mark else "normal")
            phi = ang + np.pi / 6
            a = fr["a_mid"] * 0.62e3
            fmark = ""
            if not uniform:
                fmark = " (hot)" if m == th["flats"].argmax() else (" (cold)" if m == th["flats"].argmin() else "")
            ax.text(a * np.cos(phi), a * np.sin(phi),
                    f"{fr['flats'][m]}{fmark}\n{th['flats'][m] - 273.15:.1f} °C",
                    ha="center", va="center", fontsize=8,
                    rotation=np.degrees(phi) % 180 - 90, rotation_mode="anchor")
        # s from V1 counter-clockwise, xi outward, on the first flat
        p0 = np.array([0.35, 0.65]) * fr["L_mid"]
        x0, y0, n = hp.perimeter_point(fr, p0, 0.0)
        ax.annotate("", (x0[1] * 1e3, y0[1] * 1e3), (x0[0] * 1e3, y0[0] * 1e3),
                    arrowprops=dict(arrowstyle="->", color="tab:green", lw=2))
        ax.text(x0.mean() * 1e3, y0.mean() * 1e3, "  s", color="tab:green", fontsize=11, fontweight="bold")
        xm, ym = x0.mean() * 1e3, y0.mean() * 1e3
        L = 0.25 * fr["a_mid"] * 1e3
        ax.annotate("", (xm + L * n[0, 0], ym + L * n[0, 1]), (xm, ym),
                    arrowprops=dict(arrowstyle="->", color="tab:purple", lw=2))
        ax.text(xm + 1.15 * L * n[0, 0], ym + 1.15 * L * n[0, 1], "ξ", color="tab:purple", fontsize=11)
        if "e" in res and np.abs(res["w"]).max() > 1e-9:
            e = res["e"] * 0.35 * fr["a_mid"] * 1e3
            ax.annotate("", tuple(e), (0, 0), arrowprops=dict(arrowstyle="-|>", color="0.3", lw=1.5))
            ax.text(*(1.1 * e), "bowing\ndirection e", fontsize=8, color="0.3", ha="center")
        ax.set_aspect("equal")
        lim = 1.35 * R * 1e3
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(alpha=0.3)
        wall = hot_cold_text(res, z, th).split("\n")[-1] if th["span"] >= 1e-3 else "T uniform"
        ax.set_title(f"{tag}: section at z = {z * 1e3:.1f} mm\n{wall}", fontsize=9)
    fig.suptitle(title, fontweight="bold")
    fig.text(0.5, 0.015, "s: perimeter coordinate at mid-wall, from V1 counter-clockwise (p = 0 at V1); "
             "ξ: through the wall, outward, 0 at mid-wall. Flat and vertex T: nodal means at that z.",
             ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.04 / nrow, 1, 1 - 0.04 / nrow))
    fig.savefig(path, dpi=130)
    plt.close(fig)


def fig_geometry_3d(res, zs, title, path):
    fr, d = res["frame"], res["data"]
    nd = fr["nodes"]
    X, T = d["X"], d["T"] - 273.15
    outer = hp._nearest(nd["xi"], nd["xi"].max())
    R = (fr["a_mid"] + fr["t"] / 2) / np.cos(np.pi / 6) * 1e3
    z0, z1 = fr["z0"] * 1e3, (fr["z0"] + fr["H"]) * 1e3
    # outer faces: one strip per flat and z band, coloured by the mean nodal T
    edges = np.linspace(z0, z1, 61)
    band = np.clip(np.digitize(X[:, 2] * 1e3, edges) - 1, 0, 59)
    norm = plt.Normalize(T.min(), T.max())
    quads, cols = [], []
    for m in range(6):
        a, b = fr["a0"] + m * np.pi / 3, fr["a0"] + (m + 1) * np.pi / 3
        for k in range(60):
            sel = outer & (nd["sector"] == m) & (band == k)
            if not sel.any():
                continue
            quads.append([(R * np.cos(a), R * np.sin(a), edges[k]), (R * np.cos(b), R * np.sin(b), edges[k]),
                          (R * np.cos(b), R * np.sin(b), edges[k + 1]), (R * np.cos(a), R * np.sin(a), edges[k + 1])])
            cols.append(plt.cm.coolwarm(norm(T[sel].mean())))
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(9, 10))
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Poly3DCollection(quads, facecolors=cols, edgecolors="none", alpha=0.55))
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="coolwarm"), ax=ax, shrink=0.5, pad=0.1,
                 label="T on the outer face (°C)")
    ax.set_xlim(-1.2 * R, 1.2 * R)
    ax.set_ylim(-1.2 * R, 1.2 * R)
    ax.set_zlim(z0, z1)
    for m in range(6):
        ang = fr["a0"] + m * np.pi / 3
        vx, vy = R * np.cos(ang), R * np.sin(ang)
        ax.plot([vx, vx], [vy, vy], [z0, z1], color="k", lw=1.2)
        ax.text(1.12 * vx, 1.12 * vy, z1, f"V{m + 1}", fontsize=9, fontweight="bold")
        phi = ang + np.pi / 6
        a = fr["a_mid"] * 1.25e3
        ax.text(a * np.cos(phi), a * np.sin(phi), z1 + 0.04 * (z1 - z0), fr["flats"][m], fontsize=7.5,
                ha="center", color="0.25")
    for xi in (-fr["t"] / 2, fr["t"] / 2):
        for zz in (z0, z1):
            ax.plot(*_outline(fr, xi), zz, color="k", lw=0.6)
    for tag, z in zs.items():
        x, y = _outline(fr, fr["t"] / 2)
        ax.plot(x, y, z * 1e3, color="tab:green", lw=2.5)
        ax.text(x[0] * 1.15, y[0] * 1.15, z * 1e3, f"{tag}\nz = {z * 1e3:.1f} mm",
                color="tab:green", fontsize=8)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z (mm)")
    ax.set_box_aspect((1, 1, min(3.0, max(1.0, fr["H"] / fr["D"] / 2))))
    ax.view_init(elev=22, azim=-60)
    ax.set_title(f"{title}\nouter face coloured by T; black edges = vertices V1..V6; "
                 "green = heights of the along-s figures", fontsize=10, pad=40)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --. annotation of the along-s figures --..
def annotate(png, text, flags):
    img = plt.imread(png)
    h, w = img.shape[:2]
    lines = text.split("\n")
    if flags:
        lines += ["FLAG convergence/statistics (scale at the statistical zero, < 1e-4):"]
        lines += [f"   - {f}" for f in flags]
    else:
        lines += ["no panel at the statistical zero (< 1e-4)"]
    band = 34 * (len(lines) + 1)
    dpi = 130
    fig = plt.figure(figsize=(w / dpi, (h + band) / dpi), dpi=dpi)
    ax = fig.add_axes((0, band / (h + band), 1, h / (h + band)))
    ax.imshow(img)
    ax.axis("off")
    for i, line in enumerate(lines):
        flag = line.startswith("FLAG") or line.startswith("   -")
        fig.text(0.03, (band - 30 * (i + 1)) / (h + band), line, fontsize=12,
                 color="crimson" if flag else "k", fontweight="bold" if line.startswith("FLAG") else "normal",
                 va="bottom")
    fig.savefig(png, dpi=dpi)
    plt.close(fig)


# --. main --..
def readme(sub, desc, files, heights, notes, flags, beam):
    out = [f"# {sub}", "", desc, "",
           "Geometry and coordinates: `geometry_2d.png` (section at the along-s heights) and "
           "`geometry_3d.png`. s runs at mid-wall from V1 counter-clockwise; ξ through the wall, outward.", ""]
    for tag, z in heights.items():
        out += [f"## {tag} (z = {z * 1e3:.1f} mm)", "", "```", notes[tag], "```", ""]
        if flags[tag]:
            out += [f"**FLAG convergence/statistics:** {f}  " for f in flags[tag]] + [""]
        else:
            out += ["No along-s panel at the statistical zero (< 1e-4).", ""]
    if beam:
        out += ["## Bowing figures", ""] + [f"**FLAG convergence/statistics** `*{k}.png`: {v}  "
                                            for k, v in beam.items()] + [""]
    out += ["## Files", ""] + [f"- `{f}`" for f in sorted(files)] + [""]
    return "\n".join(out)


def main():
    index = ["# results", "", "One subfolder per problem; each has `geometry_2d.png`, `geometry_3d.png` "
             "and a README with hot/cold sides and the statistical-zero flags.", ""]
    for sub, case_rel, prefixes, desc in PROBLEMS:
        dest = os.path.join(RESULTS, sub)
        os.makedirs(dest, exist_ok=True)
        fresh = {}
        for src in sorted(prefixes, key=len, reverse=True):     # longest first: IAEA_VP1_ before VP1_
            for f in sorted(glob.glob(os.path.join(RESULTS, src + "*"))):
                fresh.setdefault(f, prefixes[src] + os.path.basename(f)[len(src):])
        res = load(case_rel)
        if res is None:
            print(f"[SKIP] {sub}: no output in {case_rel}")
            continue
        _, cfg = fields_h5(case_rel)
        heights = cfg.get("along_s_z") or hp.s_heights(res)
        notes = {t: hot_cold_text(res, z, thermal_sides(res, z), t) for t, z in heights.items()}
        flags = {t: along_s_flags(res, z) for t, z in heights.items()}
        beam = beam_flags(res)
        # new along-s figures replace every older along-s figure of the same series
        for name in set(fresh.values()):
            if "along_s_" in name:
                for old in glob.glob(os.path.join(dest, name.split("along_s_")[0] + "along_s_*")):
                    os.remove(old)
        for f, name in fresh.items():
            stem = name[:-4].split("along_s_")[-1] if "along_s_" in name else None
            if stem in GROUPS:        # one png per pair of panels, at every height
                keys = GROUPS[stem]
                annotate(f, "\n".join(notes.values()),
                         [f"{t}: {x}" for t, z in heights.items() for x in along_s_flags(res, z, keys)])
            elif stem in heights:
                annotate(f, notes[stem], flags[stem])
            else:
                kind = next((k for k in beam if name.endswith(f"_{k}.png")), None)
                if kind is not None:
                    annotate(f, "bowing along z", [beam[kind]])
            shutil.move(f, os.path.join(dest, name))
            print(f"[INFO] {os.path.relpath(os.path.join(dest, name), HERE)}")
        fig_geometry_2d(res, heights, sub, os.path.join(dest, "geometry_2d.png"))
        fig_geometry_3d(res, heights, sub, os.path.join(dest, "geometry_3d.png"))
        files = [f for f in os.listdir(dest) if f != "README.md"]
        with open(os.path.join(dest, "README.md"), "w") as fh:
            fh.write(readme(sub, desc, files, heights, notes, flags, beam))
        index.append(f"- [{sub}]({sub}/README.md) — {desc}")
    with open(os.path.join(RESULTS, "README.md"), "w") as fh:
        fh.write("\n".join(index) + "\n")


if __name__ == "__main__":
    main()
