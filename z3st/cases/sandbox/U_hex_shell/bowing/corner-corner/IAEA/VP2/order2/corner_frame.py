#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP2: hyperstatic ring moment at the corners, closed form against z3st
Q1 (``../``) and Q2 (this folder, ``mechanical.order: 2``).

Above the core T1 is on the walls 1-2 and 6-1, T2 on 3-4 and 4-5, and T is
linear in y on the two walls parallel to y; no wall has a through-wall
gradient. Far from the ends the section is in generalized plane strain: the
wall is a plate with free membrane strain (1+nu) alpha theta - nu eps_zz,
the eps_zz part is linear in (x, y) and compatible, the theta part is not
(T not harmonic: Dundurs does not apply). Hexagonal frame of side R (mid-wall
corner radius), bending stiffness D_b = E t^3 / (12 (1 - nu^2)) per unit
length, inextensible sides: cut on the symmetry axis at corners 1 and 4, the
half frame leaves a horizontal gap (sqrt(3)/2) (1+nu) alpha R (T2 - T1) at
corner 4, closed by H and M (force method):

    M(y) = H y,   M_c = (2 sqrt(3)/5) D_b (1+nu) alpha (T2 - T1) / R   (corners 1, 4)
    sigma_c = 6 M_c / t^2 = (sqrt(3)/5) E alpha (T2 - T1) (t/R) / (1 - nu)

M_c / 2 at corners 2, 3, 5, 6 (IAEA numbering: 1 cold at +y, 4 hot at -y; on the
hex_post axis V2 and V5). The closed form fixes magnitude and shape; the sign
of H y is taken as that of z3st below. ``frame_stiffness`` solves the same
frame by direct stiffness (finite EA) as an independent check.

z3st: plate bending stress on each wall from the nodal u (Q2: vertex values,
exact), w = u.n the mid-wall deflection along the outward normal, k_ss =
d2 w/ds2, k_zz = d2 w/dz2 (bowing of the duct, k_beam n_y on the inclined
walls). Kirchhoff, zeta outward from mid-wall: eps_ss = eps_0 - zeta w,ss, so
the inner face (zeta = -t/2) carries +6 D_b (k_ss + nu k_zz) / t^2 and the
outer face the opposite. Sharp corners in the mesh: the inner corner is singular, so
the stress at the corner itself is mesh-dependent and K_t does not apply;
the frame value is the nominal stress.

Units: E in Pa, M in N m/m, sigma in Pa (plotted in MPa, x 1e-6), k in 1/m,
T in K, y and p in m (plotted in mm, x 1e3).

Usage (z3st conda env)::

    python corner_frame.py

Writes ``vp2_corner_along_s.png``, ``vp2_corner_along_z.png`` and
``vp2_corner_summary.json`` here.
"""

import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
VP2 = os.path.dirname(HERE)
U_HEX = os.path.abspath(os.path.join(VP2, *[os.pardir] * 4))   # U_hex_shell, home of hex_post.py
sys.path.insert(0, U_HEX)
import hex_post as hp  # noqa: E402

# runs of the same case and option: (label, case folder whose vp2_temperature the run used, fields file)
RUNS = (("z3st Q1", VP2, os.path.join(VP2, "output", "vp2_fields_3d.xdmf")),
        ("z3st Q2", HERE, os.path.join(HERE, "output", "vp2_Q2_fields_3d.xdmf")))
Z_PLOT = 3.0         # (m) along-s figure: mid-way between the core top (2.5 m) and the free end (4.0 m)
Z_RANGE = (2.5, 4.0)  # (m) above the core: T1, T2 constant in z
MPA, MM = 1e-6, 1e3   # Pa -> MPa, m -> mm
OUT = os.path.join(HERE, "vp2_corner_")


# --. closed form --..
def frame_closed_form(E, nu, alpha, t, R, dT):
    """H (N/m), M_c (N m/m) at corners 1 and 4, sigma_c (Pa) of the hexagonal frame."""
    Db = E * t**3 / (12 * (1 - nu**2))
    Mc = 2 * np.sqrt(3) / 5 * Db * (1 + nu) * alpha * dT / R
    return {"D_b": Db, "H": Mc / R, "M_c": Mc, "sigma_c": 6 * Mc / t**2}


def frame_stiffness(E, nu, alpha, t, R, theta_of_y):
    """End moments (N m/m) at corners 1..6 of the hexagonal frame by direct
    stiffness: Euler-Bernoulli sides (plate D_b, membrane E t/(1-nu^2)) with
    free axial strain (1+nu) alpha theta, theta averaged along each side."""
    Ep = E / (1 - nu**2)
    EI, EA = Ep * t**3 / 12, Ep * t
    P = R * np.c_[np.cos(np.radians([90, 30, -30, -90, -150, 150])),
                  np.sin(np.radians([90, 30, -30, -90, -150, 150]))]   # corners 1..6
    K, F, sides = np.zeros((18, 18)), np.zeros(18), []
    for i in range(6):
        j = (i + 1) % 6
        d = P[j] - P[i]
        L = np.linalg.norm(d)
        c, s = d / L
        e0 = (1 + nu) * alpha * np.mean(theta_of_y(np.linspace(P[i, 1], P[j, 1], 2001)))
        k = np.zeros((6, 6))
        k[np.ix_([0, 3], [0, 3])] = EA / L * np.array([[1, -1], [-1, 1]])
        b = [1, 2, 4, 5]
        k[np.ix_(b, b)] = EI / L**3 * np.array([[12, 6 * L, -12, 6 * L], [6 * L, 4 * L**2, -6 * L, 2 * L**2],
                                                [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L**2, -6 * L, 4 * L**2]])
        T = np.kron(np.eye(2), np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]]))
        dofs = [3 * i, 3 * i + 1, 3 * i + 2, 3 * j, 3 * j + 1, 3 * j + 2]
        K[np.ix_(dofs, dofs)] += T.T @ k @ T
        F[dofs] += T.T @ np.array([-EA * e0, 0, 0, EA * e0, 0, 0])
        sides.append((dofs, T, k))
    free = [q for q in range(18) if q not in (0, 1, 9)]   # rigid modes: u_x, u_y at corner 1, u_x at corner 4
    u = np.zeros(18)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    return np.array([(k @ (T @ u[dofs]))[2] for dofs, T, k in sides])   # moment at the start corner of each side


# --. z3st --..
def walls(d, fr):
    """Node indices of each wall, inner and outer face, as (z layers, nf + 1) grids sorted by s."""
    X = d["X"]
    a_in, a_out = fr["nodes"]["a"].min(), fr["nodes"]["a"].max()
    out = []
    for m in range(6):
        phi = fr["a0"] + m * np.pi / 3 + np.pi / 6
        n = np.array([np.cos(phi), np.sin(phi), 0.0])
        t = np.array([-np.sin(phi), np.cos(phi), 0.0])
        a = X @ n
        faces = []
        for av in (a_in, a_out):
            idx = np.where(np.abs(a - av) < 1e-7)[0]
            idx = idx[np.lexsort((X[idx] @ t, np.round(X[idx, 2], 9)))]
            nz = np.unique(np.round(X[idx, 2], 9)).size
            faces.append(idx.reshape(nz, -1))
        out.append((n, *faces))
    return out


def plate_bending(d, fr, mat, grids):
    """Per wall: z (layers), y and p (mid-wall nodes), k_ss, k_zz and the
    inner-face bending stress 6 D_b (k_ss + nu k_zz) / t^2 (Pa) at the
    interior nodes of the interior layers."""
    X, u, t = d["X"], d["u"], fr["t"]
    Db = mat["E"] * t**3 / (12 * (1 - mat["nu"] ** 2))
    p_nodes = fr["nodes"]["p"]
    res = []
    for n, gi, go in grids:
        s = 0.5 * (X[gi] @ np.r_[-n[1], n[0], 0] + X[go] @ np.r_[-n[1], n[0], 0])
        w = 0.5 * (u[gi] @ n + u[go] @ n)
        z = X[gi[:, 0], 2]
        ds = np.diff(s, axis=1)
        kss = np.diff(w, 2, axis=1) / (ds[:, :-1] * ds[:, 1:])
        kzz = np.diff(w, 2, axis=0) / (np.diff(z)[:-1] * np.diff(z)[1:])[:, None]
        kss, kzz = kss[1:-1], kzz[:, 1:-1]
        y = 0.5 * (X[gi, 1] + X[go, 1])[1:-1, 1:-1]
        p = 0.5 * (p_nodes[gi] + p_nodes[go])[1:-1, 1:-1]
        res.append({"z": z[1:-1], "y": y, "p": p, "k_ss": kss, "k_zz": kzz,
                    "sigma": 6 * Db * (kss + mat["nu"] * kzz) / t**2})
    return res


def at_z(res, z):
    """Index and value of the node layer nearest to z."""
    k = np.argmin(np.abs(res[0]["z"] - z))
    return k, res[0]["z"][k]


def measures(res, cf):
    """Along z, as FEM / closed form: mean sigma on the inclined walls
    (sign of y) and the slope d sigma / dy fitted on the walls parallel to y."""
    slope_cf = 6 * cf["H"] / cf["t"] ** 2
    incl = [r for r in res if np.all(r["y"] > 0) or np.all(r["y"] < 0)]   # inclined walls (normal 60, 120, 240, 300 deg)
    para = [r for r in res if not (np.all(r["y"] > 0) or np.all(r["y"] < 0))]   # walls parallel to y (normal 0, 180 deg)
    mean = np.mean([np.sign(r["y"].mean()) * r["sigma"].mean(axis=1) / (slope_cf * np.abs(r["y"]).mean(axis=1))
                    for r in incl], axis=0)
    slope = np.mean([(r["sigma"] * r["y"]).sum(1) / (r["y"] ** 2).sum(1) / slope_cf for r in para], axis=0)
    return res[0]["z"], mean, slope


# --. figures --..
def fig_along_s(runs, cf, fr0, hc):
    k, z = at_z(runs[0][1], Z_PLOT)
    fig, ax = plt.subplots(figsize=(16, 5.5))
    p = np.linspace(0, 6 * fr0["L_mid"], 1201)
    _, y, _ = hp.perimeter_point(fr0, p, 0.0)
    ax.plot(p * MM, 6 * cf["H"] * y / cf["t"] ** 2 * MPA, **hp.REF_STYLE["slab"],
            label="closed form (hyperstatic frame), 6 H y / t²")
    for (label, res), color in zip(runs, hp.RUN_COLORS):
        for i, r in enumerate(res):
            o = np.argsort(r["p"][k])
            ax.plot(r["p"][k][o] * MM, r["sigma"][k][o] * MPA, "o-", ms=3, color=color, lw=1.4,
                    label=f"{label}, 6 D_b (κ_ss + ν κ_zz) / t²" if i == 0 else None)
    hp.s_axes(ax, fr0, "ss", flats=True, xlabel=True, hc=hc)
    ax.set_ylabel("plate bending σ_ss at the inner face (MPa)\n(outer face: opposite sign)")
    return hp.finish(fig, "IAEA VP2: hyperstatic ring moment at the corners",
                      f"z = {z * MM:.0f} mm (above the core; the part odd along each wall depends on z, see along_z); M_c = {cf['M_c']:.3f} N m/m, "
                      f"σ_c = 6 M_c / t² = {cf['sigma_c'] * MPA:.2f} MPa at the cold and hot corners (V2, V5), "
                      "half at the others",
                      top=0.86)


def fig_along_z(runs, cf):
    fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
    for (label, res), color in zip(runs, hp.RUN_COLORS):
        z, mean, slope = measures(res, cf)
        axes[0].plot(z, mean, color=color, lw=1.4, label=label)
        axes[1].plot(z, slope, color=color, lw=1.4, label=label)
    for ax, title in zip(axes, ("mean σ on the inclined walls (normal 60°, 120°, 240°, 300°)",
                                "slope dσ/dy on the walls parallel to y (normal 0°, 180°)")):
        ax.axhline(1.0, **hp.REF_STYLE["slab"], label="closed form")
        ax.axvline(Z_PLOT, color="0.5", lw=0.8, ls=":")
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlabel("z (m); dotted: along-s figure")
        ax.grid(alpha=0.35)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("z3st / closed form (-)")
    axes[0].set_ylim(-1.0, 2.5)
    return hp.finish(fig, "IAEA VP2: hyperstatic ring moment along z",
                      "Plane strain away from the ends only: end zones of slow decay at the core top "
                      "(axial kink of T) and at the free end", top=0.86)


def main():
    matplotlib.use("Agg")
    with open(os.path.join(VP2, "post.yaml")) as fh:
        cfg = yaml.safe_load(fh)
    if cfg["option"] != "z-gradient":
        sys.exit(f"option {cfg['option']}: the closed form is written for the VP2 z-gradient case")
    mat = hp.load_material(os.path.join(VP2, cfg["material"]))
    tm = hp.import_case(cfg["temperature"].rsplit(".", 1)[0], VP2)
    dT = tm.T2_TOP - tm.T1_TOP   # (K) hot - cold corner above the core

    runs, d0 = [], None
    for label, case, path in RUNS:
        d = hp.load_z3st(path, temperature=cfg["temperature"], case_dir=case)
        if d0 is not None and (d["X"].shape != d0["X"].shape or not np.allclose(d["X"], d0["X"])
                               or np.abs(d["T"] - d0["T"]).max() > 1e-9):
            sys.exit(f"{label}: mesh or T differs from z3st Q1: not the same case option")
        fr = hp.frame(d)
        if not hp.wall_layers(fr)[1]:
            sys.exit("more than one cell through the wall: mid-wall from the two faces does not apply")
        if d0 is None:
            d0, fr0 = d, fr
        res = plate_bending(d, fr, mat, walls(d, fr))
        keep = (res[0]["z"] > Z_RANGE[0]) & (res[0]["z"] < Z_RANGE[1])
        runs.append((label, [{k: v[keep] for k, v in r.items()} for r in res]))

    R = fr0["a_mid"] / np.cos(np.pi / 6)   # (m) mid-wall corner radius = side of the mid-wall hexagon
    cf = {**frame_closed_form(mat["E"], mat["nu"], mat["alpha"], fr0["t"], R, dT), "t": fr0["t"], "R": R, "dT": dT}
    theta = lambda y: tm.T_vp2(np.array([0 * y, y, 0 * y + Z_PLOT])) - mat["T_ref"]
    M_dsm = frame_stiffness(mat["E"], mat["nu"], mat["alpha"], fr0["t"], R, theta)
    y_c = np.array([1, 0.5, -0.5, -1, -0.5, 0.5])   # corners 1..6: M = H y, y / R
    if not np.allclose(M_dsm * np.sign(M_dsm[0]), cf["M_c"] * y_c, rtol=1e-3):
        sys.exit(f"direct-stiffness frame {np.round(M_dsm, 4)} != closed form M_c y / R, M_c {cf['M_c']:.4f} N m/m")

    figs = [("along_s", fig_along_s(runs, cf, fr0, hp.hot_cold({"frame": fr0, "data": d0}, Z_PLOT))),
            ("along_z", fig_along_z(runs, cf))]
    for name, f in figs:
        f.savefig(f"{OUT}{name}.png", dpi=140)
        plt.close(f)

    summary = {"closed_form": {"R_m": R, "t_m": fr0["t"], "dT_K": dT, "D_b_Nm": cf["D_b"], "H_N_per_m": cf["H"],
                               "M_c_Nm_per_m": cf["M_c"], "sigma_c_MPa": cf["sigma_c"] * MPA},
               "direct_stiffness_M_corners_Nm_per_m": M_dsm.tolist(), "runs": {}}
    for label, res in runs:
        k, z = at_z(res, Z_PLOT)
        z_all, mean, slope = measures(res, cf)
        ana = np.concatenate([6 * cf["H"] * r["y"][k] / fr0["t"] ** 2 for r in res])
        fem = np.concatenate([r["sigma"][k] for r in res])
        summary["runs"][label] = {"z_m": float(z), "max_abs_sigma_MPa": float(np.abs(fem).max() * MPA),
                                  "max_abs_err_vs_closed_form_MPa": float(np.abs(fem - ana).max() * MPA),
                                  "inclined_mean_ratio": float(mean[k]), "parallel_slope_ratio": float(slope[k]),
                                  "inclined_mean_ratio_range_z": [float(mean.min()), float(mean.max())]}
        s = summary["runs"][label]
        print(f"{label}: z = {z:.3f} m  max |sigma| {s['max_abs_sigma_MPa']:.3f} MPa  "
              f"max |err| {s['max_abs_err_vs_closed_form_MPa']:.3f} MPa  inclined mean x{s['inclined_mean_ratio']:.4f}  "
              f"parallel slope x{s['parallel_slope_ratio']:.4f}")
    print(f"closed form: M_c {cf['M_c']:.4f} N m/m (direct stiffness {abs(M_dsm[0]):.4f}), "
          f"sigma_c {cf['sigma_c'] * MPA:.3f} MPa")
    with open(f"{OUT}summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print("written: " + ", ".join(f"{OUT}{n}.png" for n, _ in figs) + f", {OUT}summary.json")


if __name__ == "__main__":
    main()
