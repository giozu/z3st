#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1, corner effect: hexagon vs circular tube of the same mid-wall
perimeter, same T field, material, clamp and mesh layout. Compared at the same
polar angle, so the difference is the effect of the corners.

Runs (all z-gradient, bottom clamped):
  hexagon, IWGFR mesh   ../output/vp1_step1_fields_3d       nf 8, nt 1, 400 layers of 10 mm
  tube, IWGFR mesh      output/vp1_tube_fields_3d           same layout
  hexagon, fine         fine_hex/output/vp1_hex_fine_...    nt 3, 1 mm layers within 100 mm of
  tube, fine            fine/output/vp1_tube_fine_...       the core entry and exit (mesh.geo there)

Stresses are averaged over the cells through the wall (membrane part, the
bending part has zero mean) and, on a cell face, over the two layers around it.

Closed form for the tube (thin shell, axisymmetric part, N_z = 0): at a kink of
the axial T slope (core entry +s0, core exit -s0, s0 = section-mean rise / core
length) the free radial expansion R alpha T has a slope jump, and the shell
closes it with a bending layer of length 1/beta:

    sigma_theta = E alpha s0 / (4 beta) * exp(-|u|) (cos u - sin|u|),   u = beta (z - z_k)
    beta^4 = 3 (1 - nu^2) / (R^2 t^2)

membrane; compared with the perimeter mean of the FE hoop stress (the
cos(phi) part of T averages out), cell-averaged in z as DG0.

Outputs in U_hex_shell/results/:
  IAEA_VP1_tube_along_s.{pdf,png}   sigma_ss, sigma_zz along s at the along-s heights of VP1
  IAEA_VP1_tube_along_z.{pdf,png}   bow, max von Mises, corners and perimeter mean along z
  IAEA_VP1_tube_summary.json
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
VP1 = os.path.dirname(HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "..")))
import hex_post as hp  # noqa: E402

RUNS = {   # name: (xdmf, case dir of vp1_temperature.py, shape)
    "hex_c": (os.path.join(VP1, "output", "vp1_const_fields_3d.xdmf"), VP1, "hex"),
    "tube_c": (os.path.join(HERE, "output", "vp1_tube_fields_3d.xdmf"), HERE, "tube"),
    "hex_f": (os.path.join(HERE, "fine_hex", "output", "vp1_hex_fine_fields_3d.xdmf"), os.path.join(HERE, "fine_hex"), "hex"),
    "tube_f": (os.path.join(HERE, "fine", "output", "vp1_tube_fine_fields_3d.xdmf"), os.path.join(HERE, "fine"), "tube"),
}
STYLE = {"hex_f": dict(color="tab:blue", ls="-", lw=1.3, label="hexagon, fine"),
         "tube_f": dict(color="tab:orange", ls="--", lw=1.3, label="tube, fine (same polar angle)"),
         "hex_c": dict(color="0.55", ls="-", lw=0.8, label="hexagon, IWGFR mesh"),
         "tube_c": dict(color="0.55", ls="--", lw=0.8, label="tube, IWGFR mesh")}
HEIGHTS = {"inlet core": 1.5, "mid core": 2.0, "outlet core": 2.5, "tip": 4.0}   # along_s_z of VP1
IWGFR_Y = {2.5: 1.000, 3.0: 3.249, 4.0: 12.245}   # IWGFR/75 Table I (mm)
OUT = hp.RESULTS


def p_of_phi(fr, phi):
    """Perimeter coordinate of the hexagon mid-wall point at polar angle phi."""
    m = np.floor(((phi - fr["a0"]) % (2 * np.pi)) / (np.pi / 3)).astype(int) % 6
    phi_n = fr["a0"] + m * np.pi / 3 + np.pi / 6
    s_frac = np.tan(phi - phi_n) / np.tan(np.pi / 6)
    return (m + 0.5 * (s_frac + 1)) * fr["L_mid"]


def load(xdmf, case_dir, shape, fr_hex=None):
    d = hp.load_z3st(xdmf, temperature="vp1_temperature.T_vp1", case_dir=case_dir)
    Xn = d["X"][d["cells"]]
    Xc = Xn.mean(axis=1)
    d["zlo"], d["zhi"] = Xn[:, :, 2].min(1), Xn[:, :, 2].max(1)
    d["z"] = 0.5 * (d["zlo"] + d["zhi"])
    d["phi"] = np.arctan2(Xc[:, 1], Xc[:, 0])
    if shape == "hex":
        fr = hp.frame(d)
        t_dir, d["p"], d["fr"] = fr["cells"]["t"], fr["cells"]["p"], fr
    else:
        t_dir = np.stack([-np.sin(d["phi"]), np.cos(d["phi"]), 0 * d["phi"]], axis=1)
        d["p"] = p_of_phi(fr_hex, d["phi"])
    S = d["S"]
    d["ss"] = np.einsum("ci,cij,cj->c", t_dir, S, t_dir)
    d["zz"] = S[:, 2, 2]
    d["vm"] = hp.von_mises(S)
    return d


def at_height(d, z, keys=("ss", "zz")):
    """Layers touching z (two on a cell face), averaged at the same p (and through the wall)."""
    sel = (d["zlo"] <= z + 1e-9) & (d["zhi"] >= z - 1e-9)
    pu, inv = np.unique(np.round(d["p"][sel], 6), return_inverse=True)
    cnt = np.bincount(inv)
    return pu, {k: np.bincount(inv, weights=d[k][sel]) / cnt for k in keys}


def along_z(d, key, how="mean", phi0=None):
    """Per layer: mean or max over the section, or the through-wall mean at the polar angle phi0."""
    zu, inv = np.unique(np.round(d["z"], 7), return_inverse=True)
    v = d[key]
    if phi0 is not None:
        dphi = np.abs(np.angle(np.exp(1j * (d["phi"] - phi0))))
        best = np.full(len(zu), np.inf)
        np.minimum.at(best, inv, dphi)
        sel = dphi <= best[inv] + 1e-6
        return zu, np.bincount(inv[sel], weights=v[sel], minlength=len(zu)) / np.bincount(inv[sel], minlength=len(zu))
    if how == "max":
        out = np.full(len(zu), -np.inf)
        np.maximum.at(out, inv, v)
        return zu, out
    return zu, np.bincount(inv, weights=v) / np.bincount(inv)


def tube_closed_form(mat, Rm, t, m, lo, hi):
    """Perimeter-mean hoop stress (membrane), averaged over [lo, hi] in z."""
    E, nu, al = mat["E"], mat["nu"], mat["alpha"]
    beta = (3 * (1 - nu**2)) ** 0.25 / np.sqrt(Rm * t)
    s0 = (0.5 * (m.T1_TOP + m.T2_TOP) - m.T_INLET) / (m.Z_CORE_TOP - m.Z_CORE_BOT)   # K/m
    F = lambda u: np.exp(-np.abs(u)) * np.sin(u)    # primitive of exp(-|u|)(cos u - sin|u|)
    out = np.zeros_like(lo)
    for zk, sk in ((m.Z_CORE_BOT, s0), (m.Z_CORE_TOP, -s0)):
        out += E * al * sk / (4 * beta) * (F(beta * (hi - zk)) - F(beta * (lo - zk))) / (beta * (hi - lo))
    return out, beta, s0


def layers_lohi(d):
    zu, inv = np.unique(np.round(d["z"], 7), return_inverse=True)
    c = np.bincount(inv)
    return zu, np.bincount(inv, weights=d["zlo"]) / c, np.bincount(inv, weights=d["zhi"]) / c


def extent(z, v, z0, frac=0.1, half=0.09):
    """Height range around z0 where |v| > frac * max|v| (inside the fine zone of 100 mm)."""
    w = np.abs(z - z0) < half
    a = np.abs(v[w])
    keep = a > frac * a.max()
    return [float(z[w][keep].min()), float(z[w][keep].max())]


def bow(d, e):
    zu, inv = np.unique(np.round(d["X"][:, 2], 7), return_inverse=True)
    return zu, np.bincount(inv, weights=d["u"] @ e) / np.bincount(inv)


def main():
    mat = hp.load_material(os.path.join(HERE, "15_15Ti_local.yaml"))
    m = hp._import("vp1_temperature", HERE)
    D = {}
    D["hex_c"] = load(*RUNS["hex_c"])
    fr = D["hex_c"]["fr"]
    for k in ("tube_c", "hex_f", "tube_f"):
        D[k] = load(*RUNS[k], fr_hex=fr)
    Rm, t = 3 * fr["L_mid"] / np.pi, fr["t"]
    phi_hot = m.PHI_HOT                     # hot corner V5 (-y)
    phi_cold = m.PHI_HOT + np.pi            # cold corner V2 (+y)
    phi_flat = m.PHI_HOT + np.pi / 6        # mid-flat of a hot flat (normal 300 deg)
    summary = {"Rm_mm": Rm * 1e3, "t_mm": t * 1e3, "heights": {}}

    # --. along s --..
    fig, axs = plt.subplots(len(HEIGHTS), 2, figsize=(7.2, 8.8), sharex=True)
    verts = np.arange(7) * fr["L_mid"] * 1e3
    for row, (name, z) in enumerate(HEIGHTS.items()):
        prof = {k: at_height(D[k], z) for k in ("hex_c", "hex_f", "tube_f")}
        zt, lo, hi = layers_lohi(D["tube_f"])
        sel = (lo <= z + 1e-9) & (hi >= z - 1e-9)
        cf = tube_closed_form(mat, Rm, t, m, lo[sel], hi[sel])[0].mean() if z < 3.9 else 0.0
        summary["heights"][name] = {"z_m": z, "closed_ss_mean_MPa": float(cf / 1e6)}
        for k, (pu, v) in prof.items():
            summary["heights"][name][k] = {
                "ss_MPa": [float(v["ss"].min() / 1e6), float(v["ss"].max() / 1e6)],
                "ss_mean_MPa": float(v["ss"].mean() / 1e6),
                "zz_MPa": [float(v["zz"].min() / 1e6), float(v["zz"].max() / 1e6)]}
        pu, v = at_height(D["tube_f"], z)
        summary["heights"][name]["tube_f"] = {"ss_MPa": [float(v["ss"].min() / 1e6), float(v["ss"].max() / 1e6)],
                                              "ss_mean_MPa": float(v["ss"].mean() / 1e6),
                                              "zz_MPa": [float(v["zz"].min() / 1e6), float(v["zz"].max() / 1e6)]}
        for col, key in enumerate(("ss", "zz")):
            ax = axs[row, col]
            for k in ("hex_c", "hex_f", "tube_f"):
                pu, v = at_height(D[k], z)
                o = np.argsort(pu)
                ax.plot(pu[o] * 1e3, v[key][o] / 1e6, **STYLE[k])
            ax.axhline(cf / 1e6 if key == "ss" else 0.0, color="k", ls=":", lw=1.0,
                       label="tube, closed form (perimeter mean)" if key == "ss" else r"tube, closed form ($N_z = 0$)")
            for vv in verts:
                ax.axvline(vv, color="0.75", lw=0.6)
            ax.grid(alpha=0.3)
            ax.tick_params(labelsize=7)
            ax.set_title(f"{name}, z = {z * 1e3:.0f} mm", fontsize=8)
            ax.set_ylabel(r"$\sigma_{ss}$ (MPa)" if key == "ss" else r"$\sigma_{zz}$ (MPa)", fontsize=8)
    for ax in axs[-1]:
        ax.set_xticks(verts)
        ax.set_xticklabels([f"V{i % 6 + 1}\n{v:.0f}" for i, v in enumerate(verts)], fontsize=7)
        ax.set_xlabel("s along the perimeter (mm)", fontsize=8)
    axs[0, 0].legend(fontsize=6.2, loc="lower right")
    axs[0, 1].legend(fontsize=6.2, loc="lower right")
    fig.suptitle("IAEA VP1: hexagon vs tube of the same perimeter (through-wall mean, V5 hot, V2 cold)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"IAEA_VP1_tube_along_s.{ext}"), dpi=150)
    plt.close(fig)

    # --. along z --..
    e = -np.array([np.cos(phi_hot), np.sin(phi_hot), 0.0])          # towards the cold corner
    cfm = hp.VP1("vp1_temperature", mat, VP1, fr)
    fig = plt.figure(figsize=(7.4, 6.4))
    gs = fig.add_gridspec(2, 3)
    ax = fig.add_subplot(gs[0, 0])
    for k in ("hex_c", "hex_f", "tube_f"):
        zb, wb = bow(D[k], e)
        ax.plot(wb * 1e3, zb, **STYLE[k])
    zb, _ = bow(D["hex_c"], e)
    ax.plot(cfm.w(zb) * 1e3, zb, ":", color="k", lw=1.0, label=r"beam, $\kappa = \alpha\Delta T/D$")
    ax.plot(list(IWGFR_Y.values()), list(IWGFR_Y.keys()), "o", mfc="none", color="tab:red", ms=4, label="IWGFR/75 Table I")
    ax.set_xlabel("bow towards the cold corner (mm)", fontsize=7.5)
    ax.set_ylabel("z (m)", fontsize=7.5)
    ax.set_title("bow", fontsize=8)
    ax.legend(fontsize=5.8, loc="upper left")
    ax = fig.add_subplot(gs[1, 0])
    for k in ("hex_c", "hex_f", "tube_f"):
        zu, v = along_z(D[k], "vm", "max")
        ax.plot(v / 1e6, zu, **STYLE[k])
    ax.set_xlabel(r"max $\sigma_{vm}$ of the section (MPa)", fontsize=7.5)
    ax.set_ylabel("z (m)", fontsize=7.5)
    ax.set_title("magnitude along z", fontsize=8)
    for (r, c, z0, key, title) in ((0, 1, 1.5, "zz", r"$\sigma_{zz}$, core entry"),
                                   (0, 2, 2.5, "zz", r"$\sigma_{zz}$, core exit"),
                                   (1, 1, 1.5, "ss", r"$\sigma_{ss}$, core entry")):
        ax = fig.add_subplot(gs[r, c])
        for phi0, col, lab in ((phi_hot, "tab:red", "hot corner"), (phi_cold, "tab:blue", "cold corner"),
                               (phi_flat, "tab:green", "hot mid-flat")):
            zu, v = along_z(D["hex_f"], key, phi0=phi0)
            ax.plot(v / 1e6, zu, "-", color=col, lw=1.2, label=f"hexagon, {lab}")
            zu, v = along_z(D["tube_f"], key, phi0=phi0)
            ax.plot(v / 1e6, zu, "--", color=col, lw=0.9, label=f"tube, same angle")
        ax.set_ylim(z0 - 0.1, z0 + 0.1)
        ax.set_xlabel(f"{title} (MPa)", fontsize=7.5)
        ax.set_title(f"{title}, fine", fontsize=8)
        if (r, c) == (0, 1):
            ax.legend(fontsize=5.4, loc="upper left")
    ax = fig.add_subplot(gs[1, 2])
    for k in ("hex_f", "tube_f", "hex_c"):
        zu, v = along_z(D[k], "ss", "mean")
        ax.plot(v / 1e6, zu, **STYLE[k])
    zt, lo, hi = layers_lohi(D["tube_f"])
    cfz, beta, s0 = tube_closed_form(mat, Rm, t, m, lo, hi)
    ax.plot(cfz / 1e6, zt, ":", color="k", lw=1.2, label="tube, closed form")
    ax.set_ylim(1.4, 1.6)
    ax.set_xlabel(r"perimeter mean of $\sigma_{ss}$ (MPa)", fontsize=7.5)
    ax.set_title(r"$\sigma_{ss}$, perimeter mean, core entry", fontsize=8)
    ax.legend(fontsize=5.6, loc="upper right")
    for a in fig.axes:
        a.grid(alpha=0.3)
        a.tick_params(labelsize=7)
        for zk in (m.Z_CORE_BOT, m.Z_CORE_TOP):
            a.axhline(zk, color="0.6", lw=0.6, ls="--")
    fig.suptitle("IAEA VP1 along z: hexagon vs tube (same T, perimeter, clamp, mesh); through-wall mean", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"IAEA_VP1_tube_along_z.{ext}"), dpi=150)
    plt.close(fig)

    # --. numbers --..
    def corner(k, key, phi0, z0, how):
        zu, v = along_z(D[k], key, phi0=phi0)
        w = np.abs(zu - z0) < 0.09
        return float((v[w].min() if how == "min" else v[w].max()) / 1e6)
    pm = {k: along_z(D[k], "ss", "mean") for k in D}
    vmx = {k: along_z(D[k], "vm", "max") for k in D}
    above = lambda k, key: float(np.abs(D[k][key][(D[k]["z"] > 2.8) & (D[k]["z"] < 3.8)]).max() / 1e6)
    summary.update({
        "beta_1_per_m": float(beta), "1_over_beta_mm": float(1e3 / beta), "s0_K_per_m": float(s0),
        "closed_peak_MPa": float(mat["E"] * mat["alpha"] * s0 / (4 * beta) / 1e6),
        "bow_URP_mm": {k: float(bow(D[k], e)[1][-1] * 1e3) for k in D} | {"beam": float(cfm.w(np.array([4.0]))[0] * 1e3)},
        "perimeter_mean_ss_entry_MPa": {k: float(pm[k][1][np.abs(pm[k][0] - 1.5) < 0.09].max() / 1e6) for k in D}
                                       | {"closed_fine_cell": float(cfz.max() / 1e6)},
        "extent_10pct_entry_m": {"hex_f": extent(*pm["hex_f"], 1.5), "tube_f": extent(*pm["tube_f"], 1.5),
                                 "closed": extent(zt, cfz, 1.5),
                                 "hex_f_vm": extent(vmx["hex_f"][0], vmx["hex_f"][1], 1.5),
                                 "tube_f_vm": extent(vmx["tube_f"][0], vmx["tube_f"][1], 1.5)},
        "zz_MPa": {f"{k}_{c}_{zname}": corner(k, "zz", ph, z0, how)
                   for k in ("hex_f", "tube_f", "hex_c")
                   for c, ph in (("hot", phi_hot), ("cold", phi_cold), ("flat", phi_flat))
                   for zname, z0, how in (("entry", 1.5, "min" if c != "flat" else "max"),
                                          ("exit", 2.5, "max" if c != "flat" else "min"))},
        "max_vm_MPa": {k: float(vmx[k][1].max() / 1e6) for k in D},
        "floor_above_core_MPa": {k: {"ss": above(k, "ss"), "zz": above(k, "zz"), "vm": above(k, "vm")} for k in D},
    })
    with open(os.path.join(OUT, "IAEA_VP1_tube_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
