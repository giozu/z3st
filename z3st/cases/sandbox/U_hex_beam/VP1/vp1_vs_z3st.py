#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d case 1D-1: IAEA VP1, beam vs Z3ST 3D vs IWGFR/75
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Bow of the 1D-1 beam (output/fields.npz) against the Z3ST 3D run of the same
case and against IWGFR/75 example 1: linear theory and the 10 accepted
participants (../iwgfr.py).

From the Z3ST field, at each node level:
    w    section mean of the nodes towards the cold corner (= vp1_benchmark),
    psi  section rotation, plane fit u_z = c_0 + psi h (u_z = h w' for a plane
         section that stays normal to the axis, Euler-Bernoulli).
In the core the transverse T difference grows with z and the 3D section shears:
gamma = w' - psi = -c dkappa*/dz (beam1d.section.HexSection.thermal_shear), so
the section mean lags the beam by -c (kappa*(z) - kappa*(0)), a constant
-c kappa_0 above the core. The rest of Z3ST - beam is psi - w'_beam.

The Z3ST profile is read once from the xdmf of post.yaml and cached in
output/z3st_VP1_profile.npz; it is read again when the xdmf is newer.

Output: output/IAEA_VP1_beam1d_vs_z3st.png and _summary.json.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402

from beam1d import case as case_mod  # noqa: E402
from beam1d import solver  # noqa: E402

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(CASE_DIR, "output")
PREFIX = "IAEA_VP1_beam1d_vs_z3st"

CASE = case_mod.load_case(CASE_DIR)
CASE_3D = os.path.normpath(os.path.join(CASE_DIR, CASE.inp["case_3d"]))
sys.path.insert(0, CASE_3D)
import vp1_benchmark as vb  # noqa: E402  (puts ../iwgfr.py and hex_post on the path)
import vp1_temperature as vp1  # noqa: E402
import iwgfr  # noqa: E402

STATIONS = {"Core top": CASE.stations["core_top"], "LRP": CASE.stations["LRP"], "URP": CASE.stations["top"]}


def z3st_profile():
    """z, w, psi of the Z3ST run (node levels), cached next to the beam fields."""
    with open(os.path.join(CASE_3D, "post.yaml")) as f:
        cfg = yaml.safe_load(f)
    xdmf = os.path.join(CASE_3D, cfg["fields"])
    cache = os.path.join(OUT_DIR, "z3st_VP1_profile.npz")
    if os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(xdmf):
        c = np.load(cache)
        if "uz" in c:
            return c["z"], c["w"], c["psi"], c["uz"], xdmf
    print(f"[INFO] Reading the Z3ST run: {xdmf}")
    d = vb.hp.load_z3st(xdmf, temperature=cfg["temperature"], case_dir=CASE_3D)
    X, u = d["X"], d["u"]
    e_hot = np.array([np.cos(vp1.PHI_HOT), np.sin(vp1.PHI_HOT)])
    z = np.unique(np.round(X[:, 2], 9))
    w, psi, uz = np.empty_like(z), np.empty_like(z), np.empty_like(z)
    for i, zl in enumerate(z):
        on = np.abs(X[:, 2] - zl) < 1e-6
        h = X[on, :2] @ e_hot
        w[i] = -(u[on, :2] @ e_hot).mean()
        uz[i] = u[on, 2].mean()
        psi[i] = np.linalg.lstsq(np.column_stack([np.ones_like(h), h]), u[on, 2], rcond=None)[0][1]
    np.savez(cache, z=z, w=w, psi=psi, uz=uz)
    return z, w, psi, uz, xdmf


def main():
    d = np.load(os.path.join(OUT_DIR, "fields.npz"))
    # hot corner at -y ('vertical'): the bow towards the cold corner is w = u_y
    z_b, U_w = d["z"], d["U_y"]
    z3, w3, psi3, uz3, xdmf = z3st_profile()
    wb, dwb, _ = solver.evaluate(z_b, U_w, z3)
    dw3 = np.gradient(w3, z3)

    sec = CASE.section
    c = sec.thermal_shear()
    R = sec.a_m
    k3 = np.interp(z3, d["z_s"], d["ky_s"])
    dw_shear = -c * (k3 - k3[0])            # -c (kappa*(z) - kappa*(0))
    gamma_cf = -c * np.gradient(k3, z3)
    core = (z3 > vp1.Z_CORE_BOT + 0.1) & (z3 < vp1.Z_CORE_TOP - 0.1)
    gamma3 = float(np.mean((dw3 - psi3)[core]))
    print(f"[INFO] Thermal shear: c = {c:.4e} m^2 = {c / R**2:.4f} R^2 (R = {R * 1e3:.2f} mm)")
    print(f"[INFO]   gamma in the core: Z3ST {gamma3:.4e} rad, thin wall {float(np.mean(gamma_cf[core])):.4e} rad")

    theory = dict(zip(iwgfr.STATIONS, iwgfr.THEORY["ex1"]))
    band = iwgfr.band("ex1")
    rows = []
    for st, zz in STATIONS.items():
        b = float(solver.evaluate(z_b, U_w, zz)[0][0]) * 1e3
        f3 = float(np.interp(zz, z3, w3)) * 1e3
        sh = float(np.interp(zz, z3, dw_shear)) * 1e3
        uzb = float(solver.evaluate_rod(z_b, d["U_z"], zz)[0][0]) * 1e3
        uzf = float(np.interp(zz, z3, uz3)) * 1e3
        rows.append({"station": st, "z_mm": zz * 1e3, "theory_mm": theory[st], "beam_mm": b, "z3st_mm": f3,
                     "z3st_minus_beam_mm": f3 - b, "z3st_minus_beam_pct": 100.0 * (f3 - b) / b,
                     "thermal_shear_mm": sh, "rest_mm": f3 - b - sh,
                     "uz_beam_mm": uzb, "uz_z3st_mm": uzf, "uz_z3st_minus_beam_mm": uzf - uzb,
                     "avg10_mm": band[st]["average"], "avg10_min_mm": band[st]["min"],
                     "avg10_max_mm": band[st]["max"]})
        print(f"[INFO] {st:<9} z = {zz:.1f} m: theory {theory[st]:.3f}, beam {b:.4f}, Z3ST {f3:.4f} mm, "
              f"Z3ST - beam {f3 - b:+.4f} mm ({100.0 * (f3 - b) / b:+.2f} %): thermal shear {sh:+.4f}, "
              f"rest {f3 - b - sh:+.4f} mm; u_z beam {uzb:.4f}, Z3ST {uzf:.4f} mm ({uzf - uzb:+.4f})")

    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.2, 1.0], hspace=0.3, wspace=0.28)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)

    ax1.plot(wb * 1e3, z3 * 1e3, "-", color="#0072B2", lw=2.0, label="beam 1D-1 (Euler-Bernoulli, Hermite)")
    ax1.plot(w3 * 1e3, z3 * 1e3, "--", color="k", lw=1.2, label="Z3ST 3D (section mean of the nodes)")
    zs_mm = [STATIONS[s] * 1e3 for s in iwgfr.STATIONS]
    ax1.plot(iwgfr.THEORY["ex1"], zs_mm, "s", color="tab:red", mfc="none", ms=9, label="IWGFR/75 linear theory")
    for i, Y in enumerate(iwgfr.accepted("ex1")):
        ax1.plot(Y, zs_mm, "x", color="grey", ms=6, label="IWGFR/75 participants, AVERAGE(10)" if i == 0 else None)
    ax1.set_xlabel("bow w towards the cold corner (mm)")
    ax1.set_ylabel("height above support (mm)")
    ax1.set_title("Thermal bow")

    ax2.plot((w3 - wb) * 1e6, z3 * 1e3, "-", color="k", lw=1.6, label="Z3ST - beam")
    ax2.plot(dw_shear * 1e6, z3 * 1e3, "--", color="#D55E00", lw=1.4,
             label=r"thermal shear $-c\,(\kappa^* - \kappa^*_0)$")
    ax2.plot((w3 - wb - dw_shear) * 1e6, z3 * 1e3, ":", color="#009E73", lw=1.6, label="rest")
    ax2.axvline(0.0, color="#0072B2", lw=1.0)
    ax2.set_xlabel(r"Z3ST - beam ($\mu$m)")
    ax2.set_title("Difference from the beam (= closed form)")

    ax3.plot((dw3 - psi3) * 1e6, z3 * 1e3, "-", color="k", lw=1.6, label=r"Z3ST $\gamma = w' - \psi$")
    ax3.plot(gamma_cf * 1e6, z3 * 1e3, "--", color="#D55E00", lw=1.4,
             label=rf"thin wall $-c\,d\kappa^*/dz$, $c$ = {c / R**2:.3f} $R^2$")
    ax3.plot((psi3 - dwb) * 1e6, z3 * 1e3, ":", color="#009E73", lw=1.6, label=r"Z3ST $\psi$ - beam $w'$")
    g_min = float(gamma_cf.min()) * 1e6
    ax3.set_xlim(1.8 * g_min, -0.6 * g_min)   # the one-sided w' of the last levels (tip) is off scale
    ax3.set_xlabel(r"rotation ($\mu$rad)")
    ax3.set_title(r"Shear angle and section rotation (tip off scale)")

    for ax in (ax1, ax2, ax3):
        for zz in (1500.0, 2500.0):
            ax.axhline(zz, color="grey", ls=":", lw=0.8)
        ax.grid(True, alpha=0.4)
        ax.legend(loc="lower right" if ax is ax1 else "upper left", fontsize=8.5)
    ax1.text(0.02, 1500.0, " core bottom", transform=ax1.get_yaxis_transform(), va="bottom", fontsize=8, color="dimgrey")
    ax1.text(0.02, 2500.0, " core top", transform=ax1.get_yaxis_transform(), va="bottom", fontsize=8, color="dimgrey")

    axt = fig.add_subplot(gs[1, :])
    axt.axis("off")
    head = ["station", "z (mm)", "IWGFR/75\nlinear theory", "beam 1D-1", "Z3ST 3D", "Z3ST - beam\n(mm)",
            "Z3ST - beam\n(%)", "of which\nthermal shear", "rest", "AVERAGE(10)\n[min, max]"]
    cells = [[r["station"], f"{r['z_mm']:.0f}", f"{r['theory_mm']:.3f}", f"{r['beam_mm']:.3f}", f"{r['z3st_mm']:.3f}",
              f"{r['z3st_minus_beam_mm']:+.4f}", f"{r['z3st_minus_beam_pct']:+.2f} %", f"{r['thermal_shear_mm']:+.4f}",
              f"{r['rest_mm']:+.4f}", f"{r['avg10_mm']:.3f} [{r['avg10_min_mm']:.3f}, {r['avg10_max_mm']:.3f}]"]
             for r in rows]
    widths = [0.08, 0.07, 0.1, 0.09, 0.09, 0.1, 0.1, 0.1, 0.08, 0.19]
    tab = axt.table(cellText=cells, colLabels=head, colWidths=widths, loc="center", cellLoc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1.0, 2.0)
    for j in range(len(head)):
        tab[0, j].set_facecolor("#e8e8e8")
        tab[0, j].set_text_props(weight="bold")
    axt.set_title("Bow w (mm) at the stations", fontsize=11)

    fig.suptitle("IAEA VP1: beam 1D-1 vs Z3ST 3D vs IWGFR/75, example 1", fontsize=14, weight="bold")
    fig.text(0.5, 0.935, f"clamped base; same T(x, y, z), alpha, T_ref; beam {len(z_b) - 1} elements, "
             f"Z3ST {len(z3) - 1} layers; gamma in the core: Z3ST {gamma3 * 1e6:.2f}, thin wall "
             f"{float(np.mean(gamma_cf[core])) * 1e6:.2f} urad", ha="center", fontsize=10)
    path = os.path.join(OUT_DIR, f"{PREFIX}.png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] {path}")

    with open(os.path.join(OUT_DIR, f"{PREFIX}_summary.json"), "w") as f:
        json.dump({"z3st_source": os.path.relpath(xdmf, CASE_DIR), "beam_elements": int(len(z_b) - 1),
                   "thermal_shear_c_m2": c, "thermal_shear_c_over_R2": c / R**2,
                   "gamma_core_z3st_rad": gamma3, "gamma_core_thin_wall_rad": float(np.mean(gamma_cf[core])),
                   "max_abs_z3st_minus_beam_mm": float(np.max(np.abs(w3 - wb)) * 1e3), "rows": rows}, f, indent=2)


if __name__ == "__main__":
    main()
