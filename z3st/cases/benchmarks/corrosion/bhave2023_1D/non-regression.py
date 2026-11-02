#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: benchmarks/corrosion/bhave2023_1D -- Bhave et al. (2023) 1D benchmark.

C. V. Bhave et al., J. Nucl. Mater. 574 (2023) 154147. Ni-5Cr and Ni-20Cr in FLiBe
at 973 K, Hart effective Cr diffusivity (Eqs. 44-46, constants of the authors' inputs),
surface x_s = 1.7726e-3 atomic (the equilibrium value of the authors' inputs, quoted in
the paper as c_surf = 1.8e-3). Bhave show their 1D phase-field model follows Eq. 50 within 5 % (Sec. 4).

Checks, per alloy, at 50 / 500 / 1000 h: mass loss against the exact erfc-profile
integral (the model's own solution), 1 %. The depletion depth at 500 h is tracked
for the gold.
Reported only: the difference to Bhave Eq. 50 (linear in the mass fraction, so
about 0.7 % below the exact integral), the experiments of Table 3, the Ni-5Cr
depth at 500 h (18 um in Bhave Fig. 2) and the redox surface value from Bhave's
Table 1 data.

Figures in output/: mass_loss.png (z3st, Eq. 50 and experiments, per alloy; labels:
z3st / experiment) and
profiles.png (Cr profiles at 50, 500 and 1000 h, per alloy, with the depth).
"""

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.models.corrosion_model import (
    eq50_mass_loss,
    hart_diffusivity,
    mass_fraction_cr,
    profile_mass_loss,
    redox_surface_fraction,
)
from z3st.utils.non_regression import finish, metric, tracked
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
os.makedirs(OUT, exist_ok=True)
OUT_JSON = os.path.join(OUT, "non-regression.json")

T, XS, THR = 973.0, 1.7726e-3, 0.843
HOURS = np.array([50.0, 500.0, 1000.0])
MG_CM2 = 100.0  # kg/m2 -> mg/cm2

# Experiments of Table 3 (mg/cm2)
ALLOYS = {
    "ni5cr": {"label": "Ni-5Cr", "exp": [0.09, 0.27, 0.42]},
    "ni20cr": {"label": "Ni-20Cr", "exp": [0.67, 3.06, 7.12]},
}

errors = {}
res = {}  # per alloy, for the figures

for key, ref in ALLOYS.items():
    with open(os.path.join(CASE_DIR, key, "alloy.yaml")) as f:
        card = {k: float(v) for k, v in yaml.safe_load(f)["corrosion"].items()}
    D = hart_diffusivity(card, T)
    hist = np.loadtxt(os.path.join(CASE_DIR, key, "output", "corrosion.tsv"), skiprows=1, ndmin=2)
    t_h, m_h, depth_h = hist[:, 0], hist[:, 2], hist[:, 3]

    idx = [int(np.argmin(np.abs(t_h - h * 3600.0))) for h in HOURS]
    assert np.allclose(t_h[idx], HOURS * 3600.0), "50/500/1000 h missing from the time grid"
    m_model = m_h[idx]
    m_exact = np.array([profile_mass_loss(h * 3600.0, D, card["x0"], XS, card["rho"]) for h in HOURS]) * MG_CM2
    # Eq. 50 as published: mass fractions, c_surf = 1.8e-3
    m_eq50 = np.array([eq50_mass_loss(mass_fraction_cr(card["x0"]), 1.8e-3, card["rho"], D, h * 3600.0)
                       for h in HOURS]) * MG_CM2
    m_exp = np.array(ref["exp"])

    for i, h in enumerate(HOURS):
        errors[f"{key}_mass_loss_{h:.0f}h"] = metric(m_model[i], m_exact[i])
    errors[f"{key}_depth_500h"] = tracked(depth_h[idx[1]])

    print(f"\n{ref['label']}: D_eff = {D:.4e} m2/s, w0 = {mass_fraction_cr(card['x0']):.4f}")
    print("   t (h)   model    exact    Eq.50   vs Eq.50   exp    model/exp   (mg/cm2)")
    for i, h in enumerate(HOURS):
        print(f"  {h:6.0f}  {m_model[i]:7.4f}  {m_exact[i]:7.4f}  {m_eq50[i]:7.4f}  "
              f"{100 * (m_model[i] / m_eq50[i] - 1):+6.2f} %  {m_exp[i]:6.2f}  {m_model[i] / m_exp[i]:6.2f}")
    if key == "ni5cr":
        print(f"  depth at 500 h: {depth_h[idx[1]] * 1e6:.1f} um (Bhave Fig. 2: 18 um, with their D)")

    res[key] = dict(card=card, D=D, t=t_h, m=m_h, depth=depth_h, idx=idx, m_model=m_model, m_exp=m_exp)

x_redox, _ = redox_surface_fraction(T, E_F=3.3607, E0_F=2.871)
print(f"\nSurface Cr: fixed {XS:.4e} (authors' inputs, atomic; Eq. 50 uses c_surf = 1.8e-3); "
      f"redox with the same salt state: {x_redox:.4e}")


# --.. ..- .-.. .-.. --- figures --.. ..- .-.. .-.. ---
COL = {"ni5cr": "#0072B2", "ni20cr": "#D55E00"}
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
t_c = np.linspace(0.0, HOURS[-1], 400)
for (key, ref), ax, lab in zip(ALLOYS.items(), axes, ("(a)", "(b)")):
    r = res[key]
    w0 = mass_fraction_cr(r["card"]["x0"])
    ax.plot(t_c, [eq50_mass_loss(w0, 1.8e-3, r["card"]["rho"], r["D"], h * 3600.0) * MG_CM2 for h in t_c],
            color="k", ls="--", marker="none", lw=1.2, label="Bhave Eq. 50")
    ax.plot(r["t"] / 3600.0, r["m"], color=COL[key], ls="-", marker="none", lw=2, label="z3st")
    ax.plot(HOURS, r["m_exp"], color="k", ls="none", marker="D", ms=6, mfc="white", label="experiment (Zheng)")
    for h, m_mod, m_x in zip(HOURS, r["m_model"], r["m_exp"]):
        ax.annotate(f"{m_mod / m_x:.2f}×", (h, m_x), (-8 if h == HOURS[-1] else 6, -13),
                    textcoords="offset points", fontsize=8, color="0.35", ha="right" if h == HOURS[-1] else "left")
    ax.set_xlim(0, 1080)
    ax.set_ylim(0, 1.12 * max(r["m"].max(), r["m_exp"].max()))
    ax.set_xlabel(r"Time $t$ (h)")
    ax.set_ylabel(r"Cr mass loss (mg cm$^{-2}$)")
    ax.legend(loc="upper left", fontsize=8, title=f"{ref['label']}, FLiBe, 973 K", title_fontsize=8)
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "mass_loss.png"), dpi=200)
plt.close(fig)

BLUES = ("#9ec5f4", "#5a9be6", "#2a6fc0", "#0d366b")
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for (key, ref), ax, lab in zip(ALLOYS.items(), axes, ("(a)", "(b)")):
    r = res[key]
    x0 = r["card"]["x0"]
    vtus = sorted(glob.glob(os.path.join(CASE_DIR, key, "output", "fields_*.vtu")))
    for th, c in zip((1.0, 50.0, 500.0, 1000.0), BLUES):
        i = int(np.argmin(np.abs(r["t"] - th * 3600.0)))
        x, _, _, v = extract_field(vtus[i], field_name="CrFraction")
        d = x.max() - np.ravel(x)
        o = np.argsort(d)
        ax.plot(d[o] * 1e6, np.ravel(v)[o], color=c, ls="-", marker="none", lw=1.6, label=f"{th:.0f} h")
    ax.axhline(THR * x0, color="0.5", ls=":", lw=1, marker="none")
    ax.text(98, THR * x0 * 0.97, r"$0.843\,x_0$", ha="right", va="top", fontsize=8, color="0.35")
    d500 = r["depth"][r["idx"][1]] * 1e6
    ax.plot([d500], [THR * x0], color="k", ls="none", marker="o", ms=4)
    note = f"δ(500 h) = {d500:.1f} µm" + ("\n(Bhave Fig. 2: 18 µm)" if key == "ni5cr" else "")
    ax.annotate(note, (d500, THR * x0), (d500 + 8, 0.45 * x0), fontsize=8,
                arrowprops=dict(arrowstyle="-", color="0.4", lw=0.7))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.08 * x0)
    ax.set_xlabel("Depth below the salt-facing surface (µm)")
    ax.set_ylabel(r"Cr atomic fraction $x$")
    ax.legend(loc="lower right", fontsize=8, title=f"{ref['label']} ($x_0$ = {x0})", title_fontsize=8)
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "profiles.png"), dpi=200)
plt.close(fig)
print("[INFO] figures: output/mass_loss.png, output/profiles.png")

finish(errors, 1.0e-2, OUT_JSON, CASE_DIR)
