#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/corrosion/erfc_1D -- Cr depletion of a semi-infinite alloy.

Constant D = 1e-16 m2/s, x0 = 0.2, surface (xmax) held at x_surf = 0. The exact
solution at depth d = xmax - x is (Crank 1975)

    x(d, t) = x_surf + (x0 - x_surf) erf(d / (2 sqrt(D t)))

At t = 1000 h, 2 sqrt(D t) = 38 um against a 150 um slab, so the far end sees
erfc(3.9) ~ 1e-8 and the slab is semi-infinite to plotting accuracy.

One tolerance, 2 %, for all checks: the depth is read on nodes (0.3 um spacing at
20 um) and the flux comes from the P1 gradient of the surface element, both
within 2 %; the profile and the mass loss are within 0.5 %. The gold comparison
(0.1 % on each value) is the tighter guard.

Figures in output/: profiles.png (profiles against the exact solution, and their
collapse on erf(eta)), kinetics.png (mass loss and depth against sqrt(t)) and
errors.png (error profiles, and the error history of mass loss, depth and flux).
"""

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf, erfinv

from z3st.core.spine import Spine
from z3st.models.corrosion_model import erfc_profile, profile_mass_loss
from z3st.utils.non_regression import case_paths, error_metric, finish, metric
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR, _, OUT_JSON = case_paths(__file__)
OUT = os.path.join(CASE_DIR, "output")

D, X0, XS, LX, T_END, RHO, THR = 1.0e-16, 0.2, 0.0, 150.0e-6, 3.6e6, 8570.0, 0.843
TOLERANCE = 2.0e-2
H, UM, MG = 3600.0, 1.0e6, 100.0  # s/h, um/m, (mg/cm2)/(kg/m2)

vtus = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
x_coord, _, _, x_cr0 = extract_field(vtus[0], field_name="CrFraction")
x_coord, _, _, x_cr = extract_field(vtus[-1], field_name="CrFraction")
x_cr0, x_cr = np.ravel(x_cr0), np.ravel(x_cr)

hist = np.loadtxt(os.path.join(OUT, "corrosion.tsv"), skiprows=1, ndmin=2)
t_h, m_h, depth_h, flux_h = hist[:, 0], hist[:, 1], hist[:, 3], hist[:, 4]

m_ref = profile_mass_loss(T_END, D, X0, XS, RHO)  # 3.3587e-2 kg/m2
depth_ref = 2.0 * np.sqrt(D * T_END) * erfinv(THR)  # where erfc_profile = THR x0
flux_ref = D * (X0 - XS) / np.sqrt(np.pi * D * T_END)  # -D dx/dn at the surface

# Rollback roster: an adaptive bisection must restore the Cr field
in_roster = {"x_cr", "x_cr_n"} <= set(Spine._SNAPSHOT_FIELDS)

profile_err = np.max(np.abs(x_cr - erfc_profile(LX - x_coord, T_END, D, X0, XS)))
errors = {
    # t = 0 is the dt = 0 step: the field must still be x0 everywhere
    "initial_field": error_metric(np.max(np.abs(x_cr0 - X0)), np.max(np.abs(x_cr0 - X0)) / X0),
    "mass_loss_t0": error_metric(abs(m_h[0])),
    "profile_1000h": error_metric(profile_err, profile_err / X0),
    "mass_loss_1000h": metric(m_h[-1], m_ref),
    "depth_1000h": metric(depth_h[-1], depth_ref),
    "flux_1000h": metric(flux_h[-1], flux_ref),
    "snapshot_roster": error_metric(0.0 if in_roster else 1.0),
}
print(f"mass loss {m_h[-1]:.5e} kg/m2 (ref {m_ref:.5e}), depth {depth_h[-1]*1e6:.3f} um "
      f"(ref {depth_ref*1e6:.3f}), flux {flux_h[-1]:.4e} (ref {flux_ref:.4e})")

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)

# --.. ..- .-.. .-.. --- figures --.. ..- .-.. .-.. ---
TIMES = (1.0, 10.0, 100.0, 1000.0)  # h
BLUES = ("#9ec5f4", "#5a9be6", "#2a6fc0", "#0d366b")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
dfine = np.linspace(0.0, 80e-6, 400)
eta = np.linspace(0.0, 3.0, 300)
a2.plot(eta, erf(eta), color="k", ls="-", marker="none", lw=1.8, label=r"$\mathrm{erf}(\eta)$, exact")
for th, c in zip(TIMES, BLUES):
    i = int(np.argmin(np.abs(t_h - th * H)))
    x, _, _, v = extract_field(vtus[i], field_name="CrFraction")
    d, v = LX - np.ravel(x), np.ravel(v)
    o = np.argsort(d)
    d, v = d[o], v[o]
    a1.plot(dfine * UM, erfc_profile(dfine, t_h[i], D, X0, XS), color=c, ls="-", marker="none", lw=1.6,
            label=f"{th:.0f} h")
    near = d <= 80e-6
    a1.plot(d[near][::10] * UM, v[near][::10], color=c, ls="none", marker="o", ms=3.5, mfc="white")
    e = d / (2.0 * np.sqrt(D * t_h[i]))
    a2.plot(e[e <= 3.0][::6], v[e <= 3.0][::6] / X0, color=c, ls="none", marker="o", ms=3.0, mfc="none",
            label=f"z3st, {th:.0f} h")
a1.plot([], [], color="0.4", ls="none", marker="o", ms=3.5, mfc="white", label="z3st (nodes)")
a1.set_xlabel("Depth below the salt-facing surface (µm)")
a1.set_ylabel(r"Cr atomic fraction $x$")
a1.set_xlim(0, 80)
a1.legend(loc="lower right", fontsize=8, title="exact (lines)", title_fontsize=8)
a2.axhline(THR, color="0.5", ls=":", lw=1, marker="none")
a2.axvline(1.0, color="0.5", ls=":", lw=1, marker="none")
a2.text(1.04, 0.35, r"$\eta = 1$: $x/x_0 = \mathrm{erf}(1) = 0.843$" "\n" r"$\Rightarrow \delta = 2\sqrt{Dt}$", fontsize=8)
a2.set_xlabel(r"$\eta = d\,/\,(2\sqrt{Dt})$")
a2.set_ylabel(r"$x\,/\,x_0$")
a2.legend(loc="lower right", fontsize=8)
for ax, lab in ((a1, "(a)"), (a2, "(b)")):
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "profiles.png"), dpi=200)
plt.close(fig)

tf = np.geomspace(1.0, T_END, 300)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
a1.plot(np.sqrt(tf / H), [profile_mass_loss(t, D, X0, XS, RHO) * MG for t in tf], color="k", ls="-",
        marker="none", lw=1.8, label="exact")
a1.plot(np.sqrt(t_h / H)[::6], m_h[::6] * MG, color="#0072B2", ls="none", marker="o", ms=4, mfc="white",
        label="z3st")
a1.set_ylabel(r"Cr mass loss (mg cm$^{-2}$)")
a2.plot(np.sqrt(tf / H), 2.0 * np.sqrt(D * tf) * erfinv(THR) * UM, color="k", ls="-", marker="none", lw=1.8,
        label=r"exact: $2\sqrt{Dt}\,\mathrm{erf}^{-1}(0.843)$")
a2.step(np.sqrt(t_h / H), depth_h * UM, where="post", color="#D55E00", ls="-", marker="none", lw=1.4,
        label="z3st (read on nodes)")
a2.set_ylabel(r"Depletion depth $\delta$ (µm)")
for ax, lab in ((a1, "(a)"), (a2, "(b)")):
    ax.set_xlabel(r"$\sqrt{t}$ (h$^{1/2}$)")
    ax.legend(loc="upper left", fontsize=8)
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "kinetics.png"), dpi=200)
plt.close(fig)
# errors against the exact solution: along the slab, and through time
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
for th, c in zip(TIMES, BLUES):
    i = int(np.argmin(np.abs(t_h - th * H)))
    x, _, _, v = extract_field(vtus[i], field_name="CrFraction")
    d, v = LX - np.ravel(x), np.ravel(v)
    o = np.argsort(d)
    a1.plot(d[o] * UM, (v[o] - erfc_profile(d[o], t_h[i], D, X0, XS)) / X0 * 100, color=c, ls="-",
            marker="none", lw=1.4, label=f"{th:.0f} h")
a1.axhline(0, color="k", lw=0.8, marker="none")
a1.set_xlim(0, 80)
a1.set_xlabel("Depth below the salt-facing surface (µm)")
a1.set_ylabel(r"$(x - x_\mathrm{exact})\,/\,x_0$ (%)")
a1.legend(loc="upper right", fontsize=8)
k = t_h > 0
tk = t_h[k]
e_m = np.abs(m_h[k] / np.array([profile_mass_loss(t, D, X0, XS, RHO) for t in tk]) - 1) * 100
e_f = np.abs(flux_h[k] / (D * (X0 - XS) / np.sqrt(np.pi * D * tk)) - 1) * 100
dk = depth_h[k]
e_d = np.where(dk > 0, np.abs(dk / (2.0 * np.sqrt(D * tk) * erfinv(THR)) - 1) * 100, np.nan)
a2.loglog(tk / H, e_m, color="#0072B2", ls="-", marker="none", lw=1.6, label="mass loss")
a2.loglog(tk / H, e_f, color="#009E73", ls="-", marker="none", lw=1.6, label="surface flux")
a2.loglog(tk / H, e_d, color="#D55E00", ls="none", marker=".", ms=3, label="depth (nodal)")
a2.axhline(TOLERANCE * 100, color="k", ls=":", lw=1, marker="none", label=f"tolerance {TOLERANCE * 100:.0f} % (checked at 1000 h)")
a2.set_xlabel(r"Time $t$ (h)")
a2.set_ylabel("Relative error vs exact (%)")
a2.legend(loc="upper right", fontsize=8)
for ax, lab in ((a1, "(a)"), (a2, "(b)")):
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "errors.png"), dpi=200)
plt.close(fig)
print("[INFO] figures: output/profiles.png, output/kinetics.png, output/errors.png")

