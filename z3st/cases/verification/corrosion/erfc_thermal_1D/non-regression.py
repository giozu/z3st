#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/corrosion/erfc_thermal_1D -- erfc_1D with D evaluated at the
temperature field of the thermal model (uniform 973 K, Arrhenius D with
D(973 K) = 1e-16 m2/s).

The temperature must stay at 973 K, and the Cr field and the mass loss must match
the same exact erfc solution as erfc_1D: the coupling may not add any error.

Figures in output/: profiles.png (temperature, and the Cr profile against the exact
solution at 1000 h) and errors.png (error profile at 1000 h, mass-loss error history).
"""

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from z3st.models.corrosion_model import erfc_profile, profile_mass_loss
from z3st.utils.non_regression import case_paths, error_metric, finish, metric
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR, _, OUT_JSON = case_paths(__file__)
OUT = os.path.join(CASE_DIR, "output")

D, X0, XS, LX, T_END, RHO, T0 = 1.0e-16, 0.2, 0.0, 150.0e-6, 3.6e6, 8570.0, 973.0
TOLERANCE = 1.0e-2

vtus = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
x_coord, _, _, x_cr = extract_field(vtus[-1], field_name="CrFraction")
_, _, _, T = extract_field(vtus[-1], field_name="Temperature")
x_coord, x_cr, T = np.ravel(x_coord), np.ravel(x_cr), np.ravel(T)
m_end = np.loadtxt(os.path.join(OUT, "corrosion.tsv"), skiprows=1, ndmin=2)[-1, 1]

T_err = np.max(np.abs(T - T0))
profile_err = np.max(np.abs(x_cr - erfc_profile(LX - x_coord, T_END, D, X0, XS)))
errors = {
    "temperature": error_metric(T_err, T_err / T0),
    "profile_1000h": error_metric(profile_err, profile_err / X0),
    "mass_loss_1000h": metric(m_end, profile_mass_loss(T_END, D, X0, XS, RHO)),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)

# --.. ..- .-.. .-.. --- figure --.. ..- .-.. .-.. ---
UM = 1.0e6
o = np.argsort(LX - x_coord)
d = (LX - x_coord)[o]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
a1.plot(d * UM, T[o], color="#D55E00", ls="-", marker="none", lw=1.8, label="thermal model")
a1.axhline(T0, color="k", ls=":", lw=1, marker="none", label=f"{T0:.0f} K (both ends)")
a1.set_ylim(T0 - 2, T0 + 2)
a1.set_ylabel(r"Temperature $T$ (K)")
a1.legend(loc="upper right", fontsize=8)
dfine = np.linspace(0.0, LX, 400)
a2.plot(dfine * UM, erfc_profile(dfine, T_END, D, X0, XS), color="k", ls="-", marker="none", lw=1.8,
        label="exact erfc, 1000 h")
a2.plot(d[::8] * UM, x_cr[o][::8], color="#0072B2", ls="none", marker="o", ms=4, mfc="white",
        label=r"z3st, $D(T)$ from the thermal field")
a2.set_xlim(0, 100)
a2.set_ylabel(r"Cr atomic fraction $x$")
a2.legend(loc="lower right", fontsize=8)
for ax, lab in ((a1, "(a)"), (a2, "(b)")):
    ax.set_xlabel("Depth below the salt-facing surface (µm)")
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "profiles.png"), dpi=200)
plt.close(fig)
# errors against the exact solution: along the slab at 1000 h, and the mass loss through time
hist = np.loadtxt(os.path.join(OUT, "corrosion.tsv"), skiprows=1, ndmin=2)
k = hist[:, 0] > 0
tk, mk = hist[k, 0], hist[k, 1]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
a1.plot(d * UM, (x_cr[o] - erfc_profile(d, T_END, D, X0, XS)) / X0 * 100, color="#0072B2", ls="-",
        marker="none", lw=1.6, label="Cr, 1000 h")
a1.axhline(0, color="k", lw=0.8, marker="none")
a1.set_xlim(0, 100)
a1.set_xlabel("Depth below the salt-facing surface (µm)")
a1.set_ylabel(r"$(x - x_\mathrm{exact})\,/\,x_0$ (%)")
a1.legend(loc="upper right", fontsize=8)
a2.loglog(tk / 3600.0, np.abs(mk / np.array([profile_mass_loss(t, D, X0, XS, RHO) for t in tk]) - 1) * 100,
          color="#0072B2", ls="-", marker="none", lw=1.6, label="mass loss")
a2.axhline(TOLERANCE * 100, color="k", ls=":", lw=1, marker="none", label=f"tolerance {TOLERANCE * 100:.0f} % (checked at 1000 h)")
a2.set_xlabel(r"Time $t$ (h)")
a2.set_ylabel("Relative error vs exact (%)")
a2.legend(loc="upper right", fontsize=8)
for ax, lab in ((a1, "(a)"), (a2, "(b)")):
    ax.text(-0.15, 1.04, lab, transform=ax.transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "errors.png"), dpi=200)
plt.close(fig)
print("[INFO] figures: output/profiles.png, output/errors.png")

