#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: sandbox/msr_corrosion/V2_bhave1D -- Bhave et al. (2023) 1D benchmark.

C. V. Bhave et al., J. Nucl. Mater. 574 (2023) 154147. Ni-5Cr and Ni-20Cr in FLiBe
at 973 K, Hart effective Cr diffusivity (Eqs. 44-46), surface c_surf = 1.8e-3 mass
fraction. Bhave show their 1D phase-field model follows Eq. 50 within 5 % (Sec. 4).

Checks, per alloy, at 50 / 500 / 1000 h:
  - mass loss against the exact erfc-profile integral (the model's own solution): 1 %
  - mass loss against Bhave Eq. 50: 2 %
Reported only: the experiments of Table 3, the Ni-5Cr depth at 500 h (18 um in
Bhave Fig. 2) and the redox surface value from Bhave's Table 1 data.
"""

import os
import sys

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
from z3st.utils.utils_verification import pass_fail_check

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
os.makedirs(OUT, exist_ok=True)
OUT_JSON = os.path.join(OUT, "non-regression.json")

T, XS = 973.0, 2.0315e-3
HOURS = np.array([50.0, 500.0, 1000.0])
MG_CM2 = 100.0  # kg/m2 -> mg/cm2

# Bhave Eq. 50 with Hart D (mg/cm2) and the experiments of Table 3 (mg/cm2)
ALLOYS = {
    "ni5cr": {"label": "Ni-5Cr", "eq50": [0.133, 0.421, 0.596], "exp": [0.09, 0.27, 0.42]},
    "ni20cr": {"label": "Ni-20Cr", "eq50": [1.007, 3.184, 4.503], "exp": [0.67, 3.06, 7.12]},
}

errors_1, errors_2 = {}, {}
fig, ax = plt.subplots(figsize=(6, 4.5))

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
    m_eq50 = np.array(ref["eq50"])
    m_exp = np.array(ref["exp"])

    errors_1[f"{key}_exact"] = {"rel_error": list(np.abs(m_model - m_exact) / m_exact)}
    errors_2[f"{key}_eq50"] = {"rel_error": list(np.abs(m_model - m_eq50) / m_eq50)}

    print(f"\n{ref['label']}: D_eff = {D:.4e} m2/s, w0 = {mass_fraction_cr(card['x0']):.4f}")
    print("   t (h)   model    exact    Eq.50    exp     model/exp   (mg/cm2)")
    for i, h in enumerate(HOURS):
        print(f"  {h:6.0f}  {m_model[i]:7.4f}  {m_exact[i]:7.4f}  {m_eq50[i]:7.4f}  {m_exp[i]:6.2f}  "
              f"{m_model[i] / m_exp[i]:6.2f}")
    if key == "ni5cr":
        print(f"  depth at 500 h: {depth_h[idx[1]] * 1e6:.1f} um (Bhave Fig. 2: 18 um, with their D)")

    line, = ax.plot(t_h / 3600.0, m_h, ls="-", marker="none", lw=2, label=f"{ref['label']} z3st")
    t_c = np.linspace(0.0, HOURS[-1], 400)
    m_c = [eq50_mass_loss(mass_fraction_cr(card["x0"]), 1.8e-3, card["rho"], D, h * 3600.0) * MG_CM2 for h in t_c]
    ax.plot(t_c, m_c, ls="--", marker="none", color="k", lw=1, label="Eq. 50" if key == "ni5cr" else None)
    ax.plot(HOURS, m_exp, ls="none", marker="o", mfc="white", color=line.get_color(),
            label=f"{ref['label']} experiment")

x_redox, _ = redox_surface_fraction(T, E_F=3.36, E0_F=2.871)
print(f"\nSurface Cr: fixed {XS:.4e} (Bhave c_surf = 1.8e-3 mass fraction); "
      f"redox from Bhave Table 1 / App. A data: {x_redox:.4e}")

ax.set_xlabel("Time (h)")
ax.set_ylabel("Cr mass loss (mg/cm$^2$)")
ax.set_title("Bhave et al. 2023 1D benchmark, FLiBe, 973 K")
ax.grid(True, ls=":", alpha=0.6)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(CASE_DIR, "mass_loss.png"), dpi=150)
plt.close(fig)

ok_1 = pass_fail_check(errors_1, 1.0e-2, OUT_JSON, CASE_DIR)
ok_2 = pass_fail_check(errors_2, 2.0e-2, os.path.join(OUT, "non-regression_eq50.json"), CASE_DIR)
sys.exit(0 if ok_1 and ok_2 else 1)
