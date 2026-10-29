#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: sandbox/msr_corrosion/V1_erfc -- Cr depletion of a semi-infinite alloy.

Constant D = 1e-16 m2/s, x0 = 0.2, surface (xmax) held at x_surf = 0. The exact
solution at depth d = xmax - x is

    x(d, t) = x_surf + (x0 - x_surf) erf(d / (2 sqrt(D t)))

At t = 1000 h, 2 sqrt(D t) = 38 um against a 150 um slab, so the far end sees
erfc(3.9) ~ 1e-8 and the slab is semi-infinite to plotting accuracy.
"""

import glob
import os
import sys

import numpy as np
from scipy.special import erfinv

from z3st.core.spine import Spine
from z3st.models.corrosion_model import erfc_profile, profile_mass_loss
from z3st.utils.utils_extract_vtu import extract_field
from z3st.utils.utils_verification import pass_fail_check

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

D, X0, XS, LX, T_END, RHO, THR = 1.0e-16, 0.2, 0.0, 150.0e-6, 3.6e6, 8570.0, 0.843

vtus = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
x_coord, _, _, x_cr0 = extract_field(vtus[0], field_name="CrFraction")
x_coord, _, _, x_cr = extract_field(vtus[-1], field_name="CrFraction")
x_cr0, x_cr = np.ravel(x_cr0), np.ravel(x_cr)

errors = {}

# t = 0 is the dt = 0 step: the field must still be x0 everywhere
errors["initial_field"] = {"rel_error": float(np.max(np.abs(x_cr0 - X0)) / X0)}

# Final profile against erfc
x_ref = erfc_profile(LX - x_coord, T_END, D, X0, XS)
errors["profile_1000h"] = {"rel_error": float(np.max(np.abs(x_cr - x_ref)) / X0)}

# Scalar history: mass loss and depletion depth
hist = np.loadtxt(os.path.join(OUT, "corrosion.tsv"), skiprows=1, ndmin=2)
t_h, m_h, depth_h = hist[:, 0], hist[:, 1], hist[:, 3]
errors["mass_loss_t0"] = {"rel_error": float(abs(m_h[0]))}
m_ref = profile_mass_loss(T_END, D, X0, XS, RHO)  # 3.3587e-2 kg/m2
errors["mass_loss_1000h"] = {"rel_error": float(abs(m_h[-1] - m_ref) / m_ref)}
depth_ref = 2.0 * np.sqrt(D * T_END) * erfinv(THR)  # erfc_profile = THR x0
depth_err = float(abs(depth_h[-1] - depth_ref) / depth_ref)
flux_ref = D * (X0 - XS) / np.sqrt(np.pi * D * T_END)  # -D dx/dn at the surface
flux_err = float(abs(hist[-1, 4] - flux_ref) / flux_ref)
print(f"mass loss {m_h[-1]:.5e} kg/m2 (ref {m_ref:.5e}), depth {depth_h[-1]*1e6:.3f} um (ref {depth_ref*1e6:.3f})")

# Rollback roster: an adaptive bisection must restore the Cr field
in_roster = {"x_cr", "x_cr_n"} <= set(Spine._SNAPSHOT_FIELDS)
errors["snapshot_roster"] = {"rel_error": 0.0 if in_roster else 1.0}

ok = pass_fail_check(errors, 1.0e-2, OUT_JSON, CASE_DIR)

# Depth is read on nodes (0.3 um spacing at 20 um): 2 %
depth_ok = depth_err < 2.0e-2
print(f"  depth_1000h        → rel err = {depth_err:.2e}  → {'PASS' if depth_ok else 'FAIL'} (tolerance 2.0e-02)")
# Flux from the P1 gradient of the surface element lags like the profile: 2 %
flux_ok = flux_err < 2.0e-2
print(f"  flux_1000h         → rel err = {flux_err:.2e}  → {'PASS' if flux_ok else 'FAIL'} (tolerance 2.0e-02)")
ok = ok and depth_ok and flux_ok
sys.exit(0 if ok else 1)
