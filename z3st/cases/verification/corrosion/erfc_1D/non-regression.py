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
"""

import glob
import os

import numpy as np
from scipy.special import erfinv

from z3st.core.spine import Spine
from z3st.models.corrosion_model import erfc_profile, profile_mass_loss
from z3st.utils.non_regression import case_paths, error_metric, finish, metric
from z3st.utils.utils_extract_vtu import extract_field

CASE_DIR, _, OUT_JSON = case_paths(__file__)
OUT = os.path.join(CASE_DIR, "output")

D, X0, XS, LX, T_END, RHO, THR = 1.0e-16, 0.2, 0.0, 150.0e-6, 3.6e6, 8570.0, 0.843
TOLERANCE = 2.0e-2

vtus = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
x_coord, _, _, x_cr0 = extract_field(vtus[0], field_name="CrFraction")
x_coord, _, _, x_cr = extract_field(vtus[-1], field_name="CrFraction")
x_cr0, x_cr = np.ravel(x_cr0), np.ravel(x_cr)

hist = np.loadtxt(os.path.join(OUT, "corrosion.tsv"), skiprows=1, ndmin=2)
m_h, depth_h, flux_h = hist[:, 1], hist[:, 3], hist[:, 4]

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
