#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/corrosion/erfc_thermal_1D -- erfc_1D with D evaluated at the
temperature field of the thermal model (uniform 973 K, Arrhenius D with
D(973 K) = 1e-16 m2/s).

The temperature must stay at 973 K, and the Cr field and the mass loss must match
the same exact erfc solution as erfc_1D: the coupling may not add any error.
"""

import glob
import os

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
