#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: sandbox/msr_corrosion/V1_erfc_thermal -- V1_erfc with D evaluated at the
temperature field of the thermal model (uniform 973 K). The Cr field and the mass
loss must match V1_erfc, which uses a constant D = 1e-16 m2/s.
"""

import glob
import os
import sys

import numpy as np

from z3st.utils.utils_extract_vtu import extract_field
from z3st.utils.utils_verification import pass_fail_check

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
REF = os.path.join(CASE_DIR, "..", "V1_erfc", "output")

mine = sorted(glob.glob(os.path.join(OUT, "fields_*.vtu")))
ref = sorted(glob.glob(os.path.join(REF, "fields_*.vtu")))
_, _, _, x = extract_field(mine[-1], field_name="CrFraction")
_, _, _, x_ref = extract_field(ref[-1], field_name="CrFraction")
_, _, _, T = extract_field(mine[-1], field_name="Temperature")
m = np.loadtxt(os.path.join(OUT, "corrosion.tsv"), skiprows=1, ndmin=2)[-1, 1]
m_ref = np.loadtxt(os.path.join(REF, "corrosion.tsv"), skiprows=1, ndmin=2)[-1, 1]

errors = {
    "step_count": {"rel_error": float(len(mine) != len(ref))},
    "temperature": {"rel_error": float(np.max(np.abs(np.ravel(T) - 973.0)) / 973.0)},
    "cr_field": {"rel_error": float(np.max(np.abs(np.ravel(x) - np.ravel(x_ref))) / 0.2)},
    "mass_loss": {"rel_error": float(abs(m - m_ref) / m_ref)},
}
ok = pass_fail_check(errors, 1.0e-6, os.path.join(OUT, "non-regression.json"), CASE_DIR)
sys.exit(0 if ok else 1)
