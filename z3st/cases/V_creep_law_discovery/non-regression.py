#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_creep_law_discovery

Verifies the sparse identification of the creep mechanism from the FEM
stress-relaxation data (see discover.py):
  1. the selection recovers the cubic Norton term, and only that term;
  2. the identified coefficient matches the true A*sigma_ref^3 within the
     tolerance set by the 2% data noise.
"""

import os
import json
import numpy as np

from z3st.utils.utils_verification import pass_fail_check, regression_check

R_GAS = 8.314462618
CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

A_TRUE = 2.82e-24 * np.exp(-1.2e5 / (R_GAS * 600.0))
C3_TRUE = A_TRUE * (5.0e7) ** 3

with open(os.path.join(OUT, "discovery.json")) as f:
    disc = json.load(f)

selected_cubic_only = 1.0 if disc["selected"] == ["S^3"] else 0.0
c3 = disc["coefficients_1_per_s"].get("S^3", 0.0)

TOLERANCE = 1e-1

errors = {
    "selected_cubic_norton_only": {
        "numerical": selected_cubic_only,
        "reference": 1.0,
        "abs_error": float(abs(selected_cubic_only - 1.0)),
        "rel_error": float(abs(selected_cubic_only - 1.0)),
    },
    "identified_cubic_coefficient": {
        "numerical": float(c3),
        "reference": float(C3_TRUE),
        "abs_error": float(abs(c3 - C3_TRUE)),
        "rel_error": float(abs(c3 - C3_TRUE) / C3_TRUE),
    },
}

print(f"[INFO] selected: {disc['selected']}")
print(f"[INFO] c[S^3] = {c3:.4e} 1/s (true {C3_TRUE:.4e} 1/s)")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
