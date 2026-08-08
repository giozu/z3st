#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/creep_law_discovery

Verifies the sparse identification of the creep mechanism from the FEM
stress-relaxation data (see discover.py):
  1. the selection recovers the cubic Norton term, and only that term;
  2. the identified coefficient matches the true A*sigma_ref^3 within the
     tolerance set by the 2% data noise.
"""

import os
import json
import numpy as np

from z3st.utils.non_regression import finish, metric


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

# robustness over independent noise realisations (discover.py seed sweep):
# the cubic Norton term must be recovered alone in every realisation, and the
# worst-case coefficient error must stay within tolerance.
rob = disc.get("robustness", {})
n_seeds = int(rob.get("n_seeds", 0))
frac_cubic_only = (rob.get("n_cubic_only", 0) / n_seeds) if n_seeds else 0.0
max_c3_err = rob.get("max_c3_rel_error", None)

TOLERANCE = 1e-1

errors = {
    "selected_cubic_norton_only": metric(selected_cubic_only, 1.0, rel=abs(selected_cubic_only - 1.0)),
    "identified_cubic_coefficient": metric(c3, C3_TRUE),
    "robustness_cubic_only_fraction": metric(frac_cubic_only, 1.0, rel=abs(frac_cubic_only - 1.0)),
}
if max_c3_err is not None:
    errors["robustness_max_cubic_coeff_error"] = {
        "numerical": float(max_c3_err),
        "reference": 0.0,
        "abs_error": float(max_c3_err),
        "rel_error": float(max_c3_err),
    }

print(f"[INFO] selected: {disc['selected']}")
print(f"[INFO] c[S^3] = {c3:.4e} 1/s (true {C3_TRUE:.4e} 1/s)")
print(f"[INFO] robustness: cubic-only in {rob.get('n_cubic_only', 0)}/{n_seeds} "
      f"seeds, max coeff error {(max_c3_err or 0.0) * 100:.2f}%")

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
