#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: regression/fg_test_fuel

Bare fuel pellet with fission-gas behaviour (thermal + mechanical + burnup state
bus + SCIANTIX swelling eigenstrain bus). All metrics are read from
``output/history.csv`` — the per-step trajectory streamed by the case-local
``diagnostics.py`` — so the check is independent of the output format.

Analytic check (closed form):

  * ``burnup_avg_final`` — nodal-mean fuel burnup equals
    bu = Σ_k lhr_k·Δt_k / (area·ρ·HM·8.64e10) (radial form factor area-normalised
    to mean 1; ``update_state`` uses the right-endpoint rule).

The other scalars (peak burnup, temperatures, gas release) have no closed form:
recorded with ``rel_error = 0`` so the analytic gate ignores them, protected by
the gold comparison only.
"""

import os
import csv
import yaml
import numpy as np

from z3st.utils.utils_verification import pass_fail_check, regression_check
from z3st.utils.utils_load import generate_power_history

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
HISTORY = os.path.join(OUT, "history.csv")

TOLERANCE = 1e-2

# --. trajectory from the case-local diagnostics CSV --..
with open(HISTORY) as f:
    rows = list(csv.DictReader(f))
if not rows:
    raise RuntimeError(f"{HISTORY} is empty — did the run complete?")
last = rows[-1]

bu_avg = float(last["burnup_avg_MWdkgU"])
bu_max = float(last["burnup_max_MWdkgU"])
t_max_final = float(last["T_max_K"])
t_max_peak = max(float(r["T_max_K"]) for r in rows)
fgr = float(last["fgr_frac"])

# --. closed-form mean burnup from the input power history --..
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
fuel_card = inp["materials"]["fuel"]
with open(os.path.join(CASE_DIR, fuel_card)) as f:
    fuel = yaml.safe_load(f)

Ro = float(geom["Ro"])
area = np.pi * Ro**2
rho = float(fuel["rho"])
hm = float(fuel.get("heavy_metal_fraction", 0.8815))
SECONDS_PER_MWD = 8.64e10  # 86400 s/day * 1e6 W/MW

raw_n_steps = inp["n_steps"]
n_increments = raw_n_steps if isinstance(raw_n_steps, (list, tuple)) else int(raw_n_steps) - 1
times, lhrs, _ = generate_power_history(
    inp["time"], inp["lhr"], n_steps=n_increments, filename=None
)
energy_per_m = float(np.sum(np.asarray(lhrs)[1:] * np.diff(np.asarray(times))))
BU_REF = energy_per_m / (area * rho * hm * SECONDS_PER_MWD)

print(f"[INFO] final mean burnup : numerical = {bu_avg:.4f}, "
      f"closed form = {BU_REF:.4f} MWd/kgU")
print(f"[INFO] final peak burnup : {bu_max:.4f} MWd/kgU (rim)")
print(f"[INFO] T_max final / peak: {t_max_final:.2f} / {t_max_peak:.2f} K")
print(f"[INFO] fractional gas release: {fgr:.6f}")

def _regression_only(value):
    """Metric with no closed form: gold-comparison protection only."""
    return {"numerical": float(value), "reference": None,
            "abs_error": 0.0, "rel_error": 0.0}

errors = {
    "burnup_avg_final": {
        "numerical": bu_avg,
        "reference": BU_REF,
        "abs_error": float(abs(bu_avg - BU_REF)),
        "rel_error": float(abs(bu_avg - BU_REF) / BU_REF),
    },
    "burnup_max_final": _regression_only(bu_max),
    "T_max_final_K": _regression_only(t_max_final),
    "T_max_peak_K": _regression_only(t_max_peak),
    "fgr_final": _regression_only(fgr),
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
