#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: U_pwr_rod_2D

Protects the full PCMI coupling chain (thermal + mechanical + Gas gap
conductance + penalty contact + burnup state bus + swelling eigenstrain bus)
against regression. All metrics are read from ``output/history.csv`` — the
per-step trajectory streamed by the case-local ``diagnostics.py`` — so the
check is independent of the output format (the run writes a single XDMF).

One metric has a closed form and is checked analytically:

  * ``burnup_avg_final`` — the nodal-mean burnup over the fuel equals the flat
    closed form  bu = Σ_k lhr_k·Δt_k / (area·ρ·HM·8.64e10), because the
    radial form factor is area-normalised to mean 1 (the source bus preserves
    the rating) and ``update_state`` accumulates with the right-endpoint rule
    over the generated power history.

The PCMI end-state scalars (gap, contact pressure, temperatures) have no
closed form — they are recorded with ``rel_error = 0`` so the analytic
pass/fail gate ignores them, and are protected purely by the gold comparison
(``regression_check`` vs ``output/non-regression_gold.json``).
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
gap_um = float(last["gap_um"])
p_mpa = float(last["contact_pressure_MPa"])
t_max_final = float(last["T_max_K"])
t_max_peak = max(float(r["T_max_K"]) for r in rows)

# --. closed-form mean burnup from the input power history --..
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
fuel_card = inp["materials"]["cyl_1"]
with open(os.path.join(CASE_DIR, fuel_card)) as f:
    fuel = yaml.safe_load(f)

Ro = float(geom["outer_radius_1"])
area = np.pi * Ro**2
rho = float(fuel["rho"])
hm = float(fuel.get("heavy_metal_fraction", 0.8815))
SECONDS_PER_MWD = 8.64e10  # 86400 s/day * 1e6 W/MW

# Replicate the solver's stepping: __main__ generates the discretised power
# history and update_state accumulates Δbu_k = q(lhr_k)·Δt_k (right endpoint).
times, lhrs, _ = generate_power_history(
    inp["time"], inp["lhr"], n_steps=int(inp["n_steps"]) - 1, filename=None
)
energy_per_m = float(np.sum(np.asarray(lhrs)[1:] * np.diff(np.asarray(times))))
BU_REF = energy_per_m / (area * rho * hm * SECONDS_PER_MWD)

print(f"[INFO] final mean burnup : numerical = {bu_avg:.4f}, "
      f"closed form = {BU_REF:.4f} MWd/kgU")
print(f"[INFO] final peak burnup : {bu_max:.4f} MWd/kgU (rim)")
print(f"[INFO] final gap         : {gap_um:.4f} um "
      f"({'closed — PCMI active' if gap_um <= 0 else 'open'})")
print(f"[INFO] contact pressure  : {p_mpa:.4f} MPa")
print(f"[INFO] T_max final / peak: {t_max_final:.2f} / {t_max_peak:.2f} K")

if gap_um > 0:
    print("[WARNING] gap did not close by end of run — PCMI path not exercised.")

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
    "gap_final_um": _regression_only(gap_um),
    "contact_pressure_final_MPa": _regression_only(p_mpa),
    "T_max_final_K": _regression_only(t_max_final),
    "T_max_peak_K": _regression_only(t_max_peak),
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
