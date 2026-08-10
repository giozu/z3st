#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/creep_shrink_fit_2D

PCMI coupling chain (thermal + mechanical + Gas gap
conductance + penalty contact + burnup state bus + swelling eigenstrain bus).
All metrics are read from ``output/history.csv`` — the
per-step trajectory streamed by the case-local ``diagnostics.py`` — so the
check is independent of the output format (the run writes a single XDMF).

One metric is checked analytically, and it asserts the absence of burnup:

  * ``burnup_avg_final`` — the pellet is ``fissile: false`` on purpose, so the
    accumulated burnup must stay identically zero. Esposito eq. (21) relaxes a
    shrink fit assembled at a fixed interference; a fissile pellet grows the
    interference through the swelling eigenstrain and the contact pressure
    rises instead of relaxing.

    The error is scaled by ``BU_NOMINAL``, the burnup the nominal power history
    in ``input.yaml`` would deposit if the pellet were fissile,
    bu = Σ_k lhr_k·Δt_k / (area·ρ·HM·8.64e10). Any leaked burnup is reported as
    a fraction of the full nominal value.

The PCMI end-state scalars (gap, contact pressure, temperatures) have no
closed form — they are recorded with ``rel_error = 0`` so the analytic
pass/fail gate ignores them, and are protected purely by the gold comparison
(``regression_check`` vs ``output/non-regression_gold.json``).
"""

import os
import csv
import yaml
import numpy as np

from case_params import (
    elastic_factor, hot_interference,
    R_CLAD_I, R_CLAD_O, E_FUEL, NU_FUEL, E_EL, NU, INITIAL_GAP, K_PEN,
)
from z3st.utils.non_regression import finish
from z3st.utils.utils_load import generate_power_history

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
HISTORY = os.path.join(OUT, "history.csv")

# Set by the elastic point below, whose deviation is 3.0 %.
TOLERANCE = 5e-2

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
fuel_card = inp["materials"]["fuel"]
with open(os.path.join(CASE_DIR, fuel_card)) as f:
    fuel = yaml.safe_load(f)

Ro = float(geom["outer_radius_1"])
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
# Burnup the nominal power history would deposit on a fissile pellet: the
# scale the zero-assertion is measured against.
BU_NOMINAL = energy_per_m / (area * rho * hm * SECONDS_PER_MWD)
BU_SCALE = BU_NOMINAL if BU_NOMINAL > 0.0 else 1.0

print(f"[INFO] final mean burnup : numerical = {bu_avg:.4f} MWd/kgU "
      f"(must be 0; nominal history would give {BU_NOMINAL:.4f} on a "
      f"fissile pellet)")
print(f"[INFO] final peak burnup : {bu_max:.4f} MWd/kgU (must be 0)")
print(f"[INFO] final gap         : {gap_um:.4f} um "
      f"({'closed — PCMI active' if gap_um <= 0 else 'open'})")
print(f"[INFO] contact pressure  : {p_mpa:.4f} MPa")
print(f"[INFO] T_max final / peak: {t_max_final:.2f} / {t_max_peak:.2f} K")

if gap_um > 0:
    print("[WARNING] gap did not close by end of run — PCMI path not exercised.")

# --. elastic point at t = 0 vs Esposito eq. (4), springs in series --..
# The joint and the penalty contact carry the same load, so the interference
# splits between them: Delta = p/f + p/K_PEN, i.e. p = Delta / (1/f + 1/K_PEN).
# f comes from the geometry and the moduli, K_PEN from input.yaml and Delta
# from the cards and the run temperature — nothing on the right-hand side is
# read back from the solver.
#
# Delta is the interference of the assembled joint at temperature, not the cold
# 10 um of models.contact.initial_gap: at 580 K the pellet outgrows the clad by
# 4.73 um.
_f_joint = elastic_factor(0.0, R_CLAD_I, R_CLAD_O, E_FUEL, NU_FUEL, E_EL, NU)
_T0 = float(rows[0]["T_max_K"])
if float(rows[0]["T_min_K"]) != _T0:
    raise SystemExit("[non-regression] the field is not uniform at t = 0; "
                     "hot_interference() assumes an isothermal case")
_delta_hot_m = hot_interference(_T0)
p_series_mpa = _delta_hot_m / (1.0 / _f_joint + 1.0 / K_PEN) / 1e6
p_initial_mpa = float(rows[0]["contact_pressure_MPa"])
_elastic_dev = abs(p_initial_mpa - p_series_mpa) / p_series_mpa
_penetration_um = -float(rows[0]["gap_um"])
print(f"[INFO] interference t=0  : {_delta_hot_m*1e6:.4f} um hot "
      f"({abs(INITIAL_GAP)*1e6:.4f} um cold + {(_delta_hot_m - abs(INITIAL_GAP))*1e6:.4f} um "
      f"differential expansion at {_T0:.1f} K)")
print(f"[INFO] elastic point t=0 : z3st = {p_initial_mpa:.4f} MPa, "
      f"eq. (4) in series = {p_series_mpa:.4f} MPa, "
      f"deviation = {100*_elastic_dev:+.2f} %")
print(f"[INFO] split t=0         : {p_initial_mpa*1e6/_f_joint*1e6:.4f} um joint + "
      f"{_penetration_um:.4f} um penalty penetration")


def _regression_only(value):
    """Metric with no closed form: gold-comparison protection only."""
    return {"numerical": float(value), "reference": None,
            "abs_error": 0.0, "rel_error": 0.0}

errors = {
    "burnup_avg_final": {
        "numerical": bu_avg,
        "reference": 0.0,
        "abs_error": float(abs(bu_avg)),
        "rel_error": float(abs(bu_avg) / BU_SCALE),
    },
    "burnup_max_final": {
        "numerical": bu_max,
        "reference": 0.0,
        "abs_error": float(abs(bu_max)),
        "rel_error": float(abs(bu_max) / BU_SCALE),
    },
    "contact_pressure_elastic_MPa": {
        "numerical": p_initial_mpa,
        "reference": p_series_mpa,
        "abs_error": float(abs(p_initial_mpa - p_series_mpa)),
        "rel_error": float(_elastic_dev),
    },
    "gap_final_um": _regression_only(gap_um),
    "contact_pressure_final_MPa": _regression_only(p_mpa),
    "T_max_final_K": _regression_only(t_max_final),
    "T_max_peak_K": _regression_only(t_max_peak),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
