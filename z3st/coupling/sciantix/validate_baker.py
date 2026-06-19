#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST <-> SCIANTIX binding — Baker 1273 K validation (PROTOTYPE)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Drive the ctypes binding through the SCIANTIX standalone regression case
``regression/baker/test_Baker1977__1273K`` and compare the final-step engineering
outputs against that case's ``output_gold.txt`` (last row). This closes the
binding-correctness gap: same models (read from the case ``input_settings.txt``),
same initial conditions (seeded from ``input_initial_conditions.txt``), same
isothermal history (1273 K, 1e19 fiss/m^3 s) and the same 100 x 55 h time grid as
the standalone output, so the numbers must match.

Run (after building SCIANTIX as a shared lib — see README.md), from the case dir:
    cd <sciantix>/regression/baker/test_Baker1977__1273K
    SCIANTIX_LIB=/path/to/libsciantix.so \
      PYTHONPATH=<z3st>/z3st/coupling/sciantix python3 validate_baker.py
or pass the case dir explicitly: ``python3 validate_baker.py <case_dir>``.
"""

import os
import sys

import numpy as np

from sciantix_binding import SciantixSolver

# Baker 1273 K isothermal history + the exact output grid of output_gold.txt
T_K = 1273.0
FISSION_RATE = 1.0e19      # fiss / m^3 s
N_STEPS = 100              # 0 -> 5500 h in 100 uniform steps (matches the gold)
DT_S = 55.0 * 3600.0       # 55 h per step

# output_gold.txt last data row (t = 5500 h)
GOLD = {
    "fission_gas_release":         0.1320968,
    "intragranular_gas_swelling":  0.0003073483,
    "intergranular_gas_swelling":  0.0416604,
    "burnup_MWd_kgUO2":            6.719293,
}
RTOL = 1.0e-2   # 1% — the binding marches a coarse fixed grid, not SCIANTIX's adaptivity

case_dir = sys.argv[1] if len(sys.argv) > 1 else "."
os.chdir(case_dir)   # getSciantixOptions reads input_settings.txt from the CWD

# Burnup trajectory the host (Z3ST RADAR model) would supply, read from the gold's
# Burnup column. Fed into history[7]/[8] each step: a -DCOUPLING_TU SCIANTIX consumes
# it (and skips its own Burnup/EffectiveBurnup/Densification); a plain build ignores
# those slots and integrates burnup itself — so this harness validates both builds.
bu = np.loadtxt("output_gold.txt", skiprows=1)[:, 31]   # col 32 = Burnup (MWd/kgUO2)

pt = SciantixSolver()
pt.load_initial_conditions("input_initial_conditions.txt")

out = None
for k in range(N_STEPS):
    out = pt.advance(DT_S, T_K, T_K, FISSION_RATE,
                     burnup_old=bu[k], burnup_new=bu[k + 1])

print(f"\nBaker 1273 K — binding vs standalone gold (t = 5500 h, {N_STEPS} steps)\n")
print(f"{'quantity':<30}{'binding':>16}{'gold':>16}{'rel.err':>12}")
all_pass = True
for key, gold in GOLD.items():
    val = out[key]
    rel = abs(val - gold) / abs(gold) if gold else abs(val)
    status = "" if rel < RTOL else "   <-- FAIL"
    if rel >= RTOL:
        all_pass = False
    print(f"{key:<30}{val:>16.6g}{gold:>16.6g}{rel:>12.2e}{status}")

print(f"\n[VALIDATION] {'PASS' if all_pass else 'FAIL'} (rtol = {RTOL:.0e})")
sys.exit(0 if all_pass else 1)
