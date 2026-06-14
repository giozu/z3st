#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST <-> SCIANTIX binding — standalone smoke test (PROTOTYPE)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Drives one SCIANTIX point from Python with no Z3ST involved: a single grain at
a fixed fission rate, temperature ramped up over a few steps, printing gaseous
swelling / burnup / FGR each step. If this runs and the numbers move sensibly
(swelling rises with temperature/burnup), the ctypes binding is sound and the
array layout is correct.

Run (after building SCIANTIX as a shared lib — see README.md):
    SCIANTIX_LIB=/path/to/libsciantix.so python3 smoke_test.py

This deliberately mirrors a SCIANTIX standalone input_history.txt run, so its
output can be diffed against a reference standalone run for validation.
"""

import numpy as np

from sciantix_binding import SciantixSolver

# one point, constant fission rate ~ 1e19 fiss/m^3 s, T ramp 800 -> 1500 K
FISSION_RATE = 1.0e19          # fiss / m^3 s
DT = 1.0e6                     # s per step (~11.6 days)
TEMPS = np.linspace(800.0, 1500.0, 8)   # K

pt = SciantixSolver()
print(f"{'step':>4} {'T (K)':>8} {'bu':>8} {'gas_swell':>11} {'interg':>11} {'FGR':>8}")
T_old = TEMPS[0]
for k, T_new in enumerate(TEMPS):
    out = pt.advance(DT, T_old, T_new, FISSION_RATE)
    print(f"{k:>4} {T_new:8.1f} {out['burnup_MWd_kgUO2']:8.3f} "
          f"{out['gaseous_swelling']:11.3e} {out['intergranular_gas_swelling']:11.3e} "
          f"{out['fission_gas_release']:8.4f}")
    T_old = T_new

print("\n[OK] binding drove SCIANTIX end-to-end. Compare against a standalone "
      "input_history.txt run with the same (T, fission rate, dt) to validate.")
