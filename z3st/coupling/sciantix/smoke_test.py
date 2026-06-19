#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST <-> SCIANTIX binding — standalone smoke test (PROTOTYPE)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Drive one SCIANTIX point (no Z3ST): fixed fission rate, T ramped over a few
steps, printing swelling / burnup / FGR. Swelling rising with T/burnup means the
ctypes binding and array layout are sound.

Run from a case dir holding input_settings.txt + input_initial_conditions.txt:
    SCIANTIX_LIB=/path/to/libsciantix.so python3 smoke_test.py

ICs must be seeded here: coupling mode skips SCIANTIX's Initialization(), so
without it the intergranular model returns nan and burnup diverges (README §3).
"""

import numpy as np

from sciantix_binding import SciantixSolver

# one point, constant fission rate ~ 1e19 fiss/m^3 s, T ramp 800 -> 1500 K
FISSION_RATE = 1.0e19          # fiss / m^3 s
DT = 1.0e6                     # s per step (~11.6 days)
TEMPS = np.linspace(800.0, 1500.0, 8)   # K

pt = SciantixSolver()
pt.load_initial_conditions("input_initial_conditions.txt")
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
