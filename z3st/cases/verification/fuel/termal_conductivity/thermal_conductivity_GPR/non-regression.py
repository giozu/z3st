#!/usr/bin/env python3
"""Run the shared Magni/GPR conductivity verification in this case directory."""

import os

CASE_DIR = os.path.dirname(__file__)
SHARED = os.path.join(CASE_DIR, "..", "thermal_conductivity_magni", "non-regression.py")

with open(SHARED, "r", encoding="utf-8") as f:
    code = compile(f.read(), __file__, "exec")

exec(code, {"__file__": __file__, "__name__": "__main__"})
