#!/usr/bin/env python3
"""IAEA VP1, T-dependent 15-15Ti, Di Marcello set: the non-regression of ../
(beam with eta(T) and E(T) of the card) run in this case directory."""

import os

CASE_DIR = os.path.dirname(__file__)
SHARED = os.path.join(CASE_DIR, "..", "non-regression.py")

with open(SHARED, "r", encoding="utf-8") as f:
    code = compile(f.read(), __file__, "exec")

exec(code, {"__file__": __file__, "__name__": "__main__"})
