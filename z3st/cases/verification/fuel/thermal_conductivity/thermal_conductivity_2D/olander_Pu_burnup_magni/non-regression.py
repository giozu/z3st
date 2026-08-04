#!/usr/bin/env python3
import os

SHARED = os.path.join(os.path.dirname(__file__), "..", "case_non_regression.py")
with open(SHARED, "r", encoding="utf-8") as f:
    exec(compile(f.read(), __file__, "exec"), {"__file__": __file__, "__name__": "__main__"})
