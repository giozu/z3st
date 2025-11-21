"""
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
Author: Giovanni Zullo
Version: 0.1.0 (2025)
--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""

import importlib, sys
from pathlib import Path

_pkg_root = Path(__file__).resolve().parent.parent
for _name in ("core", "models", "materials", "utils"):
    sys.path.insert(0, str(_pkg_root / _name))
    globals()[_name] = importlib.import_module(_name)

__all__ = ["core", "models", "materials", "utils"]
__version__ = "0.1.0"
