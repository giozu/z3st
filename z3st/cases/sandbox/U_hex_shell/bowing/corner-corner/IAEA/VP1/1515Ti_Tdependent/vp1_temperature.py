# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 temperature of ../vp1_temperature.py, under the same names, so the
thermal Dirichlet (temperature: vp1_temperature.T_vp1) resolves to it from here.
"""

import importlib.util
import os

_PARENT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "vp1_temperature.py")
_spec = importlib.util.spec_from_file_location("vp1_temperature_parent", _PARENT)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)

globals().update({k: v for k, v in vars(_base).items() if not k.startswith("__")})
