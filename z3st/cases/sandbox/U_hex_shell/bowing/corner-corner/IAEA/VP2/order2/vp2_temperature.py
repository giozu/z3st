# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP2 temperature for the Q2 run: ../vp2_temperature.py unchanged.

Same name as the parent module, so the material card (eigenstrain:
vp2_temperature.eigenstrain) and hex_post (temperature) resolve to this
file when run from this folder.
"""

import importlib.util
import os

_PARENT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vp2_temperature.py")
_spec = importlib.util.spec_from_file_location("vp2_temperature_parent", _PARENT)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)

globals().update({k: v for k, v in vars(_base).items() if not k.startswith("__")})
