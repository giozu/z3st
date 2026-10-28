# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
IAEA VP1 step 1, no axial gradient: ../vp1_temperature.py with Z_UNIFORM set,
i.e. the core-top section (T1 = 500 C, T2 = 550 C) over the whole height.

Same names as the parent module, so the material card (eigenstrain:
vp1_temperature.eigenstrain), hex_post (temperature, closed form vp1:) and
compare_ansys resolve to this file when run from this folder.
"""

import importlib.util
import os

_PARENT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vp1_temperature.py")
_spec = importlib.util.spec_from_file_location("vp1_temperature_parent", _PARENT)
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)
_base.Z_UNIFORM = True   # read at call time by _ramp in the parent module

globals().update({k: v for k, v in vars(_base).items() if not k.startswith("__")})
