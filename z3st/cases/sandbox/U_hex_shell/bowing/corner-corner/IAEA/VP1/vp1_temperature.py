# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Imposed temperature field of IAEA VP1 and the matching thermal eigenstrain.

MOOSE VP1 reads the temperature from an input variable ``temp`` and does not
solve for it; here the thermal model is off and the field enters the
mechanical solve as a material eigenstrain (``eigenstrain:
vp1_temperature.eigenstrain`` in the material card).

T varies linearly across the section, from the cold corner (T1) to the
opposite hot corner (T2), along the direction PHI_HOT (coordinate h, h = -y
for 'vertical'). T_vp1 is the only place the field is defined; it accepts UFL
coordinates and numpy arrays alike. s is the coordinate along the flat at
mid-wall, as in hex_post.
"""

import os

import numpy as np
import ufl
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "geometry.yaml")) as _f:
    _geo = yaml.safe_load(_f)


D = float(_geo["D"]) - 2.0 * float(_geo["t"]) / np.sqrt(3.0) # (m)
# corner-to-corner depth at mid-wall, the bending direction: 0.150 m in VP1.
# geometry.yaml D is the outer circumscribed diameter; the mid-wall corner is t/sqrt(3) inside it.

# Hot-corner direction (in rad). 'vertical': corners on +-y as in the MOOSE VP1 input (temp_func, hot at y = -0.075 m); 
# 'flat': corners on +-x, hot at -x as in the VTB figure.
PHI_HOT = -np.pi / 2 if _geo.get("orientation", "flat") == "vertical" else np.pi

# VP1 temperatures (VTB hex_duct_linear, eqs. T1 and T2), converted to K
Z_CORE_BOT = 1.5   # (m) core bottom
Z_CORE_TOP = 2.5   # (m) core top
T_INLET = 400.0 + 273.15   # (K) whole section below the core
T1_TOP = 500.0 + 273.15    # (K) cold corner, from the core top up
T2_TOP = 550.0 + 273.15    # (K) hot corner, from the core top up

# Step 1: the core-top section over the whole height (f = 1, radial dT only, no z ramp).
# False restores the VP1 axial ramp.
Z_UNIFORM = False


def _ramp(z):
    """0 below the core, linear across it, 1 above; valid for UFL and numpy.
    Identically 1 when Z_UNIFORM is set."""
    if Z_UNIFORM:
        return 1.0 + 0.0 * z
    f = (z - Z_CORE_BOT) / (Z_CORE_TOP - Z_CORE_BOT)
    if isinstance(f, ufl.core.expr.Expr):
        return ufl.max_value(0.0, ufl.min_value(1.0, f))
    return np.clip(f, 0.0, 1.0)


def T1(z):
    """Cold-corner temperature (K)."""
    return T_INLET + (T1_TOP - T_INLET) * _ramp(z)


def T2(z):
    """Hot-corner temperature (K)."""
    return T_INLET + (T2_TOP - T_INLET) * _ramp(z)


def T_vp1(x):
    """Temperature (K) at x = (x, y, z); linear across the section from the cold corner T1 to the hot corner T2.
    It supports as input both ufl.SpatialCoordinate and numpy.ndarray, so it is called by both the solver and the notebook."""
    h = x[0] * np.cos(PHI_HOT) + x[1] * np.sin(PHI_HOT)  # -D/2 at the cold corner, +D/2 at the hot one
    return T1(x[2]) + (T2(x[2]) - T1(x[2])) * (h / D + 0.5)

# This approach is necessary because a BC on a line/ vertix is not supported in z3st, only on surfaces in 3D geometries.
# The temperature is imposed as a material eigenstrain instead of a Dirichlet BC, if it was imposed on the vertexes, the solver would solve it along the perimeter instead of a linear field.
def eigenstrain(T, material, model=None, dim=3):
    """alpha * (T_vp1(x) - T_ref) * I, with T_ref the stress-free temperature of the card."""
    x = ufl.SpatialCoordinate(model.mesh)
    return material["alpha"] * (T_vp1(x) - material["T_ref"]) * ufl.Identity(dim)
