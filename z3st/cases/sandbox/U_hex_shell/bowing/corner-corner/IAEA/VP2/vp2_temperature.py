# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Imposed temperature field of IAEA VP2 (IWGFR/75, App. 2, EXAMPLE 2: single
wrapper, non-linear thermal gradient) and the matching thermal eigenstrain.

As VP1 (../VP1/vp1_temperature.py) except the section: at the mid-wall
corners, 1, 2 and 6 (cold half) are at T1 and 3, 4 and 5 (hot half) at T2,
and T varies linearly between corners along each wall (Fig. 1 of IWGFR/75).
Along the hot-corner direction h (corner 1 at -D/2, corner 4 at +D/2;
h = -y for 'vertical') this is T1 for h <= -D/4, linear for |h| < D/4 (the
two walls parallel to the axis) and T2 for h >= D/4. The axial profile of
T1 and T2 is the one of VP1. s is the coordinate along the flat at
mid-wall, as in hex_post.

T_vp2 is the only place the field is defined; it accepts UFL coordinates and
numpy arrays alike.
"""

import os

import numpy as np
import ufl
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(HERE, "geometry.yaml")) as _f:
    _geo = yaml.safe_load(_f)


D = float(_geo["D"]) - 2.0 * float(_geo["t"]) / np.sqrt(3.0) # (m)
# corner-to-corner depth at mid-wall, the bending direction: 0.150 m in VP1/VP2.
# geometry.yaml D is the outer circumscribed diameter; the mid-wall corner is t/sqrt(3) inside it.

# Hot-corner direction (in rad), as in VP1. 'vertical': corners on +-y as in the MOOSE VP1 input (temp_func, hot at y = -0.075 m);
# 'flat': corners on +-x, hot at -x as in the VTB figure.
PHI_HOT = -np.pi / 2 if _geo.get("orientation", "flat") == "vertical" else np.pi

# VP2 temperatures (as VP1, IWGFR/75 Fig. 1), converted to K
Z_CORE_BOT = 1.5   # (m) core bottom
Z_CORE_TOP = 2.5   # (m) core top
T_INLET = 400.0 + 273.15   # (K) whole section below the core
T1_TOP = 500.0 + 273.15    # (K) cold corners 1, 2, 6, from the core top up
T2_TOP = 550.0 + 273.15    # (K) hot corners 3, 4, 5, from the core top up


def _clip01(v):
    """min(max(v, 0), 1); valid for UFL and numpy."""
    if isinstance(v, ufl.core.expr.Expr):
        return ufl.max_value(0.0, ufl.min_value(1.0, v))
    return np.clip(v, 0.0, 1.0)


def _ramp(z):
    """0 below the core, linear across it, 1 above."""
    return _clip01((z - Z_CORE_BOT) / (Z_CORE_TOP - Z_CORE_BOT))


def T1(z):
    """Cold-corner temperature (K), corners 1, 2, 6."""
    return T_INLET + (T1_TOP - T_INLET) * _ramp(z)


def T2(z):
    """Hot-corner temperature (K), corners 3, 4, 5."""
    return T_INLET + (T2_TOP - T_INLET) * _ramp(z)


def section_shape(h):
    """0 on the cold half (h <= -D/4), 1 on the hot half (h >= D/4), linear between."""
    return _clip01(2.0 * h / D + 0.5)


def T_vp2(x):
    """Temperature (K) at x = (x, y, z): T1 + (T2 - T1) section_shape(h), with h
    the coordinate along the hot-corner direction.
    It supports as input both ufl.SpatialCoordinate and numpy.ndarray, so it is called by both the solver and the scripts."""
    h = x[0] * np.cos(PHI_HOT) + x[1] * np.sin(PHI_HOT)  # -D/2 at the cold corner 1, +D/2 at the hot corner 4
    return T1(x[2]) + (T2(x[2]) - T1(x[2])) * section_shape(h)


# Imposed as a material eigenstrain, as in VP1 (no line/vertex BC in z3st).
def eigenstrain(T, material, model=None, dim=3):
    """alpha * (T_vp2(x) - T_ref) * I, with T_ref the stress-free temperature of the card."""
    x = ufl.SpatialCoordinate(model.mesh)
    return material["alpha"] * (T_vp2(x) - material["T_ref"]) * ufl.Identity(dim)
