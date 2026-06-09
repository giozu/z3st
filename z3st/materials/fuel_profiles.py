# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: radial power form factors for fissile materials (the source bus)
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Radial power form factors f(r, bu) for fissile materials.

A fissile material card may name one of these via ``radial_profile:
materials.fuel_profiles.<name>``; ``spine.set_power`` resolves it, evaluates it on
the fuel dofs, normalises it to mean 1 (so it only *redistributes* the linear
heat rate radially, never changes its integral) and multiplies the volumetric
source by it.

Signature (the source-bus contract)::

    f(coords, burnup, material, model=None) -> ndarray   # one multiplier per dof

- ``coords``   : (N, 3) dof coordinates.
- ``burnup``   : (N,)  current local burnup [MWd/kgU] — lets f depend on bu.
- ``material`` : the material card dict (read parameters from it).
- ``model``    : the solver/spine (regime, geometry, ...).

These are the *parametric* stand-ins of the "parametric first, TUBRNP behind it"
plan: a mechanistic TUBRNP profile (Pu-239 build-up from U-238 resonance capture,
flux depression) drops in here later behind exactly this interface, with no
change to set_power or the burnup bus.
"""

import numpy as np


def _radius(coords, model):
    """In-plane radius from the dof coordinates, respecting the regime. For
    axisymmetric / 2-D cylindrical meshes the first coordinate *is* r (the mesh
    lives in the r-z plane); in 3-D it is sqrt(x^2 + y^2)."""
    regime = getattr(model, "regime", None)
    if regime == "3d":
        return np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    return coords[:, 0]


def rim_peaking(coords, burnup, material, model=None):
    """Parametric rim-peaking radial form factor — a stand-in for the TUBRNP
    Pu-rim profile.

        f(r) = 1 + A (r / R)^p

    flat through the pellet interior and rising steeply at the surface, where
    resonance capture in U-238 breeds Pu-239 and the local rating peaks. ``R`` is
    the pellet outer radius (taken as max r over the fuel). Parameters come from
    the material card with sensible defaults:

        radial_peak_amplitude  A  (default 3.0) — rim peak height above the core
        radial_peak_exponent   p  (default 8.0) — how tightly peaked at the rim

    set_power normalises f to mean 1, so the area-average rating (and hence total
    power) is unchanged; only its radial distribution is shaped.
    """
    r = _radius(coords, model)
    R = r.max() if r.max() > 0 else 1.0
    A = float(material.get("radial_peak_amplitude", 3.0))
    p = float(material.get("radial_peak_exponent", 8.0))
    return 1.0 + A * (r / R) ** p
