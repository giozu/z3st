# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: radial and axial power form factors for materials with a heat source
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Radial f(r, bu) and axial f(z) power form factors for materials with a heat source (e.g. fissile materials). 
These are *multiplicative* shaping factors that modulate the linear heat rate distribution within the fuel, redistributing it radially and axially.

A material card may name these via ``radial_profile: materials.fuel_profiles.<name>`` and/or ``axial_profile: materials.fuel_profiles.<name>``;
``spine.set_power`` resolves them, evaluates them on the fuel dofs, multiplies
them together, normalises the composite to mean 1 (so the shaping only
*redistributes* the linear heat rate, never changes its integral) and multiplies
the volumetric source by it.

Signature::

    f(coords, burnup, material, model=None) -> ndarray   # one multiplier per dof

- ``coords``   : (N, 3) dof coordinates.
- ``burnup``   : (N,)  current local burnup [MWd/kgU] — lets f depend on bu.
- ``material`` : the material card dict (read parameters from it).
- ``model``    : the solver/spine (regime, geometry, ...).

These are parametric stand-ins: a mechanistic TUBRNP profile (Pu-239 build-up
from U-238 resonance capture, flux depression) can replace them later behind the
same interface, with no change to set_power or the burnup bus.
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


def _axial_coord(coords, model):
    """Axial coordinate from the dof coordinates, respecting the regime: in
    axisymmetric (r-z) meshes the second coordinate is z; in 3-D it is the
    third."""
    regime = getattr(model, "regime", None)
    if regime == "3d":
        return coords[:, 2]
    return coords[:, 1]


def chopped_cosine(coords, burnup, material, model=None):
    """Chopped-cosine axial form factor (Todreas & Kazimi, 1-D axial problem):

        f(z) = cos(pi (z - z_mid) / L')

    with L the active (meshed) fuel height and L' >= L the *extrapolated*
    length. z_mid and L are inferred from the fuel dof extent; L' comes from
    the material card:

        axial_extrapolated_length  L'  (default 1.1 * L)

    set_power normalises the composite form factor to mean 1, so the prescribed
    LHR stays the segment average and the axial peaking factor emerges as

        max f_norm = 1 / [ (2 L' / pi L) sin(pi L / 2 L') ].

    L' < L would make the cosine negative beyond the extrapolated boundary, so
    the profile is clamped at zero (physically: no fission outside it).
    """
    z = _axial_coord(coords, model)
    z_min, z_max = float(z.min()), float(z.max())
    L = z_max - z_min if z_max > z_min else 1.0
    z_mid = 0.5 * (z_min + z_max)
    L_prime = float(material.get("axial_extrapolated_length", 1.1 * L))
    f = np.cos(np.pi * (z - z_mid) / L_prime)
    return np.clip(f, 0.0, None)


def tabulated_axial(coords, burnup, material, model=None):
    """Tabulated axial form factor — piecewise-linear interpolation of user
    points (e.g. node-wise peaking factors from a core-physics calculation),
    read from the material card:

        axial_table_z: [z1, z2, ...]   # (m) elevations, same axis as the mesh
        axial_table_f: [f1, f2, ...]   # (-) relative power at each elevation

    Outside the table range the end values are held (``np.interp`` clamping).
    set_power normalises the composite form factor to mean 1, so only the
    *shape* of the table matters — its absolute scale is irrelevant (peaking
    factors and raw kW/m readings give the same source).
    """
    z = _axial_coord(coords, model)
    z_tab = np.asarray(material["axial_table_z"], dtype=float)
    f_tab = np.asarray(material["axial_table_f"], dtype=float)
    if len(z_tab) != len(f_tab) or len(z_tab) < 2:
        raise ValueError(
            "axial_table_z and axial_table_f must be equal-length lists with >= 2 points"
        )
    order = np.argsort(z_tab)
    return np.interp(z, z_tab[order], f_tab[order])


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
