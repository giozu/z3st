"""Magni et al. MA-MOX thermal-conductivity correlation.

The material card opts in with::

    k: materials.magni_mox_thermal.k

Composition parameters are read from the material card when available:
``Pu``, ``Am``, ``Np`` as atomic fractions, ``x`` as deviation from
stoichiometry (or ``OM``/``O_M``), ``p`` as porosity fraction, and ``burnup``
or the model burnup field in GWd/tHM (numerically MWd/kgHM).
"""

import numpy as np

try:
    import ufl
except ModuleNotFoundError:  # Allows offline NumPy diagnostics without FEniCSx.
    ufl = None


# Magni et al. / INSPYRE coefficients. Pu, Am, Np are atomic fractions.
A0 = 0.01926
Ax = 1.06e-6
APu = 2.63e-8
AAm = 0.596
ANp = 2.22e-14
B0 = 2.39e-4
BPu = 1.37e-13
BAm = 5.47e-4
BNp = 2.48e-14
D = 5.27e9
E = 17109.5
K_INF = 1.755
PHI = 128.75


def _card_value(material, *keys, default=0.0):
    if material is None:
        return default
    for key in keys:
        if key in material:
            return float(material[key])
    return default


def _fraction_value(material, *keys, default=0.0):
    value = _card_value(material, *keys, default=default)
    if value is None:
        return value
    # Material cards normally use fractions (0.20). Accept percentages (20)
    # as a convenience because the paper reports applicability in at.%.
    return value / 100.0 if abs(value) > 1.0 else value


def _stoich_deviation(material):
    if material is None:
        return 0.0
    if "x" in material:
        return float(material["x"])
    for key in ("OM", "O_M", "oxygen_to_metal"):
        if key in material:
            return 2.0 - float(material[key])
    return 0.0


def _porosity(material):
    p = _card_value(material, "p", "porosity", default=None)
    if p is not None:
        return p / 100.0 if abs(p) > 1.0 else p
    rho_rel = _card_value(material, "relative_density", "density_fraction", default=None)
    if rho_rel is not None:
        return 1.0 - rho_rel
    return 0.0


def _as_fraction_array(value):
    arr = np.asarray(value, dtype=float)
    return np.where(np.abs(arr) > 1.0, arr / 100.0, arr)


def _burnup(material, model):
    if model is not None and getattr(model, "burnup", None) is not None:
        return model.burnup
    return _card_value(material, "burnup", "bu", default=0.0)


def k(T, material=None, model=None):
    """Return k(T, x, Pu, Am, Np, p, burnup) as a UFL expression."""
    if ufl is None:
        raise ModuleNotFoundError("materials.magni_mox_thermal.k requires ufl")
    pu = _fraction_value(material, "Pu", "pu", "Pu_fraction", "plutonium", default=0.0)
    am = _fraction_value(material, "Am", "am", "Am_fraction", "americium", default=0.0)
    np_ = _fraction_value(material, "Np", "np", "Np_fraction", "neptunium", default=0.0)
    x = _stoich_deviation(material)
    p = _porosity(material)
    bu = _burnup(material, model)

    a = A0 + Ax * x + APu * pu + AAm * am + ANp * np_
    b = B0 + BPu * pu + BAm * am + BNp * np_
    k0 = (1.0 / (a + b * T) + D / (T * T) * ufl.exp(-E / T)) * (1.0 - p) ** 2.5
    return K_INF + (k0 - K_INF) * ufl.exp(-bu / PHI)


def k_numpy(T, Pu=0.0, Am=0.0, Np=0.0, x=0.0, p=0.0, burnup=0.0):
    """Numpy/scalar version used by GPR fitting and diagnostics."""
    T = np.asarray(T, dtype=float)
    Pu = _as_fraction_array(Pu)
    Am = _as_fraction_array(Am)
    Np = _as_fraction_array(Np)
    p = _as_fraction_array(p)
    a = A0 + Ax * x + APu * Pu + AAm * Am + ANp * Np
    b = B0 + BPu * Pu + BAm * Am + BNp * Np
    k0 = (1.0 / (a + b * T) + D / (T * T) * np.exp(-E / T)) * (1.0 - p) ** 2.5
    return K_INF + (k0 - K_INF) * np.exp(-burnup / PHI)


def dk_dT_numpy(T, Pu=0.0, Am=0.0, Np=0.0, x=0.0, p=0.0, burnup=0.0):
    """Derivative of ``k_numpy`` with respect to temperature."""
    T = np.asarray(T, dtype=float)
    Pu = _as_fraction_array(Pu)
    Am = _as_fraction_array(Am)
    Np = _as_fraction_array(Np)
    p = _as_fraction_array(p)
    a = A0 + Ax * x + APu * Pu + AAm * Am + ANp * Np
    b = B0 + BPu * Pu + BAm * Am + BNp * Np
    por = (1.0 - p) ** 2.5
    d_lattice = -b / (a + b * T) ** 2
    d_electronic = D * np.exp(-E / T) * (E / T**4 - 2.0 / T**3)
    return por * (d_lattice + d_electronic) * np.exp(-burnup / PHI)


def applicability():
    return {
        "T_K": (500.0, 2700.0),
        "x": (0.0, 0.05),
        "Pu_atomic_fraction": (0.0, 0.45),
        "Am_atomic_fraction": (0.0, 0.10),
        "Np_atomic_fraction": (0.0, 0.12),
        "porosity_fraction": (0.0, 0.10),
        "burnup_GWd_tHM": (0.0, 130.0),
    }


__all__ = [
    "k",
    "k_numpy",
    "dk_dT_numpy",
    "applicability",
]
