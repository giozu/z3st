"""Data-driven wrapper for the Magni MA-MOX thermal conductivity law.

This exposes the same ``__call__`` and ``value_and_grad`` contract used by the
NN/GPR thermal-conductivity hooks, so the Newton thermal solver can evaluate the
Magni baseline through ``dolfinx-external-operator``.
"""

import numpy as np

from z3st.materials.magni_mox_thermal import dk_dT_numpy, k_numpy


class MagniConductivity:
    """Magni conductivity for a fixed material composition."""

    def __init__(self, *, Pu=0.0, Am=0.0, Np=0.0, x=0.0, p=0.0, burnup=0.0):
        self.Pu = float(Pu)
        self.Am = float(Am)
        self.Np = float(Np)
        self.x = float(x)
        self.p = float(p)
        self.burnup = float(burnup)

    def __call__(self, T_array):
        T = np.asarray(T_array, dtype=float)
        return k_numpy(
            T,
            Pu=self.Pu,
            Am=self.Am,
            Np=self.Np,
            x=self.x,
            p=self.p,
            burnup=self.burnup,
        ).reshape(T.shape)

    def value_and_grad(self, T_array):
        T = np.asarray(T_array, dtype=float)
        k = self(T)
        dk = dk_dT_numpy(
            T,
            Pu=self.Pu,
            Am=self.Am,
            Np=self.Np,
            x=self.x,
            p=self.p,
            burnup=self.burnup,
        ).reshape(T.shape)
        return k, dk


def _card_value(card, material, *keys, default=0.0):
    for src in (card, material or {}):
        for key in keys:
            if key in src:
                return src[key]
    return default


def load_from_card(card, material=None, base_dir=None):
    """Build a :class:`MagniConductivity` from a material ``k`` card."""
    om = _card_value(card, material, "OM", "O_M", "oxygen_to_metal", default=None)
    x = _card_value(card, material, "x", default=None)
    if x is None:
        x = 0.0 if om is None else 2.0 - float(om)

    return MagniConductivity(
        Pu=_card_value(card, material, "Pu", "pu", "Pu_fraction", "plutonium", default=0.0),
        Am=_card_value(card, material, "Am", "am", "Am_fraction", "americium", default=0.0),
        Np=_card_value(card, material, "Np", "np", "Np_fraction", "neptunium", default=0.0),
        x=x,
        p=_card_value(card, material, "p", "porosity", default=0.0),
        burnup=_card_value(card, material, "burnup", "bu", default=0.0),
    )


__all__ = ["MagniConductivity", "load_from_card"]
