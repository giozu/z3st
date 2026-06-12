# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: temperature-dependent thermal conductivity for UO2 fuel
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
UO2 thermal conductivity as a symbolic k(T) material card.

A fuel card opts in via ``k: materials.fuel_thermal.k`` — spine resolves the
function and substitutes k(T) into the thermal weak form (the conductivity is
lagged at the previous staggered iterate, so the linear thermal solve becomes
a Picard iteration that converges with the staggered loop).

The correlation is the phonon + electronic form of Fink (J. Nucl. Mater. 279,
2000) for 95 % TD UO2,

    k(T) = 1 / (0.0452 + 2.46e-4 T) + 3.5e9 / T^2 * exp(-16361 / T)   [W/(m.K)]

valid 298-3120 K: ~8.7 at 300 K, ~5.3 at 580 K, ~3.0 at 1200 K, ~2.4 at
1600 K. Burnup degradation (Lucuta factors) is a later refinement behind the
same hook once the k(T) contract carries the state bus.
"""

import ufl


def k(T):
    """Fink (2000) UO2 conductivity, 95 % TD, as a UFL expression in T (K)."""
    lattice = 1.0 / (0.0452 + 2.46e-4 * T)
    electronic = 3.5e9 / (T * T) * ufl.exp(-16361.0 / T)
    return lattice + electronic
