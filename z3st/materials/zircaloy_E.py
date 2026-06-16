# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: temperature-dependent Young's modulus for Zircaloy-4 cladding
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Zircaloy-4 Young's modulus as a symbolic E(T) material card.

A card opts in via ``E: materials.zircaloy_E.E``
- spine resolves the function and builds lmbda/G as UFL expressions in the live temperature field (see
spine.initialize_fields), exactly like the symbolic k(T) hook.

This stub returns a *constant* 99.3 GPa, written
as a UFL expression in T so the symbolic-E code path is exercised
while the result is numerically identical to the scalar-E baseline.
A correlation, e.g. the linear-softening form is

    E(T) = 1.088e11 - 5.475e7 * T            (Pa)   (~99 GPa at 300 K)

"""


def E(T):
    """Zircaloy-4 Young's modulus E(T) (Pa) as a UFL expression in T (K)"""
    return 99.3e9 + 0.0 * T
