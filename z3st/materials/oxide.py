import ufl

# Grain-boundary-weakened fracture toughness (pJ/micron²). Shared by the UFL
# assembly path `Gc(mesh)` and the numpy post-processing path `Gc_numpy(y)`, so
# a verification script cannot drift from the material model (e.g. an accidental
# unit factor).
GC_GB = 0.1          # toughness at the grain boundary (y = 0)
GC_BULK = 100.0      # bulk toughness
GC_HALF_WIDTH = 10e-3  # micron — GB -> bulk transition width


def _gc_profile(y, tanh):
    """Gc(|y|) with a tanh GB->bulk transition. ``tanh`` is ``ufl.tanh`` or
    ``np.tanh`` so the same formula serves UFL and numpy callers."""
    return GC_GB + (GC_BULK - GC_GB) * tanh(abs(y) / GC_HALF_WIDTH)


def k(T):
    """
    Compute the thermal conductivity.

    Returns:
        k: Thermal conductivity as a UFL expression
    """

    return 2.5e+6 # $pW/(micron-K)$.

def Gc(mesh):
    """
    Compute the fracture toughness.

    Returns:
        Gc: Fracture toughness as a UFL expression (pJ/micron²)
    """
    y = ufl.SpatialCoordinate(mesh)[1]
    return _gc_profile(y, ufl.tanh)


def Gc_numpy(y):
    """The same toughness profile evaluated on numpy y-coordinates
    (pJ/micron²), for verification / post-processing."""
    import numpy as np
    return _gc_profile(np.asarray(y, dtype=float), np.tanh)