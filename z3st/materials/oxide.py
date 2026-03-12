import ufl

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
        Gc: Fracture toughness as a UFL expression
    """
    coords = ufl.SpatialCoordinate(mesh)
    y = coords[1] 

    Gc_gb = 0.1     # pJ/micron²
    Gc_bulk = 100.0 # pJ/micron²
    
    half_width = 10e-3 # micron

    transition = ufl.tanh(abs(y) / half_width)
    
    return (Gc_gb + (Gc_bulk - Gc_gb) * transition)