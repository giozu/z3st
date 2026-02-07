import ufl

def k(T):
    """
    Compute the thermal conductivity.

    Returns:
        k: Thermal conductivity as a UFL expression
    """

    return 2.5

def Gc(mesh):
    """
    Compute the fracture toughness.

    Returns:
        Gc (J/m²): Fracture toughness as a UFL expression
    """
    coords = ufl.SpatialCoordinate(mesh)
    y = coords[1] 

    Gc_gb = 2.0     
    Gc_bulk = 100.0 
    
    half_width = 1e-3 

    transition = ufl.tanh(abs(y) / half_width)
    
    return (Gc_gb + (Gc_bulk - Gc_gb) * transition) * 1e-6