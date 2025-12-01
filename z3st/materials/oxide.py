import ufl

def k(T):
    """
    Compute the thermal conductivity.

    Returns:
        k (W/m-K): Thermal conductivity as a UFL expression
    """
    # return 1.0 * ufl.exp(-0.01 * T)
    # return 10 + 0.02*T + 5e-4*T**2

    return 2.5 / (T*1e-6 + 1.0)