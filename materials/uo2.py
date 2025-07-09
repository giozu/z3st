import ufl

def k(T, bu=0.0, gad=0.0):
    """
    Compute the thermal conductivity of 95% theoretical density (TD) UO₂ and UO₂-Gd₂O₃ fuel
    using the Nuclear Fuel Industries (NFI) model modified to account for burnup and 
    gadolinia content effects.

    This model is described in:
        Ohira and Itagaki, 1997 (NFI)
        Lanning et al., 2005 (modification for burnup and Gd content)
        Reference: PNNL-35702

    The conductivity is given by:

        k_95(T, Bu, gad) = 1 / [A + gad + B*T + f(Bu) + (1 - 0.9 * exp(-0.04 * Bu)) * g(Bu) * h(T)] + C/T**2 * exp(-D/T)

    Where:
        - T      : Temperature (K)
        - Bu     : Burnup (GWd/MTU)
        - gad    : Gadolinium weight fraction (unitless)

        - A = 0.0452           (m·K/W)
        - B = 2.46e-4          (m·K/W/K)
        - C = 3.5e9            (W·K/m)
        - Q = 6380.0           (K)
        - D = 16361.0          (K)

        - f(Bu) = 0.00187 * Bu
        - g(Bu) = 0.038 * Bu^0.28
        - h(T)  = (T^2 / D) * exp(1 / [1 + 396 * exp(Q / T)])

    Returns:
        k (W/m-K): Thermal conductivity as a UFL expression
    """

    # Constants
    A = 0.0452       # m·K/W
    B = 2.46e-4      # m·K/W/K
    C = 3.5e9        # W·K/m
    Q = 6380.0       # K
    D = 16361.0      # K

    # Sub-functions
    f_bu = 0.00187 * bu
    g_bu = 0.038 * bu**0.28
    h_T = (T**2 / D) * ufl.exp(1 / (1 + 396 * ufl.exp(Q / T)))

    # Conductivity formula
    denominator = A + B * T
    return 1.0 / denominator
