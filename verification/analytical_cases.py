import numpy as np

def lame_solutions(r, a=2.5, b=3.0, pi=1e6, po=1e7, nu=0.3):
    """
    Analytical Lamé solution for a thick-walled cylindrical shell under internal and external pressure.
    The formulas are preserved from the original implementation.

    Args:
        r (array-like): Radial positions (m).
        a (float): Inner radius (m).
        b (float): Outer radius (m).
        pi (float): Inner pressure (Pa).
        po (float): Outer pressure (Pa).
        nu (float): Poisson's ratio.

    Returns:
        tuple: (sigma_r, sigma_theta, sigma_z_plane_strain, sigma_z_plane_stress).
    """
    r = np.asarray(r)

    # Original formula for the constant K
    K = (-po*b**2 + pi*a**2) / (b**2 - a**2)

    # Original formulas for stresses
    sigma_r = -pi * (a/r)**2 + K * (1 - (a/r)**2)
    sigma_theta = +pi * (a/r)**2 + K * (1 + (a/r)**2)
    sigma_z_plane_strain = 2 * nu * K * np.ones_like(r)
    sigma_z_plane_stress = np.zeros_like(r)

    return sigma_r, sigma_theta, sigma_z_plane_strain, sigma_z_plane_stress

def radial_displacement(r, a=2.5, b=3.0, pi=1e6, po=1e7, nu=0.3, E=2e11):
    """
    Calculates the radial displacement u_r for a thick-walled cylinder.
    The formula is preserved from the original implementation.

    Args:
        r (array-like): Radial positions (m).
        a (float): Inner radius (m).
        b (float): Outer radius (m).
        pi (float): Inner pressure (Pa).
        po (float): Outer pressure (Pa).
        nu (float): Poisson's ratio.
        E (float): Young's modulus (Pa).

    Returns:
        np.ndarray: Radial displacement u_r (m).
    """
    r = np.asarray(r)
    
    # Original formula for the constant K
    K = (-po*b**2 + pi*a**2) / (b**2 - a**2)
    
    # Original formula for displacement
    u_r = (1 + nu) / E * ((a/r)**2 * (pi + K) + K * (1 - 2*nu)) * r

    return u_r

def radial_strain(r, a=2.5, b=3.0, pi=1e6, po=1e7, nu=0.3, E=2e11):
    """
    Calculates the radial strain ε_rr (du_r/dr) for a thick-walled cylinder.
    The formula is preserved from the original implementation.

    Args:
        r (array-like): Radial positions (m).
        a (float): Inner radius (m).
        b (float): Outer radius (m).
        pi (float): Inner pressure (Pa).
        po (float): Outer pressure (Pa).
        nu (float): Poisson's ratio.
        E (float): Young's modulus (Pa).

    Returns:
        np.ndarray: Radial strain ε_rr.
    """
    r = np.asarray(r)

    # Original formula for the constant K
    K = (-po*b**2 + pi*a**2) / (b**2 - a**2)

    # Original formula for strain
    strain_r = (1 + nu) / E * (-(a/r)**2 * (pi + K) + K * (1 - 2*nu))

    return strain_r

def thermal_stresses(r, alpha=1.45e-5, E=170e9, T_hot=891.45, T_cold=573, nu=0.3, R_s=0.004):
    """
    Analytical solution for thermal stress in a cylinder with radial temperature gradient.
    The formulas are preserved from the original implementation.

    Args:
        r (array-like): Radial positions (m).
        alpha (float): Coefficient of thermal expansion (1/K).
        E (float): Young's modulus (Pa).
        T_hot (float): Hot side temperature (K).
        T_cold (float): Cold side temperature (K).
        nu (float): Poisson's ratio.
        R_s (float): Outer radius (m).

    Returns:
        tuple: (sigma_r, sigma_theta, sigma_z).
    """
    r = np.asarray(r)
    
    delta_T = T_hot - T_cold
    thermal_stress_term = alpha * E * delta_T / (4 * (1 - nu))

    # Original formulas for thermal stresses
    sigma_r = thermal_stress_term * (r**2 / R_s**2 - 1)
    sigma_theta = thermal_stress_term * (3 * r**2 / R_s**2 - 1)
    sigma_z = sigma_r + sigma_theta

    return sigma_r, sigma_theta, sigma_z

def thermal_shield_slab_analytical_temperature(x, Ti=494, T0=491, q0=3.46e6,
                                               mu=24.0, K=48.1, L=0.045):
    """
    Analytical temperature profile in a slab with exponential volumetric heat generation.
    The formula is preserved from the original implementation.

    Args:
        x (array-like): Position along the slab thickness (m).
        Ti (float): Inner surface temperature (K).
        T0 (float): Outer surface temperature (K).
        q0 (float): Volumetric heat generation (W/m^3).
        mu (float): Attenuation coefficient (1/m).
        K (float): Thermal conductivity (W/m·K).
        L (float): Slab thickness (m).

    Returns:
        np.ndarray: Temperature distribution T(x) (K).
    """
    x = np.asarray(x)

    # Original formula
    term1 = Ti + (T0 - Ti) * x / L
    term2 = (q0 / (mu**2 * K)) * (
        (x / L) * (np.exp(-mu * L) - 1) - (np.exp(-mu * x) - 1)
    )

    return term1 + term2

def full_cylinder_analytical_temperature(r, T0=573, K=2.5, R=0.004, LHR=1e4):
    """
    Analytical temperature profile in a cylindrical element with volumetric heat generation.
    The formula is preserved from the original implementation. Renamed 'x' to 'r'.

    Args:
        r (array-like): Radial position (m).
        T0 (float): Outer surface temperature (K).
        K (float): Thermal conductivity (W/m·K).
        R (float): Outer radius (m).
        LHR (float): Linear Heat Rate (W/m).

    Returns:
        np.ndarray: Temperature distribution T(r) (K).
    """
    r = np.asarray(r)

    # Original formula, with 'x' renamed to 'r' for clarity
    return LHR / (4 * np.pi * K) * (-r**2 / R**2 + 1) + T0