import numpy as np

def lame_solutions(r, a=2.5, b=3.0, pi=1e6, po=1e7, nu=0.3):
    """
    Analytical Lamé solution for a thick-walled cylindrical shell under internal and external pressure.

    Parameters:
        r (array-like): Radial positions (m)
        a (float): Inner radius (m)
        b (float): Outer radius (m)
        pi (float): Inner pressure (Pa)
        po (float): Outer pressure (Pa)
        nu (float): Poisson's ratio

    Returns:
        tuple: (sigma_r, sigma_theta, sigma_z_plane_strain, sigma_z_plane_stress)
    """

    r = np.asarray(r)

    K = (-po*b**2 + pi*a**2) / (b**2 - a**2)

    sigma_r = -pi * (a/r)**2 + K * (1 - (a/r)**2) 
    sigma_theta = +pi * (a/r)**2 + K * (1 + (a/r)**2) 
    sigma_z_plane_strain = 2 * nu * K * np.ones_like(r)
    sigma_z_plane_stress = np.zeros_like(r)

    return sigma_r, sigma_theta, sigma_z_plane_strain, sigma_z_plane_stress

def radial_displacement(r, a=2.5, b=3.0, pi=1e6, po=1e7, nu=0.3, E=2e11):

    r = np.asarray(r)
    K = (-po*b**2 + pi*a**2) / (b**2 - a**2)
    u_r = (1 + nu) / E * ((a/r)**2 * (pi + K) + K * (1- 2*nu)) * r

    return u_r

def radial_strain(r, a=2.5, b=3.0, pi=1e6, po=1e7, nu=0.3, E=2e11):

    r = np.asarray(r)
    K = (-po*b**2 + pi*a**2) / (b**2 - a**2)
    strain_r = (1 + nu) / E * (- (a/r)**2 * (pi + K) + K * (1- 2*nu))

    return strain_r

def thermal_stresses(r, alpha=1.45e-5, E=170e9, T_hot=891.45, T_cold=573, nu=0.3, R_s=0.004):
    """
    Analytical solution for thermal stress in a cylinder with radial temperature gradient.

    Parameters:
        r (array-like): Radial positions (m)
        alpha (float): Coefficient of thermal expansion (1/K)
        E (float): Young's modulus (Pa)
        T_hot (float): Hot side temperature (K)
        T_cold (float): Cold side temperature (K)
        nu (float): Poisson's ratio
        R_s (float): Outer radius (m)

    Returns:
        tuple: (sigma_r, sigma_theta, sigma_z)
    """
    r = np.asarray(r)
    delta_T = T_hot - T_cold
    thermal_stress_term = alpha * E * delta_T / (4 * (1 - nu))

    sigma_r = thermal_stress_term * (r**2 / R_s**2 - 1)
    sigma_theta = thermal_stress_term * (3 * r**2 / R_s**2 - 1)
    sigma_z = sigma_r + sigma_theta

    print(f"[INFO] Thermal stress term: {thermal_stress_term*1e-6:.2f} MPa")

    return sigma_r, sigma_theta, sigma_z

def thermal_shield_slab_analytical_temperature(x):
    """
    Analytical temperature profile in a slab with volumetric heat generation
    and exponential attenuation (e.g., thermal shield or gamma heating).

    Parameters:
        x (array-like): Position along the slab thickness (m)

    Returns:
        np.ndarray: Temperature distribution T(x) (K)
    """
    Ti = 494         # Inner surface temperature (K)
    T0 = 491         # Outer surface temperature (K)
    q0 = 3.46e6      # Volumetric heat generation (W/m^3)
    mu = 24.0        # Attenuation coefficient (1/m)
    K = 48.1         # Thermal conductivity (W/m·K)
    L = 0.045        # Slab thickness (m)

    x = np.asarray(x)

    term1 = Ti + (T0 - Ti) * x / L
    term2 = (q0 / (mu**2 * K)) * (
        (x / L) * (np.exp(-mu * L) - 1) - (np.exp(-mu * x) - 1)
    )

    return term1 + term2

def full_cylinder_analytical_temperature(x):
    """
    Analytical temperature profile in a cylindrical element with volumetric heat generation

    Parameters:
        x (array-like): Position along the slab thickness (m)

    Returns:
        np.ndarray: Temperature distribution T(x) (K)
    """
    T0 = 573         # Outer surface temperature (K)
    K = 2.5          # Thermal conductivity (W/m·K)
    R = 0.004        # Slab thickness (m)
    LHR = 1e4

    x = np.asarray(x)

    return LHR / (4 * np.pi * K) * (- x**2 / R**2 + 1) + T0
