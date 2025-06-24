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
    A = (pi * a**2 - po * b**2) / (b**2 - a**2)
    B = (po - pi) * (a * b)**2 / (b**2 - a**2)

    sigma_r = A + B / r**2
    sigma_theta = A - B / r**2
    sigma_z_plane_strain = 2 * nu * A * np.ones_like(r)
    sigma_z_plane_stress = np.zeros_like(r)

    return sigma_r, sigma_theta, sigma_z_plane_strain, sigma_z_plane_stress


def thermal_stresses(r, alpha=1.45e-5, E=170e9, T_hot=1050.63, T_cold=573, nu=0.3, R_s=0.004):
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
