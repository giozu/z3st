import numpy as np
import matplotlib.pyplot as plt

def sigma_rr(r, rho, b, omega, a, nu):
    return ((3 + nu) / 8) * rho * b**2 * omega**2 * \
           (1 + (a**2 - r**2) / b**2 - (a**2 / r**2))

def sigma_theta_theta(r, rho, b, omega, a, nu):
    return ((3 + nu) / 8) * rho * b**2 * omega**2 * \
           (1 + ((3 + nu) * a**2 - (1 + 3 * nu) * r**2) / ((3 + nu) * b**2) + (a**2 / r**2))

# Define parameters
rho = 7800      # density (kg/m³), e.g., steel
omega = 1e3     # angular velocity (rad/s)
a = 0.005       # inner radius
b = 0.004       # outer radius
nu = 0.3        # poisson's ratio

# Define range of r values
r = np.linspace(a, b, 400)

# Compute stress values
sigma_rr_values = sigma_rr(r, rho, b, omega, a, nu)
sigma_theta_theta_values = sigma_theta_theta(r, rho, b, omega, a, nu)

# Plot results
plt.figure(figsize=(8, 6))
plt.plot(r*1e3, sigma_rr_values*1e-6, label="$\\sigma_{r}$", color='b')
plt.plot(r*1e3, sigma_theta_theta_values*1e-6, label="$\\sigma_{\\theta}$", color='r')
plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
plt.xlabel("radial position r (mm)")
plt.ylabel("stress (MPa)")
plt.legend()
plt.grid()
plt.show()

# Compute sigma_theta_max
sigma_theta_max = rho * b**2 * omega**2
print(f"sigma_theta_max: {sigma_theta_max*1e-6} MPa")