import numpy as np
import matplotlib.pyplot as plt

# Material and problem parameters
R_inner = 0.05  # Inner radius [m]
R_outer = 0.1   # Outer radius [m]
k = 200          # Thermal conductivity [W/mK]
q_prime = 1e6    # Heat generation rate [W/m^3]
alpha = 1e-5     # Thermal expansion coefficient [1/K]
E = 200e9        # Young's modulus [Pa]
nu = 0.3         # Poisson's ratio
T_ref = 300      # Reference temperature [K]

# Radial positions
r = np.linspace(R_inner, R_outer, 100)

# Analytical temperature distribution (parabolic profile)
T = T_ref + (q_prime / (4 * k)) * (R_outer**2 - r**2)

# Lame parameters
lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))
G = E / (2 * (1 + nu))

# Thermal stress sigma_theta (hoop stress)
sigma_theta = (3 * lmbda + 2 * G) * alpha * ((T[0] + T[-1])/2 - T)

# Plot results
plt.figure(figsize=(8, 5))
plt.subplot(2, 1, 1)
plt.plot(r, T, label='Temperature [K]', color='r')
plt.xlabel('Radius [m]')
plt.ylabel('Temperature [K]')
plt.legend()
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(r, sigma_theta / 1e6, label='Hoop stress $\sigma_\theta$ [MPa]', color='b')
plt.xlabel('Radius [m]')
plt.ylabel('Stress [MPa]')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
