import numpy as np
import matplotlib.pyplot as plt

# Define problem parameters
a = 0.05  # Inner radius (m)
b = 0.15  # Outer radius (m)
pi = 1e6  # Inner pressure (Pa)
po = 15e6  # Outer pressure (Pa)

nu = 0.3

# Radial positions for evaluation
r = np.linspace(a, b, 500)

# Compute constants for the stress distribution
A = (pi * a**2 - po * b**2) / (b**2 - a**2)
B = (a**2 * b**2 * (po - pi)) / (b**2 - a**2)

# Radial stress
sigma_r = A + B / r**2

# Hoop (tangential) stress
sigma_theta = A - B / r**2

# Axial (longitudinal) stress
sigma_z = 2 * nu * A

# Plot results
plt.figure(figsize=(10, 6))

# Normalized radial stress plot
plt.plot(r, sigma_r, label='Radial stress $\sigma_r$', color='blue')

# Normalized hoop stress plot
plt.plot(r, sigma_theta, label='Hoop stress $\sigma_t$', color='red')

# Normalized axial stress plot
plt.axhline(sigma_z, label='Axial stress $\sigma_z$', color='green', linestyle='--')

# Add labels and legend
plt.title('Stress distribution in a thick-walled cylinder', fontsize=14)
plt.xlabel('Radial position (m)', fontsize=12)
plt.ylabel('Stress', fontsize=12)
plt.grid(True)
plt.legend(fontsize=12)

# Show plot
plt.show()
