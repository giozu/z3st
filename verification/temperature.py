import numpy as np
import matplotlib.pyplot as plt

def temperature_analytical(lhr, R_in, R_out, T_outer, k, r):
    """
    Calculates the temperature distribution in a cylinder in steady-state with
    a boundary condition on the outer surface and volumetric heat generation.
    
    Parameters:
    q_third (float): Volumetric heat generation [W/m^3].
    R_in (float): Inner radius of the cylinder [m].
    R_out (float): Outer radius of the cylinder [m].
    T_outer (float): Temperature on the outer surface of the cylinder [K].
    k (float): Thermal conductivity of the material [W/mK].
    r (array): Radial distance from the center of the cylinder.
    
    Returns:
    array: Temperature T(r) as a function of radial distance.
    """
    # Analytical solution for the temperature gradient with volumetric heat generation
    T = T_outer + (lhr / (4 * np.pi * k)) * ((R_out**2 - r**2)/(R_out**2 - R_in**2) + (R_in**2)/(R_out**2 - R_in**2) * 2 * np.log(r/R_out))
    return T

# Parameters of the cylinder
R_in = 0.5  # Inner radius [m]
R_out = 1.50  # Outer radius [m]
T_outer = 700  # Temperature on the outer surface [K]
k = 2.5  # Thermal conductivity [W/mK]

lhr = 15000 # W/m
q_third = lhr / (np.pi * (R_out**2 - R_in**2))

# Radial distance range
r = np.linspace(R_in, R_out, 100)  # Radial distance ranging from R_in to R_out

# Calculate the temperature as a function of r
T = temperature_analytical(lhr, R_in, R_out, T_outer, k, r)

delta_T = (lhr / (4 * np.pi * k)) * (1 + (R_in**2)/(R_out**2 - R_in**2) * 2 * np.log(R_in/R_out))
print(f"LHR = {lhr/1000:.2f} kW/m")
print(f"q_third = {q_third/1000:.2f} kW/m^3")
print(f"Temperature difference between the inner and outer surfaces: {delta_T:.2f} K")

# Display the result
plt.plot(r, T, label="Analytical Temperature", color='b')
plt.xlabel("Radial Distance (m)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Distribution in a Cylinder with Volumetric Heat Generation and Inner Radius")
plt.grid(True)
plt.legend()
plt.show()
