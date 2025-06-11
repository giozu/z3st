import numpy as np
import matplotlib.pyplot as plt

# Correcting the error: B1 must be assigned before use

# New data provided
L1 = 0.2
L2 = 0.21
L3 = 0.41
k1 = 2.5
k2 = 28.7
h = 1
T0 = 350
T3 = 290

# B1 = T0
B1 = T0

# Proportionality coefficient between the slopes
r = k1 / k2

# Equation to determine A1
# k1*A1 = h * [A1*(L1 - r*L2 + r*L3) + B1 - T3]
# Correcting the sign of the denominator in the formula

# Correcting the sign of the numerator to reflect the heat flow from hot (T0) to cold (T3)
numerator = h * (T3 - B1)  # Negative, so A1 becomes negative as it should be
denominator = k1 + h * (L1 - r * L2 + r * L3)
A1 = numerator / denominator

# Recalculate constants
A2 = r * A1
B2 = T3 - A2 * L3

# Generate points
x1 = np.linspace(0, L1, 100)
x2 = np.linspace(L2, L3, 100)
T1 = A1 * x1 + B1
T2 = A2 * x2 + B2

print(f"heat flux = {k1 * (T0 - (A1 * L1 + B1)) / L1} W/m2")

# Plotting
plt.figure()
plt.plot(x1, T1, label="Region 1", linewidth=2)
plt.plot(x2, T2, label="Region 2", linewidth=2)
plt.axvline(x=L1, linestyle='--', color='gray', label='Interface')
plt.axvline(x=L2, linestyle='--', color='gray')
plt.xlabel("x [m]")
plt.ylabel("Temperature [K]")
plt.title("1D Steady-State Conduction with Interfacial Resistance (Corrected)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()