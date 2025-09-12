import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("output/load_disp.csv")

# Plot vertical reaction force (Ry) vs. avg displacement
plt.figure(figsize=(7,5))
plt.plot(df["avg_disp"], df["Ry"], label="Reaction force (Ry)", lw=2)

# Optionally, also plot the total reaction force magnitude
plt.plot(df["avg_disp"], df["Rmag"], "--", label="|Reaction|", lw=2)

plt.xlabel("Average displacement at loaded face [m]", fontsize=12)
plt.ylabel("Reaction force [N]", fontsize=12)
plt.title("Load–displacement curve", fontsize=14)
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# --- Young modulus ---
# Convert to numeric (ignore errors like headers)
df["Ry"] = pd.to_numeric(df["Ry"], errors="coerce")
df["avg_disp"] = pd.to_numeric(df["avg_disp"], errors="coerce")

# Drop NaN rows if any
df = df.dropna(subset=["Ry", "avg_disp"])

# Geometry
A0 = 0.0004
l0 = 0.1

# Elastic modulus curve
E_curve = (df["Ry"].to_numpy()[1:] / df["avg_disp"].to_numpy()[1:]) * (l0 / A0)

for i, E in enumerate(E_curve, start=1):
    print(f"Step {i}: E = {E:.3e} Pa")

# Approximate E as average slope over first increments
E_est = E_curve[:5].mean()
print(f"Estimated elastic modulus: {E_est:.2e} Pa")