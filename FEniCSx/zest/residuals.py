import matplotlib.pyplot as plt
import pandas as pd

with open("residuals_u_fuel.tsv", "r") as f:
    lines = f.readlines()

residuals = []
segment_indices = []
segment_id = 0

for line in lines:
    if line.startswith("#"):
        if residuals:
            segment_id += 1
        continue
    parts = line.strip().split()
    if len(parts) == 2:
        iter_idx, res = int(parts[0]), float(parts[1])
        residuals.append((segment_id, iter_idx, res))

df = pd.DataFrame(residuals, columns=["Segment", "Iteration", "Residual"])
df["GlobalIndex"] = range(1, len(df) + 1)

plt.figure(figsize=(10, 6))
plt.plot(df["GlobalIndex"], df["Residual"], label="Residual", color="tab:blue")
plt.axhline(1e-6, color="red", linestyle="--", label="Tolerance (1e-6)")

for seg in df["Segment"].unique()[1:]:
    idx = df[df["Segment"] == seg]["GlobalIndex"].min()
    plt.axvline(idx, color="gray", linestyle="--", linewidth=0.8)

plt.xlabel("Iteration (global index)")
plt.ylabel("Residual")
plt.yscale("log")
plt.title("Mechanical Solver Residuals per Iteration (Fuel)")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.tight_layout()
plt.show()
