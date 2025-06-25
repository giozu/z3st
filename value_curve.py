import matplotlib.pyplot as plt
import numpy as np

# Features
features = [
    "Language",
    "Thermo-mechanical model",
    "Multi-material support",
    "Arbitrary geometries",
    "Computational time",
    "Integrated post-processing",
    "Flexible input",
    "Modularity/extensibility",
    "Coupling frameworks", 
    "Conventional applications",
    "License"
]

# Scores
scores = {
    "ZEST":         [5, 4, 5, 5, 2, 4, 5, 5, 5, 2, 5],
    "BISON":        [2, 2, 2, 2, 5, 2, 3, 2, 2, 4, 2],
    "TRANSURANUS":  [1, 2, 2, 1, 5, 2, 1, 3, 2, 5, 1],
    # "ABAQUS":       [2, 4, 5, 5, 5, 5, 2, 2, 2, 3, 1],
    "FAST":         [1, 1, 1, 1, 5, 2, 1, 1, 1, 5, 1],
    "OFFBEAT":      [2, 2, 3, 3, 3, 1, 1, 3, 3, 4, 5],
    # "FALCON":       [1, 3, 2, 3, 5, 2, 1, 2, 2, 5, 1]
}

x = np.arange(len(features))
plt.figure(figsize=(16, 9))

# Line styles
for name, vals in scores.items():
    plt.plot(
        x, vals, marker='o', markersize=6, linewidth=2, label=name
    )

# Y-axis
y_labels = {
    1: "Very Basic",
    2: "Basic",
    3: "Good",
    4: "Very Good",
    5: "Excellent"
}
plt.yticks(list(y_labels.keys()), list(y_labels.values()), fontsize=12)

# X-axis
plt.xticks(x, features, rotation=60, ha="right", fontsize=11)

# Labels and grid
plt.ylim(0.8, 5.2)
plt.grid(True, linestyle='--', alpha=0.6)
# plt.title("Value Curve Comparison of ZEST and Other Fuel Performance Codes", fontsize=16)
plt.ylabel("Capability Level", fontsize=13)

# Legend
plt.legend(title="Code", title_fontsize=12, fontsize=11, loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig("output/value_curve_comparison_detailed.png", dpi=300)
plt.show()
