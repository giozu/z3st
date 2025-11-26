"""
Post-processing for attenuation map (thermal shield)
→ Reads non-regression.json results, normalizes σθθ, and produces two plots:
   1) σ_T vs μRᵢ (grouped by (Rₒ-Rᵢ)/Rᵢ)
   2) σ_T vs Rₒ/Rᵢ (grouped by μRᵢ)

Author: Giovanni Zullo
"""

import json
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# CONFIG
# =============================================================================
ROOT = Path.cwd()

# =============================================================================
# GATHER RESULTS
# =============================================================================
data_rows = []

for sub in sorted(ROOT.glob("ba_*_mua_*")):
    json_path = sub / "output/non-regression.json"
    mat_path = sub / "vessel_steel.yaml"

    if not json_path.exists():
        print(f"[WARN] Missing: {json_path}")
        continue
    if not mat_path.exists():
        print(f"[WARN] Missing: {mat_path}")
        continue

    try:
        # --- Read from non-regression.json ---
        data = json.loads(json_path.read_text())
        sigma_tt = data["results"]["sigma_th_tt"]["numerical"]

        # --- Read material ---
        mat = yaml.safe_load(open(mat_path))
        E = float(mat["E"])
        nu = float(mat["nu"])
        alpha = float(mat["alpha"])
        k = float(mat["k"])
        q0 = float(mat["gamma_heating"])
        mu = float(mat["mu_gamma"])

        # --- Parameters from folder name ---
        parts = sub.name.split("_")
        ba = float(parts[1]) / 100.0      # Ro/Ri
        mua = float(parts[-1])            # mu*Ri
        slenderness = 0.5 * (1 + ba) / (ba - 1) # 0.5 * (Ro + Ri) / (Ro - Ri)

        # --- Normalization ---
        sigma_T = sigma_tt / (alpha * E * q0 / ((1 - nu) * k * mu**2))

        data_rows.append({
            "Ro/Ri": ba,
            "slenderness": slenderness,
            "muRi": mua,
            "sigma_tt": sigma_tt,
            "sigma_T": sigma_T,
            "case": sub.name
        })
        print(f"→ {sub.name}: σθθ={sigma_tt:.3e} Pa → σ_T={sigma_T:.3f}")

    except Exception as e:
        print(f"[ERROR] Failed to process {sub}: {e}")

# =============================================================================
# CREATE DATAFRAME + CSV
# =============================================================================
if not data_rows:
    print("[WARN] No valid results found.")
    exit()

df = pd.DataFrame(data_rows)
df.sort_values(by=["slenderness", "muRi"], inplace=True)
csv_path = ROOT / "attenuation_map_results.csv"
df.to_csv(csv_path, index=False)
print(f"[OK] Saved CSV → {csv_path}")

# =============================================================================
# PLOT 1 – σ_T vs μRᵢ (grouped by slenderness)
# =============================================================================
plt.figure(figsize=(8, 6))
for rori, group in df.groupby("Ro/Ri"):
    group = group.sort_values("muRi")
    plt.plot(group["muRi"], group["sigma_T"], marker="o", label=f"Ro/Ri={rori:.3f}")

ba_11 = pd.read_csv("b_a_11.txt", header=None)
plt.plot(ba_11.iloc[:, 0], ba_11.iloc[:, 1], "--", label=r"$Ro/Ri = 1.10$ (W)")

ba_12 = pd.read_csv("b_a_12.txt", header=None)
plt.plot(ba_12.iloc[:, 0], ba_12.iloc[:, 1], "--", label=r"$Ro/R_i = 1.20$ (W)")

plt.xlabel(r"$\mu R_i$ (adim.)")
plt.ylabel(r"$\sigma_T = \sigma_{\theta\theta} / [\alpha E q_0 / ((1-\nu)k\mu^2)]$")
plt.title("Normalized thermal stress vs attenuation (grouped by Ro/Ri)")
plt.grid(True)
plt.legend(title="Slenderness")
plt.tight_layout()
plot_path1 = ROOT / "attenuation_map_sigmaT_vs_muRi.png"
plt.savefig(plot_path1, dpi=200)
print(f"[OK] Saved plot → {plot_path1}")
plt.show()

# =============================================================================
# PLOT 2 – σ_T vs Rₒ/Rᵢ (grouped by μRᵢ)
# =============================================================================
plt.figure(figsize=(8, 6))
for mua, group in df.groupby("muRi"):
    group = group.sort_values("Ro/Ri")
    plt.plot(group["Ro/Ri"], group["sigma_T"], marker="s", label=fr"$\mu R_i$ = {mua:.0f}")

mua10 = pd.read_csv("mua_10.txt", header=None)
plt.plot(mua10.iloc[:, 0], mua10.iloc[:, 1], "--", label=r"$\mu R_i = 10$ (GS)")

mua20 = pd.read_csv("mua_20.txt", header=None)
plt.plot(mua20.iloc[:, 0], mua20.iloc[:, 1], "--", label=r"$\mu R_i = 20$ (GS)")

mua50 = pd.read_csv("mua_50.txt", header=None)
plt.plot(mua50.iloc[:, 0], mua50.iloc[:, 1], "--", label=r"$\mu R_i = 50$ (GS)")

plt.xlabel(r"$R_o/R_i$ (adim.)")
plt.ylabel(r"$\sigma_T = \sigma_{\theta\theta} / [\alpha E q_0 / ((1-\nu)k\mu^2)]$")
plt.title("Normalized thermal stress vs Rₒ/Rᵢ (grouped by μRᵢ)")
plt.grid(True)
plt.legend(title=r"$\mu R_i$")
plt.tight_layout()
plot_path2 = ROOT / "attenuation_map_sigmaT_vs_RoRi.png"
plt.savefig(plot_path2, dpi=200)
print(f"[OK] Saved plot → {plot_path2}")
plt.show()