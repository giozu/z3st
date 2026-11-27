import re
import numpy as np
import matplotlib.pyplot as plt

# 1. Legge il log di Z3ST ed estrae la deformazione equivalente di creep ad ogni step
eps_creep = []
with open("log_z3st_H500.md", "r") as f: # Qui modifico il nome del file log a seconda del contenuto di idrogeno (H = 0 oppure H = 500)
    for line in f:
        m = re.search(r"max equivalent creep strain = ([0-9.eE+-]+)", line)
        if m:
            eps_creep.append(float(m.group(1)) * 100.0)  # conversione in percentuale (%)

eps_creep = np.array(eps_creep)
n = len(eps_creep)

# Ricostruisce tempo e temperatura (come in Fig. 8 del report D2.1: t da 40s a 127.5s, T da 823K a 1173K)
t = np.linspace(40.0, 40.0 + 0.5 * (n - 1), n)
T = 823.0 + 4.0 * (t - 40.0)

# 2. Crea il grafico doppio asse identico alla Figura 8 del report D2.1
fig, ax1 = plt.subplots(figsize=(8, 5))

color1 = "tab:green"
ax1.set_xlabel("t (s)", fontsize=12)
ax1.set_ylabel(r"$\varepsilon_{eff}$ (%)", color=color1, fontsize=12)
ax1.plot(t, eps_creep, color=color1, linewidth=2.5, label="Z3ST (H = 0 wppm)")
ax1.tick_params(axis="y", labelcolor=color1)
ax1.set_ylim(0, 100)
ax1.grid(True, linestyle="--", alpha=0.5)

ax2 = ax1.twinx()
color2 = "tab:red"
ax2.set_ylabel("T (K)", color=color2, fontsize=12)
ax2.plot(t, T, color=color2, linestyle="--", linewidth=2, label="Temperatura")
ax2.tick_params(axis="y", labelcolor=color2)
ax2.set_ylim(823, 1173)

plt.title("Transizione di fase $\\alpha \\to \\beta$ dello Zircaloy ($\\sigma_{eff} = 30$ MPa, $\\dot{T} = 4$ K/s)")
fig.tight_layout()
plt.savefig("figura8_z3st.png", dpi=200)
print(f"Grafico salvato in figura8_z3st.png (punti letti: {n})")