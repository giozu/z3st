#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---

import os
import numpy as np
import matplotlib.pyplot as plt
from z3st.utils.utils_extract_vtu import *
from z3st.utils.utils_plot import *
from z3st.utils.utils_verification import *

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(__file__)
VTU_FILE = os.path.join(CASE_DIR, "output", "fields.vtu")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

# Parametri fisici fissi
Ti, To = 350.0, 350.0
P_applied = 2.0e6  # Pa
TOLERANCE = 1.0e-1 

# --.. ..- .-.. .-.. --- execution --.. ..- .-.. .-.. ---
if not os.path.exists(VTU_FILE):
    print(f"[ERROR] VTU file not found at {VTU_FILE}")
    exit(1)

# Estrazione info mesh automatica
list_fields(VTU_FILE)
# Assumiamo y come altezza (z) e x come raggio (r)
y_target = 0.05 
mask_tol = 0.005

# --- 1. Temperature Extraction ---
x_raw, y_raw, z_raw, T_all = extract_field(VTU_FILE, field_name="Temperature")
mask_T = np.abs(y_raw - y_target) < mask_tol
sort_idx_T = np.argsort(x_raw[mask_T])

x_T = x_raw[mask_T][sort_idx_T]
T_num = T_all[mask_T][sort_idx_T]
T_ref = np.full_like(x_T, Ti)

# --- 2. Stress Extraction & Merge ---
xs_1, ys_1, zs_1, S1 = extract_field(VTU_FILE, field_name="Stress_cyl_inner (cells)")
xs_2, ys_2, zs_2, S2 = extract_field(VTU_FILE, field_name="Stress_cyl_outer (cells)")

# Pulizia: scartiamo i valori nulli (fuori dal dominio del materiale)
mask_act1 = np.any(np.abs(S1) > 1e-3, axis=1)
mask_act2 = np.any(np.abs(S2) > 1e-3, axis=1)

xs_all = np.concatenate([xs_1[mask_act1], xs_2[mask_act2]])
ys_all = np.concatenate([ys_1[mask_act1], ys_2[mask_act2]])
S_all = np.concatenate([S1[mask_act1], S2[mask_act2]])

# Filtro sulla sezione centrale
mask_s = np.abs(ys_all - y_target) < mask_tol
sort_idx_s = np.argsort(xs_all[mask_s])

x_s = xs_all[mask_s][sort_idx_s]
# Indici tensore 9 comp: 0=xx (radiale), 4=yy (assiale)
sigma_xx = S_all[mask_s, 0][sort_idx_s]
sigma_yy = S_all[mask_s, 4][sort_idx_s]
sigma_th_ref = np.full_like(x_s, -P_applied)

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
Pa_to_MPa = 1e-6
plt.figure(figsize=(10, 7))
ax1 = plt.gca()

ax1.plot(x_s, sigma_xx * Pa_to_MPa, "b-o", label=r"Num. $\sigma_{xx}$ (Radial)", markersize=4)
ax1.plot(x_s, sigma_yy * Pa_to_MPa, "r-s", label=r"Num. $\sigma_{yy}$ (Axial)", markersize=4)
ax1.plot(x_s, sigma_th_ref * Pa_to_MPa, "k--", label="Ref. Pressure", alpha=0.8)

ax1.set_xlabel("Radius x (m)")
ax1.set_ylabel("Stress (MPa)")
ax1.grid(True, linestyle="--", alpha=0.5)

ax2 = ax1.twinx()
ax2.plot(x_T, T_num, "gD", label="Temperature", markersize=3, alpha=0.4)
ax2.set_ylabel("Temperature (K)")
ax2.set_ylim([Ti-2, Ti+2])

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc="best")

plt.title(f"Z3ST Verification: Coaxial Contact Profile at y={y_target}m")
plt.tight_layout()
plt.savefig(os.path.join(CASE_DIR, "output", "stress_comparison.png"), dpi=300)

# --.. ..- .-.. .-.. --- metrics --.. ..- .-.. .-.. ---
L2_T = float(np.sqrt(np.mean((T_num - T_ref) ** 2)))
# Errore relativo sullo stress radiale rispetto alla pressione applicata
err_sigma_xx = np.sqrt(np.mean((sigma_xx - sigma_th_ref) ** 2)) / np.abs(np.mean(sigma_th_ref))

errors = {
    "L2_error_T": {"numerical": L2_T, "reference": 0.0, "rel_error": L2_T/Ti},
    "L2_error_sigma_xx": {"numerical": float(err_sigma_xx), "reference": 0.0, "rel_error": float(err_sigma_xx)},
}

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)
print("\n[INFO] Non-regression completed.\n")