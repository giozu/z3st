#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/plasticity/crystal_fcc_001

non-regression script
---------------------
Multi-slip crystal plasticity, uniaxial tension along a cube axis.
The stress-strain curve and the cumulative slip are compared with the
material-point integration of the same scheme (numpy).

"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np

from z3st.materials.crystal_plasticity import saturation_stress, uniaxial_reference
from z3st.utils.non_regression import finish, load_case, metric, tracked
from z3st.utils.utils_extract_vtu import extract_field

# --.. ..- .-.. .-.. --- configuration --.. ..- .-.. .-.. ---
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUTPUT_DIR, "non-regression.json")
TOLERANCE = 1e-5

geom, inp, mat = load_case(CASE_DIR)
dt = float(inp["time"][-1]) / (int(inp["n_steps"]) - 1)

# --.. ..- .-.. .-.. --- results --.. ..- .-.. .-.. ---
eps, sig, slip = [], [], []
for vtu in sorted(glob.glob(os.path.join(OUTPUT_DIR, "fields_*.vtu"))):
    eps.append(np.mean(extract_field(vtu, "Strain (cells)")[3][:, 8]))
    stress = extract_field(vtu, "Stress (cells)")[3]
    sig.append(np.mean(stress[:, 8]))
    slip.append(np.mean(extract_field(vtu, "CumulativePlasticStrain")[3]))
eps, sig, slip = np.array(eps), np.array(sig), np.array(slip)

# --.. ..- .-.. .-.. --- reference --.. ..- .-.. .-.. ---
sig_ref, slip_ref = uniaxial_reference(mat, eps, dt)
sig_sat = saturation_stress(mat, eps[-1] / float(inp["time"][-1]))

curve_err = np.max(np.abs(sig - sig_ref)) / np.max(np.abs(sig_ref))
lateral = np.max(np.abs(np.delete(stress, 8, axis=1))) / sig[-1]

print(f"[RESULTS] sigma_zz final: {sig[-1]/1e6:.3f} MPa (reference {sig_ref[-1]/1e6:.3f} MPa)")
print(f"          exact saturation stress: {sig_sat/1e6:.3f} MPa")
print(f"          cumulative slip final: {slip[-1]:.6e} (reference {slip_ref[-1]:.6e})")

# --.. ..- .-.. .-.. --- plotting --.. ..- .-.. .-.. ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(eps * 100, sig_ref / 1e6, "-", color="#0072B2", label="material point, same scheme")
ax.plot(eps * 100, sig / 1e6, "o", color="#D55E00", ms=4, label="Z3ST")
ax.axhline(sig_sat / 1e6, color="#009E73", ls="-.", label="exact saturation stress")
ax.set_xlabel(r"strain $\varepsilon_{zz}$ (%)")
ax.set_ylabel(r"stress $\sigma_{zz}$ (MPa)")
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "stress_strain_curve.png"), dpi=200)

# --.. ..- .-.. .-.. --- non-regression metrics --.. ..- .-.. .-.. ---
errors = {
    "sigma_zz_curve": metric(sig[-1], sig_ref[-1], rel=curve_err),
    "cumulative_slip_final": metric(slip[-1], slip_ref[-1]),
    "other_stress_components": metric(lateral, 0.0),
    "saturation_gap": tracked((sig_sat - sig[-1]) / sig_sat),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
