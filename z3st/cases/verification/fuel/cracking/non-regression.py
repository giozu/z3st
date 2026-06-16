#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: verification/fuel/cracking

Unit test of the isotropic-softening cracking model (models/cracking_model.py). 
The model rescales the elastic constants from the number of macro-cracks.
"""

import os
import re
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_verification import pass_fail_check, regression_check

CASE = os.path.dirname(__file__)
OUT = os.path.join(CASE, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")
TOLERANCE = 5e-3

# --. correlation/scaling constants --..
mat = yaml.safe_load(open(os.path.join(CASE, "fuel.yaml")))
nu0 = float(mat["nu"])
LHR0, N0, N_INF, TAU = 5.0e3, 1.0, 12.0, 21.0e3
f_nu = (2.0 / 3.0) * (2.0 - nu0) / (2.0 + nu0) / (1.0 - nu0)


def n_of(lhr):
    return 0.0 if lhr < LHR0 else N0 + (N_INF - N0) * (1.0 - np.exp(-(lhr - LHR0) / TAU))


def E_ratio_of(n):
    return f_nu ** n


def nu_of(n):
    return nu0 / (2.0**n + (2.0**n - 1.0) * nu0)

log = open(os.path.join(CASE, "log_z3st.md")).read()
rows = re.findall(
    r"\[cracking\]\s+\w+:\s+LHR_max\s+=\s+([\d.]+)\s+kW/m.*?"
    r"n\s+=\s+([\d.]+)\s+cracks,\s+E_iso/E\s+=\s+([\d.]+),\s+nu_iso\s+=\s+([\d.]+)",
    log,
)
if not rows:
    raise RuntimeError("no '[cracking]' lines found in log_z3st.md")
lhr_kW = np.array([float(r[0]) for r in rows])      # this is LHR_max (kW/m)
n_mod = np.array([float(r[1]) for r in rows])
Er_mod = np.array([float(r[2]) for r in rows])
nu_mod = np.array([float(r[3]) for r in rows])

ip = int(np.argmax(lhr_kW))
lhr_peak = lhr_kW[ip] * 1e3
i20 = int(np.argmin(np.abs(lhr_kW - 20.0)))

print(f"[INFO] parsed {len(rows)} cracking steps; LHR_max peak = {lhr_kW[ip]:.1f} kW/m, "
      f"final = {lhr_kW[-1]:.1f} kW/m")
print(f"[INFO] n(peak)     model={n_mod[ip]:.4f}  formula={n_of(lhr_peak):.4f}")
print(f"[INFO] E_iso/E     model={Er_mod[ip]:.6f}  formula={E_ratio_of(n_of(lhr_peak)):.6f}")
print(f"[INFO] nu_iso      model={nu_mod[ip]:.6f}  formula={nu_of(n_of(lhr_peak)):.6f}")
print(f"[INFO] n(20 kW/m)  model={n_mod[i20]:.4f}  formula(paper~6.6)={n_of(20.0e3):.4f}")
print(f"[INFO] irreversibility: n final={n_mod[-1]:.4f} vs n peak={n_mod[ip]:.4f}")


def _entry(num, ref):
    return {"numerical": float(num), "reference": float(ref),
            "abs_error": float(abs(num - ref)), "rel_error": float(abs(num - ref) / abs(ref))}


errors = {
    "n_at_peak": _entry(n_mod[ip], n_of(lhr_peak)),
    "E_iso_ratio_at_peak": _entry(Er_mod[ip], E_ratio_of(n_of(lhr_peak))),
    "nu_iso_at_peak": _entry(nu_mod[ip], nu_of(n_of(lhr_peak))),
    "n_at_20kW_paper": _entry(n_mod[i20], n_of(20.0e3)),
    "irreversibility_n_held": _entry(n_mod[-1], n_mod[ip]),
}

# --. figure: n and E_iso/E vs LHR --..
lhr_line = np.linspace(0.0, 45.0e3, 300)
n_line = np.array([n_of(x) for x in lhr_line])
fig, ax1 = plt.subplots(figsize=(7, 5))
ax1.plot(lhr_line / 1e3, n_line, "k-", lw=2, alpha=0.7, label="n(LHR) model")
ax1.plot(lhr_kW, n_mod, "o", mfc="none", mec="C0", mew=1.8, label="Z3ST n")
ax1.set_xlabel("rod-average LHR (kW/m)"); ax1.set_ylabel("number of cracks n")
ax1.set_ylim(0, 13)
ax2 = ax1.twinx()
ax2.plot(lhr_line / 1e3, E_ratio_of(n_line), "r--", lw=2, alpha=0.7, label="E_iso/E model")
ax2.plot(lhr_kW, Er_mod, "s", mfc="none", mec="r", mew=1.8, label="Z3ST E_iso/E")
ax2.set_ylabel("E_iso / E"); ax2.set_ylim(0, 1)
ax1.set_title(f"Isotropic-softening cracking (nu0={nu0}): model vs correlation")
ax1.grid(True, ls=":", alpha=0.5)
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="center right", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "cracking_correlation.png"), dpi=150)
print("[INFO] cracking_correlation.png saved")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE)
regression_check(errors, CASE)
print("\n[INFO] non-regression completed.\n")
