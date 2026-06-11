#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST non-regression script --.. ..- .-.. .-.. ---
"""
Z3ST case: V_creep_relaxation_verification

Verifies implicit Norton creep on stres relaxation at constant total strain

Closed form:  sigma(t) = [sigma_0^(1-n) + (n-1)*E*A*t]^(-1/(n-1)),
with sigma_0 = E*u0/Lz and A = A0*exp(-Q/RT).

Checks:
  1. sigma_end vs the scalar backward-Euler replica of the same recursion, 
     sigma_{k+1} + E*dt*A*sigma_{k+1}^n = sigma_k,
     solved per step by a scalar Newton.
  2. sigma_end vs the exact closed form
     shrinks as O(dt) with n_steps (re-run with more steps to see it drop).

One figure: sigma(t), Z3ST per-step vs scalar BE vs the exact solution.
"""

import os
import glob
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from z3st.utils.utils_verification import *
from z3st.utils.utils_extract_xdmf import *

R_GAS = 8.314462618

CASE_DIR = os.path.dirname(__file__)
OUT = os.path.join(CASE_DIR, "output")
OUT_JSON = os.path.join(OUT, "non-regression.json")

XDMF_FILE = os.path.join(OUT, "fields.xdmf")

# --. case parameters --..
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    geom = yaml.safe_load(f)
Lz = float(geom["Lz"])

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    inp = yaml.safe_load(f)
t_total = float(max(inp["time"]))
n_steps = int(inp["n_steps"]) - 1

with open(os.path.join(CASE_DIR, next(iter(inp["materials"].values())))) as f:
    mat = yaml.safe_load(f)
E = float(mat["E"])
A0, n, Q = float(mat["creep_A0"]), float(mat["creep_n"]), float(mat["creep_Q"])
T = float(mat["T_initial"])

with open(os.path.join(CASE_DIR, "boundary_conditions.yaml")) as f:
    bcs = yaml.safe_load(f)
u0 = float(next(b["displacement"] for b in bcs["mechanical"]["clad"]
                if b["type"] == "Dirichlet_y"))

# --. references --..
A = A0 * np.exp(-Q / (R_GAS * T))
sigma0 = E * u0 / Lz

def sigma_exact(t):
    return (sigma0 ** (1.0 - n) + (n - 1.0) * E * A * np.asarray(t)) ** (-1.0 / (n - 1.0))

def be_replica(num_steps):
    """Scalar backward-Euler recursion sigma_{k+1} + E*dt*A*sigma_{k+1}^n = sigma_k."""
    dt = t_total / num_steps
    s = sigma0
    out = [s]
    for _ in range(num_steps):
        x = s
        for _ in range(60):
            g = x + E * dt * A * x**n - s
            gp = 1.0 + E * dt * A * n * x ** (n - 1.0)
            x = x - g / gp
        s = x
        out.append(s)
    return np.array(out)

SIGMA_BE = be_replica(n_steps)
SIGMA_BE_END = float(SIGMA_BE[-1])
SIGMA_EXACT_END = float(sigma_exact(t_total))
TOLERANCE = 1e-2

# --. numerical result: mean sigma_zz over the (uniform) bar --..
_, _, _, S = extract_field_xdmf(XDMF_FILE, field_name="Stress", step_index=-1)
sigma_end = float(np.mean(np.asarray(S)[:, 8]))

be_defect = abs(SIGMA_BE_END - SIGMA_EXACT_END) / SIGMA_EXACT_END
print(f"[INFO] sigma_0 = {sigma0/1e6:.2f} MPa, A(T) = {A:.4e}, n = {n:.0f}, "
      f"{n_steps} steps")
print(f"[INFO] sigma_end: Z3ST = {sigma_end/1e6:.4f} MPa, "
      f"scalar BE = {SIGMA_BE_END/1e6:.4f} MPa, exact = {SIGMA_EXACT_END/1e6:.4f} MPa")
print(f"[INFO] backward-Euler defect vs exact at this dt: {be_defect:.2%} "
      f"(shrinks O(dt) with more steps)")

errors = {
    "sigma_relaxed_vs_backward_euler": {
        "numerical": sigma_end,
        "reference": SIGMA_BE_END,
        "abs_error": float(abs(sigma_end - SIGMA_BE_END)),
        "rel_error": float(abs(sigma_end - SIGMA_BE_END) / SIGMA_BE_END),
    },
    "time_discretisation_defect": {
        "numerical": float(abs(sigma_end - SIGMA_EXACT_END) / SIGMA_EXACT_END),
        "reference": float(be_defect),
        "abs_error": float(abs(abs(sigma_end - SIGMA_EXACT_END) / SIGMA_EXACT_END - be_defect)),
        "rel_error": float(
            abs(abs(sigma_end - SIGMA_EXACT_END) / SIGMA_EXACT_END - be_defect) / be_defect
        ),
    },
}

# --. figure: relaxation history --..
try:
    times = np.linspace(0.0, t_total, n_steps + 1)
    sig_steps = []
    for k in range(n_steps + 1):
        Sk = extract_field_xdmf(XDMF_FILE, "Stress", step_index=k, return_coords=False)
        sig_steps.append(float(np.mean(np.asarray(Sk)[:, 8])))
    t_line = np.linspace(0.0, t_total, 300)
    plt.figure(figsize=(7, 5))
    plt.plot(t_line / 86400.0, sigma_exact(t_line) / 1e6, "k-", lw=2.5, alpha=0.7,
             label="Exact  [σ₀^{1−n} + (n−1)EAt]^{−1/(n−1)}")
    plt.plot(times / 86400.0, SIGMA_BE / 1e6, "s--", ms=5, color="tab:blue",
             alpha=0.7, label=f"scalar backward Euler ({n_steps} steps)")
    if len(sig_steps) == len(times):
        plt.plot(times / 86400.0, np.array(sig_steps) / 1e6, "o", ms=8, mfc="none",
                 mec="r", mew=2, label="Z3ST (implicit, AD tangent)")
    plt.xlabel("time (days)")
    plt.ylabel("axial stress σ_zz (MPa)")
    plt.title("Norton creep stress relaxation at constant strain\n"
              f"σ₀ = {sigma0/1e6:.0f} MPa, T = {T:.0f} K, n = {n:.0f}")
    plt.grid(True, ls=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "stress_relaxation.png"), dpi=150)
    print("[INFO] stress_relaxation.png saved")
except Exception as e:
    print(f"[WARNING] relaxation plot skipped: {type(e).__name__}: {e}")

pass_fail_check(errors, TOLERANCE, OUT_JSON, CASE_DIR)
regression_check(errors, CASE_DIR)

print("\n[INFO] non-regression completed.\n")
