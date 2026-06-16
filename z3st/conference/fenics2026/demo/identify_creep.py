#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST demo: creep-law identification --.. ..- .-.. .-.. ---
"""
Gradient-based identification of Norton creep parameters from relaxation data.

A first, preliminary step toward EUCLID-style automated constitutive-law
discovery (Flaschel et al., CMAME 2022) in Z3ST: today this identifies the
*parameters* (A, n) of a known law from data; library-based sparse regression
over candidate energies is the roadmap.

The physics is the verified stress-relaxation case
(z3st/cases/verification/fuel/creep_relaxation): a bar held at constant total
strain relaxes by Norton creep, sigma' = -E*A*sigma^n, integrated with the
same implicit backward-Euler recursion the FEM solver uses,

    sigma_{k+1} + E*dt*A*sigma_{k+1}^n = sigma_k .

The demo message, in one sentence: *the time integrator is differentiable* —
forward-mode automatic differentiation (dual numbers, implemented below in
~40 lines, no external AD library) propagates d(sigma)/d(theta) through every
Newton-corrected implicit step, and Gauss-Newton least squares recovers the
parameters from noisy synthetic data.

Note on EUCLID: the published EUCLID codes (github.com/EUCLID-code) are
GPL-3.0; nothing here is derived from them. This file is an independent,
from-scratch demonstration of the underlying idea.

Run:  python3 identify_creep.py          (numpy + matplotlib only, ~1 s)
Out:  baked/creep_identification.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PNG = os.path.join(HERE, "baked", "creep_identification.png")

# ---------------------------------------------------------------------------
# forward-mode AD: dual numbers with a 2-slot gradient (d/dlogA, d/dn)
# ---------------------------------------------------------------------------
class Dual:
    """value + gradient w.r.t. the parameter vector theta."""

    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v = float(v)
        self.g = np.asarray(g, dtype=float)

    @staticmethod
    def lift(x, nparams=2):
        return x if isinstance(x, Dual) else Dual(x, np.zeros(nparams))

    def __add__(self, o):
        o = Dual.lift(o)
        return Dual(self.v + o.v, self.g + o.g)

    __radd__ = __add__

    def __sub__(self, o):
        o = Dual.lift(o)
        return Dual(self.v - o.v, self.g - o.g)

    def __rsub__(self, o):
        return Dual.lift(o) - self

    def __mul__(self, o):
        o = Dual.lift(o)
        return Dual(self.v * o.v, self.v * o.g + o.v * self.g)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Dual.lift(o)
        return Dual(self.v / o.v, (self.g * o.v - self.v * o.g) / o.v**2)

    def __pow__(self, o):
        # x**y with both x and y possibly dual:  d = x^y (y' ln x + y x'/x)
        o = Dual.lift(o)
        val = self.v**o.v
        return Dual(val, val * (o.g * np.log(self.v) + o.v * self.g / self.v))

    def exp(self):
        val = np.exp(self.v)
        return Dual(val, val * self.g)


# ---------------------------------------------------------------------------
# the physics: verified relaxation case, verification/fuel/creep_relaxation
# ---------------------------------------------------------------------------
R_GAS = 8.314462618
E = 99.3e9            # Pa        (zircaloy_creep.yaml)
A0_TRUE = 2.82e-24    # Pa^-n/s
N_TRUE = 3.0
Q = 1.2e5             # J/mol
T = 600.0             # K
SIGMA0 = E * 5.0e-5 / 0.05      # = E*u0/Lz = 99.3 MPa
T_TOTAL = 1.0e8       # s  (~3.2 years)

A_TRUE = A0_TRUE * np.exp(-Q / (R_GAS * T))           # ~1.0e-34

# Identify theta = (log10 ghat, n) with the law written around a reference
# stress,  sigma' = -E * ghat * (sigma/sigma_ref)^n,  ghat = A*sigma_ref^n.
# This decorrelates the rate amplitude from the exponent (with raw (A, n) the
# least-squares valley is nearly degenerate and Gauss-Newton crawls).
SIGMA_REF = 5.0e7      # Pa, mid-range of the relaxation data
GHAT_TRUE = A_TRUE * SIGMA_REF**N_TRUE
THETA_TRUE = np.array([np.log10(GHAT_TRUE), N_TRUE])


def sigma_exact(t, A=A_TRUE, n=N_TRUE):
    """Closed form: sigma(t) = [sigma0^(1-n) + (n-1) E A t]^(-1/(n-1))."""
    return (SIGMA0 ** (1.0 - n) + (n - 1.0) * E * A * np.asarray(t)) ** (-1.0 / (n - 1.0))


def relax_be(theta, n_steps, keep_every=1):
    """Backward-Euler relaxation history, differentiable in theta=(log10 ghat, n).

    Each implicit step is solved by Newton *in dual arithmetic*, so the
    parameter sensitivities ride through the solve (implicit function theorem
    for free). Returns the sampled history as a list of Duals.
    """
    lggh = Dual(theta[0], np.array([1.0, 0.0]))
    n = Dual(theta[1], np.array([0.0, 1.0]))
    gh = (lggh * np.log(10.0)).exp()
    dt = T_TOTAL / n_steps

    s = Dual.lift(SIGMA0)
    out = [s]
    for k in range(n_steps):
        x = s
        for _ in range(40):
            g = x + E * dt * gh * (x / SIGMA_REF) ** n - s
            gp = 1.0 + E * dt * gh * n * (x / SIGMA_REF) ** (n - 1.0) / SIGMA_REF
            x = x - g / gp
            if abs(g.v) < 1e-6 * SIGMA0:
                break
        s = x
        if (k + 1) % keep_every == 0:
            out.append(s)
    return out


# ---------------------------------------------------------------------------
# synthetic experiment: the closed form + multiplicative noise
# ---------------------------------------------------------------------------
N_STEPS = 400          # BE grid of the differentiable model (discretisation
                       # defect ~0.3%, well below the 2% data noise)
KEEP = 8               # data live on every 8th step -> 50 points + t=0
rng = np.random.default_rng(42)
t_data = np.linspace(0.0, T_TOTAL, N_STEPS // KEEP + 1)
sigma_data = sigma_exact(t_data) * (1.0 + 0.02 * rng.standard_normal(t_data.size))
sigma_data[0] = SIGMA0          # the initial stress is known exactly


# ---------------------------------------------------------------------------
# Gauss-Newton least squares on theta = (log10 A, n)
# ---------------------------------------------------------------------------
def residuals_and_jacobian(theta):
    # relative residuals: the synthetic noise is multiplicative, so this is
    # the statistically consistent weighting (and it values the late-time
    # tail, whose shape is what identifies the exponent n)
    hist = relax_be(theta, N_STEPS, keep_every=KEEP)
    r = np.array([(h.v - d) / d for h, d in zip(hist, sigma_data)])
    J = np.array([h.g / d for h, d in zip(hist, sigma_data)])
    return r, J


# coarse scan to seed the gradient refinement (the landscape has local
# valleys in n if the start is decades off), then Gauss-Newton polishes.
# Scan range: relaxation rate ghat = A*sigma_ref^n from ~centuries to ~hours
# (log10 ghat in [-14, -8] 1/s), exponent n in the physical range 1.5-5.5.
SCAN_LG = np.linspace(-14.0, -8.0, 7)
SCAN_N = np.linspace(1.5, 5.5, 5)
THETA_INIT = min(
    (np.array([lg, n]) for lg in SCAN_LG for n in SCAN_N),
    key=lambda th: 0.5 * (lambda r: r @ r)(residuals_and_jacobian(th)[0]),
)
theta = THETA_INIT.copy()
history = []
lam = 1e-3
r, J = residuals_and_jacobian(theta)
loss = 0.5 * r @ r
for it in range(30):
    history.append((theta.copy(), loss))
    H = J.T @ J + lam * np.eye(2)
    step = np.linalg.solve(H, -J.T @ r)
    if (norm := np.linalg.norm(step)) > 1.5:   # trust region in theta units
        step *= 1.5 / norm
    try:
        r_new, J_new = residuals_and_jacobian(theta + step)
        loss_new = 0.5 * r_new @ r_new
        if not np.isfinite(loss_new):
            raise FloatingPointError
    except (OverflowError, FloatingPointError, ZeroDivisionError):
        lam *= 10.0
        continue
    if loss_new < loss:                  # accept, relax damping
        theta, r, J, loss = theta + step, r_new, J_new, loss_new
        lam = max(lam / 3.0, 1e-12)
    else:                                # reject, increase damping
        lam *= 10.0
    if np.linalg.norm(step) < 1e-10:
        break
history.append((theta.copy(), loss))

n_id = theta[1]
ghat_id = 10.0 ** theta[0]
A_id = ghat_id / SIGMA_REF**n_id                    # back to Norton A
print(f"[INFO] data: {t_data.size} points, 2% noise, seed 42")
print(f"[INFO] true        : n = {N_TRUE:.4f}  ghat = {GHAT_TRUE:.4e}  A = {A_TRUE:.4e}")
print(f"[INFO] identified  : n = {n_id:.4f}  ghat = {ghat_id:.4e}  A = {A_id:.4e}   "
      f"({len(history)-1} Gauss-Newton iterations)")
print(f"[INFO] error       : n {abs(n_id-N_TRUE)/N_TRUE:.2%},  "
      f"rate-at-sigma_ref {abs(ghat_id/GHAT_TRUE-1.0):.2%},  "
      f"log10(A) {abs(np.log10(A_id/A_TRUE)):.2f} dec "
      f"(A and n are log-correlated — the EUCLID papers regularise exactly this)")

# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

t_line = np.linspace(0.0, T_TOTAL, 400)
days = 86400.0
A_init = 10.0 ** THETA_INIT[0] / SIGMA_REF ** THETA_INIT[1]
ax1.plot(t_line / days, sigma_exact(t_line) / 1e6, "k-", lw=2.5, alpha=0.65,
         label=f"true law  (n = {N_TRUE:.0f})")
ax1.plot(t_line / days, sigma_exact(t_line, A_init, THETA_INIT[1]) / 1e6,
         "--", color="grey", lw=1.8,
         label=f"coarse-scan seed  (n = {THETA_INIT[1]:.1f})")
ax1.plot(t_line / days, sigma_exact(t_line, A_id, n_id) / 1e6, "r-.", lw=2.2,
         label=f"identified  (n = {n_id:.2f})")
ax1.plot(t_data / days, sigma_data / 1e6, "o", ms=7, mfc="none", mec="tab:blue",
         mew=1.8, label="noisy data (2%)")
ax1.set_xlabel("time (days)")
ax1.set_ylabel("axial stress σ_zz (MPa)")
ax1.set_title("Norton creep relaxation — parameters recovered from data")
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend()

losses = [h[1] for h in history]
ax2.semilogy(range(len(losses)), losses, "o-", color="tab:red")
ax2.set_xlabel("Gauss-Newton iteration")
ax2.set_ylabel("least-squares loss  ½‖r‖²  (MPa²)")
ax2.set_title("Convergence — gradients by forward-mode AD\n"
              "through the implicit backward-Euler solver")
ax2.grid(True, ls=":", alpha=0.6)
ax2.annotate(f"identified (true):\n"
             f"n = {n_id:.3f}  ({N_TRUE:.0f})\n"
             f"rate at 50 MPa = {ghat_id:.2e}  ({GHAT_TRUE:.2e}) 1/s",
             xy=(0.35, 0.62), xycoords="axes fraction",
             bbox=dict(boxstyle="round", fc="lightyellow", ec="tab:red"))

fig.suptitle("Toward automated constitutive-law discovery — "
             "gradient-based identification (EUCLID-style, preliminary)",
             fontsize=12)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=150)
print(f"[INFO] figure saved: {OUT_PNG}")
