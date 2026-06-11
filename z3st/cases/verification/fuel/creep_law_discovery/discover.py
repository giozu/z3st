#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- Z3ST case: verification/fuel/creep_law_discovery --.. ..- .-.. .-.. ---
"""
Sparse identification of the creep mechanism from finite-element data.

The data-generating "experiment" is this case's own FEM simulation: the
verified axisymmetric stress-relaxation problem (a zircaloy bar held at
constant total strain, ../creep_relaxation) re-run with 500
implicit time steps so that the time-discretisation defect of the data
(~0.2%) sits well below the measurement noise. The mean axial stress history
sigma_zz(t), observed at 51 equally spaced times, is extracted from the XDMF
output and perturbed with 2% multiplicative noise. The inverse model is an
independent material-point integrator on a different time grid (400 steps),
so data generation and identification share no discretisation.

The candidate creep-rate library, in the dimensionless stress S = sigma/sigma_ref,

    eps_dot_cr(sigma) = sum_k c_k phi_k(S),   c_k >= 0,
    phi in { S, S^2, S^3, S^5, sinh(S) },

spans diffusional (linear), power-law (Norton n = 2, 3, 5) and Garofalo
(sinh) mechanisms. The true mechanism in the data is the cubic Norton term.

Identification: the relaxation ODE sigma' = -E * eps_dot_cr(sigma) is
integrated by backward Euler; forward-mode automatic differentiation (dual
numbers, implemented below) propagates d sigma / d log c_k through every
Newton-corrected implicit step; Gauss-Newton least squares fits the full
library; mechanisms are then eliminated backwards in order of their share of
the accumulated creep strain, each reduced library refitted, and the final
model chosen by the one-standard-error parsimony rule (the sparsest model on
the elimination path whose loss lies within one standard error of the best) -
sparse selection in the spirit of EUCLID and STLSQ, implemented independently
(no EUCLID GPL-3.0 code is used).

Run:   python3 discover.py        (after the upstream case has produced
                                   output/fields.xdmf; see ./Allrun)
Out:   output/creep_law_discovery.png, output/discovery.json
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
DATA_CSV = os.path.join(OUT, "fem_stress_history.csv")
OUT_PNG = os.path.join(OUT, "creep_law_discovery.png")
OUT_JSON = os.path.join(OUT, "discovery.json")

# --. physics of the data-generating case (zircaloy_creep.yaml) --..
R_GAS = 8.314462618
E = 99.3e9
A0_TRUE, N_TRUE, Q, T = 2.82e-24, 3.0, 1.2e5, 600.0
A_TRUE = A0_TRUE * np.exp(-Q / (R_GAS * T))
SIGMA0 = E * 5.0e-5 / 0.05
T_TOTAL = 1.0e8
N_FEM_STEPS = 500                      # 501 stored states (input.yaml n_steps)
OBS_EVERY = 10                         # observe every 10th state -> 51 points

SIGMA_REF = 5.0e7                      # reference stress for the library
C3_TRUE = A_TRUE * SIGMA_REF**3        # true coefficient of the cubic term

NOISE = 0.02
SEED = 7
N_BE = 400                             # inverse-model time grid


# ---------------------------------------------------------------------------
# data: mean sigma_zz(t) from the FEM output of the upstream case
# ---------------------------------------------------------------------------
def extract_fem_history():
    """Extract (t_i, sigma_zz_i) at the observation times; cache as CSV."""
    if os.path.isfile(DATA_CSV):
        d = np.loadtxt(DATA_CSV, delimiter=",", skiprows=1)
        return d[:, 0], d[:, 1]
    from z3st.utils.utils_extract_xdmf import extract_field_xdmf
    xdmf = os.path.join(OUT, "fields.xdmf")
    if not os.path.isfile(xdmf):
        raise FileNotFoundError(
            f"{xdmf} not found - run the FEM forward problem first (./Allrun)")
    obs = range(0, N_FEM_STEPS + 1, OBS_EVERY)
    times = np.array([T_TOTAL * k / N_FEM_STEPS for k in obs])
    sig = []
    for k in obs:
        Sk = extract_field_xdmf(xdmf, "Stress", step_index=k, return_coords=False)
        sig.append(float(np.mean(np.asarray(Sk)[:, 8])))
    sig = np.array(sig)
    np.savetxt(DATA_CSV, np.column_stack([times, sig]), delimiter=",",
               header="time_s,sigma_zz_Pa", comments="")
    return times, sig


t_data, sigma_fem = extract_fem_history()
rng = np.random.default_rng(SEED)
sigma_data = sigma_fem * (1.0 + NOISE * rng.standard_normal(sigma_fem.size))
sigma_data[0] = sigma_fem[0]           # the initial state is known exactly


# ---------------------------------------------------------------------------
# forward-mode AD: dual numbers with a K-slot gradient
# ---------------------------------------------------------------------------
class Dual:
    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v = float(v)
        self.g = np.asarray(g, dtype=float)

    def _lift(self, x):
        return x if isinstance(x, Dual) else Dual(x, np.zeros_like(self.g))

    def __neg__(self):
        return Dual(-self.v, -self.g)

    def __add__(self, o):
        o = self._lift(o)
        return Dual(self.v + o.v, self.g + o.g)

    __radd__ = __add__

    def __sub__(self, o):
        o = self._lift(o)
        return Dual(self.v - o.v, self.g - o.g)

    def __rsub__(self, o):
        return self._lift(o) - self

    def __mul__(self, o):
        o = self._lift(o)
        return Dual(self.v * o.v, self.v * o.g + o.v * self.g)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._lift(o)
        return Dual(self.v / o.v, (self.g * o.v - self.v * o.g) / o.v**2)

    def __pow__(self, p):           # fixed (float) exponent only
        return Dual(self.v**p, p * self.v ** (p - 1.0) * self.g)

    def exp(self):
        val = np.exp(self.v)
        return Dual(val, val * self.g)

    def sinh(self):
        return Dual(np.sinh(self.v), np.cosh(self.v) * self.g)

    def cosh(self):
        return Dual(np.cosh(self.v), np.sinh(self.v) * self.g)


# ---------------------------------------------------------------------------
# candidate library  phi_k(S)  and its derivative (for the Newton direction)
# ---------------------------------------------------------------------------
LIBRARY = [
    ("S",      lambda S: S,         lambda S: 1.0 if not isinstance(S, Dual) else Dual(1.0, np.zeros_like(S.g))),
    ("S^2",    lambda S: S**2,      lambda S: 2.0 * S),
    ("S^3",    lambda S: S**3,      lambda S: 3.0 * S**2),
    ("S^5",    lambda S: S**5,      lambda S: 5.0 * S**4),
    ("sinh S", lambda S: S.sinh() if isinstance(S, Dual) else np.sinh(S),
               lambda S: S.cosh() if isinstance(S, Dual) else np.cosh(S)),
]
K_FULL = len(LIBRARY)


def relax_be(z, active):
    """BE relaxation history, differentiable in z_j = ln c_{active[j]}.

    Returns the history sampled at the observation times (the inverse grid
    N_BE is an integer multiple of the 50 observation intervals, but differs
    from the 500-step grid of the data-generating FEM run).
    """
    k_act = len(active)
    cs = [Dual(z[j], np.eye(k_act)[j]).exp() for j in range(k_act)]
    dt = T_TOTAL / N_BE
    keep = N_BE // (N_FEM_STEPS // OBS_EVERY)

    s = Dual(SIGMA0, np.zeros(k_act))
    out = [s]
    for k in range(N_BE):
        x = s
        for _ in range(40):
            rate = sum(c * LIBRARY[i][1](x / SIGMA_REF) for c, i in zip(cs, active))
            drate = sum(c * LIBRARY[i][2](x / SIGMA_REF) for c, i in zip(cs, active))
            g = x + (E * dt) * rate - s
            gp = 1.0 + (E * dt / SIGMA_REF) * drate
            x = x - g / gp
            if abs(g.v) < 1e-7 * SIGMA0:
                break
        s = x
        if (k + 1) % keep == 0:
            out.append(s)
    return out


def residuals_and_jacobian(z, active):
    hist = relax_be(z, active)
    r = np.array([(h.v - d) / d for h, d in zip(hist, sigma_data)])
    J = np.array([h.g / d for h, d in zip(hist, sigma_data)])
    return r, J


def fit(z0, active, n_iter=40):
    """Damped Gauss-Newton on z = ln c (trust region + LM damping)."""
    z = np.asarray(z0, dtype=float).copy()
    r, J = residuals_and_jacobian(z, active)
    loss = 0.5 * r @ r
    lam = 1e-3
    for _ in range(n_iter):
        step = np.linalg.solve(J.T @ J + lam * np.eye(len(z)), -J.T @ r)
        if (nrm := np.linalg.norm(step)) > 1.5:
            step *= 1.5 / nrm
        try:
            r_new, J_new = residuals_and_jacobian(z + step, active)
            loss_new = 0.5 * r_new @ r_new
            if not np.isfinite(loss_new):
                raise FloatingPointError
        except (OverflowError, FloatingPointError, ZeroDivisionError):
            lam *= 10.0
            continue
        if loss_new < loss:
            z, r, J, loss = z + step, r_new, J_new, loss_new
            lam = max(lam / 3.0, 1e-12)
        else:
            lam *= 10.0
        if np.linalg.norm(step) < 1e-10:
            break
    return z, loss


def strain_shares(z, active):
    """Share of the accumulated creep strain carried by each active term."""
    cs = np.exp(z)
    hist = relax_be(z, active)
    S_traj = np.array([h.v for h in hist]) / SIGMA_REF
    phi = np.array([[float(LIBRARY[i][1](S)) for S in S_traj] for i in active])
    incr = (cs[:, None] * phi[:, :-1] + cs[:, None] * phi[:, 1:]) / 2.0
    totals = incr.sum(axis=1)
    return totals / totals.sum()


# ---------------------------------------------------------------------------
# stage 1: seed each candidate alone by a coarse 1D scan
# ---------------------------------------------------------------------------
print(f"[INFO] data: {t_data.size} observations from the {N_FEM_STEPS}-step FEM "
      f"forward problem, {NOISE:.0%} noise, seed {SEED}")
z_single = np.empty(K_FULL)
for i in range(K_FULL):
    grid = np.linspace(np.log(1e-16), np.log(1e-6), 25)
    losses = []
    for zg in grid:
        try:
            r, _ = residuals_and_jacobian([zg], [i])
            losses.append(0.5 * r @ r)
        except (OverflowError, FloatingPointError, ZeroDivisionError):
            losses.append(np.inf)
    z_single[i] = grid[int(np.argmin(losses))]

# stage 2: joint fit of the full library
z_full0 = z_single - np.log(K_FULL)            # split the rate equally
active = list(range(K_FULL))
z_full, loss_full = fit(z_full0, active)
shares_full = strain_shares(z_full, active)
print("[INFO] full-library fit:")
for i, (name, *_), in enumerate(LIBRARY):
    print(f"         c[{name:6s}] = {np.exp(z_full[i]):.3e} 1/s   "
          f"strain share = {shares_full[i]:.1%}")

# stage 3: backward elimination ranked by strain share, then the
# one-standard-error parsimony rule (sparsest model whose loss lies within
# one standard error of the best loss on the path)
path = [(list(active), z_full.copy(), loss_full)]
z_cur, act_cur = z_full.copy(), list(active)
while len(act_cur) > 1:
    shares = strain_shares(z_cur, act_cur)
    drop = act_cur[int(np.argmin(shares))]
    keep = [i for i in act_cur if i != drop]
    z0 = np.array([z_cur[act_cur.index(i)] for i in keep])
    z_cur, loss_cur = fit(z0, keep)
    act_cur = keep
    path.append((list(act_cur), z_cur.copy(), loss_cur))
    print(f"[INFO] eliminated {LIBRARY[drop][0]:6s} -> "
          f"{[LIBRARY[i][0] for i in act_cur]}, loss = {loss_cur:.5f}")

loss_best = min(p[2] for p in path)
se = loss_best * np.sqrt(2.0 / sigma_data.size)   # SE of the loss at the floor
admissible = [p for p in path if p[2] <= loss_best + se]
active, z, loss_final = min(admissible, key=lambda p: len(p[0]))
print(f"[INFO] one-SE rule: best loss {loss_best:.5f}, SE {se:.5f} -> "
      f"sparsest admissible model")

c_id = np.exp(z)
names_id = [LIBRARY[i][0] for i in active]
print(f"[INFO] selected mechanism(s): {names_id}")
if active == [2]:
    n_eq = 3.0
    print(f"[INFO] identified  : c[S^3] = {c_id[0]:.4e} 1/s  "
          f"(true {C3_TRUE:.4e} 1/s, error {abs(c_id[0]/C3_TRUE-1.0):.2%})")
    print("[INFO] the cubic Norton mechanism is recovered; the residual "
          "coefficient error absorbs the data noise and the time-"
          "discretisation defect of the FEM data (model-form discrepancy).")

# ---------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------
os.makedirs(OUT, exist_ok=True)
result = {
    "library": [name for name, *_ in LIBRARY],
    "selected": names_id,
    "coefficients_1_per_s": {n: c for n, c in zip(names_id, c_id.tolist())},
    "c3_true_1_per_s": C3_TRUE,
    "full_fit_coefficients": dict(zip([n for n, *_ in LIBRARY],
                                      np.exp(z_full).tolist())),
    "full_fit_strain_shares": dict(zip([n for n, *_ in LIBRARY],
                                       shares_full.tolist())),
    "loss_full": float(loss_full),
    "loss_final": float(loss_final),
    "selection": "one-standard-error rule on the backward-elimination path",
    "elimination_path": [
        {"active": [LIBRARY[i][0] for i in a], "loss": float(l)}
        for a, _, l in path
    ],
    "noise": NOISE,
    "seed": SEED,
}
with open(OUT_JSON, "w") as f:
    json.dump(result, f, indent=4)
print(f"[INFO] results written: {OUT_JSON}")

# figure: (a) relaxation history, (b) library coefficients before/after
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
days = 86400.0


def curve(zz, act, n_fine=400):
    hist = relax_be(zz, act)
    return np.array([h.v for h in hist])


t_grid = t_data
sig_id = curve(z, active)
ax1.plot(t_data / days, sigma_data / 1e6, "o", ms=6, mfc="none",
         mec="tab:blue", mew=1.5, label="FEM data + 2% noise")
ax1.plot(t_grid / days, sig_id / 1e6, "r-", lw=2.2,
         label=f"identified law ({' + '.join(names_id)})")
tt = np.linspace(0, T_TOTAL, 400)
sig_true = (SIGMA0 ** (1 - N_TRUE) + (N_TRUE - 1) * E * A_TRUE * tt) ** (-1 / (N_TRUE - 1))
ax1.plot(tt / days, sig_true / 1e6, "k--", lw=1.8, alpha=0.7,
         label="true Norton law (closed form)")
ax1.set_xlabel("time (days)")
ax1.set_ylabel("axial stress (MPa)")
ax1.set_title("Stress relaxation: FEM data and identified law")
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend()

x = np.arange(K_FULL)
c_full = np.exp(z_full)
c_sel = np.array([c_id[active.index(i)] if i in active else 0.0
                  for i in range(K_FULL)])
w = 0.38
ax2.bar(x - w / 2, c_full, w, color="lightsteelblue", label="full-library fit")
ax2.bar(x + w / 2, np.where(c_sel > 0, c_sel, np.nan), w, color="tab:red",
        label="after mechanism selection")
ax2.axhline(C3_TRUE, color="k", ls="--", lw=1.2, label="true cubic coefficient")
ax2.set_yscale("log")
ax2.set_xticks(x, [n for n, *_ in LIBRARY])
ax2.set_xlabel("candidate mechanism")
ax2.set_ylabel("coefficient (1/s)")
ax2.set_title("Library coefficients and sparse selection")
ax2.grid(True, axis="y", ls=":", alpha=0.6)
ax2.legend()

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
print(f"[INFO] figure saved: {OUT_PNG}")
