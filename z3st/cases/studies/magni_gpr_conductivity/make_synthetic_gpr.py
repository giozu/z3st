#!/usr/bin/env python3
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: build a synthetic GPR conductivity checkpoint for verification
# Author: Giovanni Zullo
# Version: 0.3.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""Fit the GPR conductivity hook to a known analytic residual.

This is the GPR counterpart of the neural-network case's `train_knet.py`: the
checkpoint is regenerated from scratch before every run, so nothing binary is
committed and any clone reproduces it bit for bit.

The samples are drawn from `synthetic_residual.py`, not from measurements.
What is verified is the machinery -- kernel evaluation, de-standardisation,
the Newton tangent, the external-operator plumbing -- against a residual whose
value and derivative are known in closed form.

Assimilating real measurements is a separate workflow; see `fit_gpr.py`.

Usage (invoked by each GPR case's Allrun):
    python3 make_synthetic_gpr.py
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))

from fit_gpr import FEATURES, fit_gp, standardize          # noqa: E402
from synthetic_residual import T_HI, T_LO, residual        # noqa: E402
from z3st.materials.magni_mox_thermal import k_numpy       # noqa: E402

# Sampling grid. It has to bracket the composition and temperature range the
# GPR cases actually visit (Pu ~0.20 with a small Olander tilt, Am 0, x 0.02,
# p 0.05 and 0.12), with margin, so the solver never evaluates the fit outside
# the fitted region. Am and x are varied even though the residual does not
# depend on them: a fit that invents a dependence there would be a real defect.
GRID = {
    "Temp": np.linspace(T_LO, T_HI, 8),
    "Pu": np.array([0.15, 0.20, 0.30]),
    "Am": np.array([0.00, 0.04]),
    "x": np.array([0.00, 0.04]),
    "p": np.array([0.03, 0.09, 0.15]),
}

# The samples are noise-free by construction, so the noise term only has to
# keep the Cholesky factorisation well conditioned. The lengthscale is wide
# because the target residual is smooth: it was selected by sweeping against
# the analytic truth, and 4.0 minimises the tangent error (see README).
LENGTHSCALE = 4.0
NOISE = 1.0e-3


def build_samples():
    mesh = np.meshgrid(*[GRID[name] for name in FEATURES], indexing="ij")
    X = np.column_stack([m.ravel() for m in mesh])
    cols = {name: X[:, i] for i, name in enumerate(FEATURES)}
    r = residual(cols["Temp"], Pu=cols["Pu"], Am=cols["Am"],
                 x=cols["x"], p=cols["p"])
    k_magni = k_numpy(cols["Temp"], Pu=cols["Pu"], Am=cols["Am"],
                      x=cols["x"], p=cols["p"])
    return X, r, k_magni * np.exp(r)


def main():
    out = Path(__file__).with_name("output")
    out.mkdir(parents=True, exist_ok=True)

    X, r, k_data = build_samples()
    print(f"synthetic samples : {len(X)} over {len(FEATURES)} features")
    print(f"  T   [K]         : {X[:, 0].min():.0f} .. {X[:, 0].max():.0f}")
    print(f"  residual        : {r.min():+.4f} .. {r.max():+.4f} "
          f"(k correction {100*(np.exp(r).min()-1):+.1f}% .. "
          f"{100*(np.exp(r).max()-1):+.1f}%)")
    print(f"  k   [W/m/K]     : {k_data.min():.3f} .. {k_data.max():.3f}")

    Xn, x_mean, x_scale = standardize(X)
    y_mean = float(r.mean())
    y_scale = float(r.std() or 1.0)
    yn = (r - y_mean) / y_scale
    model = fit_gp(Xn, yn, lengthscale=LENGTHSCALE, noise=NOISE)

    path = out / "magni_gpr_model.npz"
    np.savez(
        path,
        X_train=Xn,
        alpha=model["alpha"],
        L=model["L"],
        x_mean=x_mean,
        x_scale=x_scale,
        y_mean=y_mean,
        y_scale=y_scale,
        lengthscales=model["lengthscales"],
        signal_variance=model["signal_variance"],
        noise_variance=model["noise_variance"],
        feature_names=np.array(FEATURES),
    )
    print(f"saved synthetic checkpoint: {path}")


if __name__ == "__main__":
    main()
