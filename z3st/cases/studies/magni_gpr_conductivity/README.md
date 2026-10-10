# GPR update of the Magni conductivity

This folder holds the data-assimilation workflow, separate from the solver. The
GPR hook corrects the published Magni MA-MOX correlation with a Gaussian
process on the log residual

    r = log(k_data / k_magni)      ->      k = k_magni * exp(r)

There are two workflows: synthetic verification and assimilation of real data.

## 1. Synthetic verification (run by the suite)

`make_synthetic_gpr.py` fits the hook to a residual defined in closed form by
`synthetic_residual.py`. The residual is synthetic, smooth and bounded (a
correction of about 10 %), so the checks have an exact answer.

The neural-network case follows the same pattern: its `train_knet.py` fits the
NN hook to the analytic law `k(T) = 1/(a + b*T)`.

Checks of the machinery:

- kernel evaluation and de-standardisation of features and target
- the Newton tangent `dk/dT` against the analytic derivative
- `value_and_grad` against a finite difference of `__call__`
- the fit tracks the `Pu` and `p` dependence of the residual and shows none on
  `Am` or `x`, which the residual does not contain

`verify_machinery.py` runs those checks and exits non-zero on failure. The
`Allrun` of the three GPR verification cases under
`verification/fuel/thermal_conductivity/` calls both scripts.
`studies/gpr_uq_margin/Allrun` calls `make_synthetic_gpr.py` only. The checkpoint
is rebuilt before every run and is not committed (`*.npz` is gitignored).

Fit settings: `LENGTHSCALE = 4.0`, 8 temperatures, noise-free samples, chosen
by a sweep against the analytic truth. Relative error about 2e-3 on `k` and
7e-3 on the tangent for the compositions the cases run. Shorter lengthscales
give larger errors: the error is dominated by the isotropic kernel spanning five
standardised dimensions, not by resolution in T.

## 2. Assimilation of real data

`fit_gpr.py` fits the hook to measurements. It takes `--csv`. No dataset is
distributed with this repository.

A Gaussian process is non-parametric: the training points are part of the
model. A checkpoint fitted to real measurements carries them. The inputs are
recoverable by inverting the stored standardisation, the targets by recomputing
the kernel matrix from `X_train` and the stored hyperparameters. Do not commit a
checkpoint fitted to data you are not free to publish. To distribute a
correction without the points, fit the GP mean with a parametric form and
distribute its coefficients.

## Material card

```yaml
k:
  type: gpr
  model: ../../../../studies/magni_gpr_conductivity/output/magni_gpr_model.npz
  mode: mean

Pu: 0.20
Am: 0.00
Np: 0.00
x: 0.02
p: 0.05
```

The GPR feature vector is `Temp, Pu, Am, x, p`. `burnup` is accepted on the
card and by the hook but is not a feature, so the correction is a fresh-fuel
one on top of the full Magni formula.

`mode: mean` is the deterministic path, used by the verification cases.
`mode: affine` adds `xi * sigma` for UQ sweeps (used by `studies/gpr_uq_margin`).
It needs the stored Cholesky factor `L`, unused otherwise. Its Newton tangent
uses the mean derivative.

As for the neural-network card, a relative `model` path is resolved from the
run directory. Adjust it when copying the card into another case.
