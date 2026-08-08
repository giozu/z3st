# GPR update of the Magni conductivity

This folder keeps the data-assimilation workflow separate from the solver: the
GPR hook corrects the published Magni MA-MOX correlation with a Gaussian
process on the log residual

    r = log(k_data / k_magni)      ->      k = k_magni * exp(r)

There are two distinct workflows here, and it matters which one you are in.

## 1. Verification — synthetic, and what the test suite runs

`make_synthetic_gpr.py` fits the hook to a residual defined in closed form by
`synthetic_residual.py`. Nothing here is experimental and nothing here is a
physical claim: it is a deliberately smooth, bounded fiction (a ~10%
correction) whose only purpose is to give the checks an exact answer to
compare against.

This is the same pattern as the neural-network case, whose `train_knet.py`
fits the NN hook to the analytic law `k(T) = 1/(a + b*T)` rather than to
measurements.

What it verifies is the **machinery**, not the physics of any correction:

- kernel evaluation and de-standardisation of features and target
- the Newton tangent `dk/dT` against the analytic derivative
- `value_and_grad` against a finite difference of `__call__` — this catches a
  tangent that is internally consistent but wired to the wrong operand
- that the fit tracks the residual's `Pu` and `p` dependence and does **not**
  invent one on `Am` or `x`, which the residual does not contain

`verify_machinery.py` runs those checks and exits non-zero on failure. Both
scripts are invoked by each GPR case's `Allrun`, so the checkpoint is rebuilt
from scratch before every run and **no binary artefact is committed** —
`*.npz` is gitignored.

Fit settings (`LENGTHSCALE = 4.0`, 8 temperatures, noise-free samples) were
chosen by sweeping against the analytic truth. They reach ~2e-3 relative error
on `k` and ~7e-3 on the tangent for the compositions the cases run. Shorter
lengthscales are worse here: the target residual is smooth, and the error is
dominated by the isotropic kernel spanning five standardised dimensions rather
than by resolution in T.

## 2. Assimilation of real data

`fit_gpr.py` is the reference implementation for fitting the hook to
measurements. It takes `--csv` explicitly; **no dataset is distributed with
this repository**, so running it requires obtaining source data separately.

One caveat, stated here because it is a property of the method rather than of
any particular dataset: a Gaussian process is non-parametric, so **the
training points are part of the model**. Any checkpoint fitted to real
measurements carries those measurements — the inputs are recoverable by
inverting the stored standardisation, and the targets by recomputing the
kernel matrix from `X_train` and the stored hyperparameters. Do not commit a
checkpoint fitted to data you are not free to publish. To ship a correction
without the points, distil the GP mean into a parametric form and distribute
the coefficients instead.

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
card and by the hook but is not currently a feature, so the correction is a
fresh-fuel one on top of the full Magni formula.

`mode: mean` is the deterministic path and the only one the cases exercise.
`mode: affine` adds `xi * sigma` for UQ sweeps; it needs the stored Cholesky
factor `L`, which is otherwise unused, and its Newton tangent deliberately
uses the mean derivative.

As for the neural-network card, a relative `model` path is resolved from the
run directory — adjust it when copying the card into another case.
