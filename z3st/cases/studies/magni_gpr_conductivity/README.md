# Magni + GPR conductivity update

This folder keeps the data-assimilation workflow separate from the solver.

Workflow:

1. Use `materials.magni_mox_thermal.k_numpy` as the physics baseline.
2. Fit a Gaussian process on the log residual:

   `r = log(k_data / k_magni)`

3. Save a lightweight NumPy checkpoint:

   `output/magni_gpr_model.npz`

4. Use the checkpoint in a material card:

```yaml
k:
  type: gpr
  model: cases/studies/magni_gpr_conductivity/output/magni_gpr_model.npz
  mode: mean

Pu: 0.20
Am: 0.03
Np: 0.0
x: 0.02
p: 0.05
```

The current GPR feature vector is:

`Temp, Pu, Am, x, p`

The CSV currently has no explicit `Np` or burnup column, so the update is a
fresh-fuel correction on top of the full Magni formula.

Note: as for the existing neural-network conductivity card, a relative `model`
path is resolved from the run directory. Adjust it when copying the material
card into another case.
