# verification/fuel/creep_law_discovery: sparse identification of the creep mechanism

Identifies the creep mechanism of the verified stress-relaxation problem from
noisy FEM data, by sparse selection over a library of candidate creep laws.
It is an inverse (constitutive-identification) case in the style of EUCLID
automated model discovery (Flaschel et al., CMAME 2022). The implementation is
independent: the published EUCLID codes are GPL-3.0 and are not used.

Excluded from the local suite (`suite_exclude.txt`). It carries a gold.

## Setup

1. Forward problem (data generation): the axisymmetric Norton-creep stress
   relaxation of `../creep_relaxation`, re-run here with `n_steps: 501`. The
   time-discretisation defect of the data is about 0.2 %, below the 2 % noise.
   Solved by Z3ST (dolfinx).
2. Observations: the mean axial stress at 51 equally spaced times, perturbed
   with 2% multiplicative Gaussian noise (fixed seed).
3. Inverse model: a material-point backward-Euler integrator of the
   relaxation ODE sigma' = -E * sum_k c_k phi_k(sigma/sigma_ref) on a
   different time grid (400 steps), with the candidate library

       phi in { S, S^2, S^3, S^5, sinh S },   S = sigma/sigma_ref,

   spanning diffusional, Norton (n = 2, 3, 5) and Garofalo mechanisms.
   The true mechanism in the data is the cubic Norton term.
4. Identification: forward-mode automatic differentiation (dual numbers)
   propagates the parameter sensitivities through every Newton-corrected
   implicit step; damped Gauss-Newton fits log-coefficients; mechanisms are
   eliminated backwards by their share of the accumulated creep strain and
   the final model is chosen by the one-standard-error parsimony rule.

## Result

From `output/discovery.json`:
- The cubic Norton term is selected alone for 10 of 10 noise seeds.
- Its coefficient is within 1.8 % of the true A*sigma_ref^3 over the 10 seeds
  (1.63 % at the reference seed 7).
- In the full five-term fit, S^5 = 1.79e-13 1/s against S^3 = 1.26e-11 1/s
  (about 1.9 decades below, 2.4 % of the creep strain). S, S^2 and sinh S are
  about 6 decades below S^3.

## Run

```bash
./Allrun          # gmsh + z3st (forward FEM, ~35 min) + discover.py + checks
# or, if output/fields.xdmf or the cached CSV already exists:
python3 discover.py
python3 non-regression.py
```

Outputs: `output/creep_law_discovery.png` (paper figure),
`output/discovery.json` (selection, coefficients, elimination path),
`output/fem_stress_history.csv` (cached observations).
