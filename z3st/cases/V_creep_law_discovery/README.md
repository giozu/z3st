# V_creep_law_discovery — sparse identification of the creep mechanism

Identifies the creep mechanism of the verified stress-relaxation problem from
noisy FEM data, by sparse selection over a library of candidate creep laws.
This is the framework's first inverse / constitutive-identification case, a
preliminary step toward EUCLID-style automated model discovery (Flaschel et
al., CMAME 2022 — independent implementation; the published EUCLID codes are
GPL-3.0 and are not used).

## Setup

1. Forward problem (data generation): the axisymmetric Norton-creep stress
   relaxation of `../V_creep_relaxation_verification`, re-run here with 500
   implicit time steps so the time-discretisation defect of the data (~0.2%)
   is well below the noise. Solved by the production FEM code (dolfinx).
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

The cubic Norton mechanism is selected alone (10/10 noise seeds), with the
coefficient recovered within 2% of the true A*sigma_ref^3 (1.63% at the
reference seed). The spurious library terms land 4-6 decades below.

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
