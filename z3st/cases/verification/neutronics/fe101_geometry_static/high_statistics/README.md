# High-statistics static FE101 geometry sensitivity

Independent static eigenvalue calculations, no depletion. The original CLEAN
notebook and the first A/B series are preserved. See ../source_provenance.json
for the source notebook/hash and ../prepare.py for reproducible XML preparation.

| Case | Active fuel R [cm] | Gap [micrometres] | Density [g/cm3] | Seed |
|---|---:|---:|---:|---:|
| A | 1.791 | none, legacy Al | 6.3 | 101 |
| B | 1.791 | 130 | 6.3 | 202 |
| C | 1.7975 | 65 | 6.254519099119661 | 303 |
| D | 1.804 | 0 | 6.209528929307131 | 404 |

In B/C/D the active clad spans 1.804–1.880 cm. B/C use pure He4 at
1e-10 g/cm3 and 294 K. D uses a shared fuel/clad interface surface. Changes
are limited to B1 active fuel/cladding; all axial boundaries, end regions,
other FE, compositions, temperatures and nuclear data are retained.
The active target fuel mass is 2257.581024841552 g in all four cases.

Each run uses 5000 particles, 220 batches, 20 inactive: 1.1 million histories,
1 million active, 200 active batches. Other settings are retained. Five equal-
area radial tally bins use r/R=sqrt(j/5), j=0..5. XMLs remain unchanged during
transport. Logs, exit codes, statepoint.220.h5 and results.json are in case_A
through case_D. Runs use the existing OpenMC 0.16.0 environment and ENDF/B-VIII.0
cross-section path in each manifest.

The original ../read_results.py reads each final statepoint. analyze.py augments
the reader JSON with normalized coordinates, fission density, uncertainties
and comparisons A-B, C-B, D-B in comparison.json. The B1 fission rate is scaled
to 250 kW whole-reactor power using the global heating-local tally, not a
constant 200 MeV/fission. No depletion model/operator/integrator is executed.

Raw tally standard errors and keff errors come from OpenMC. Covariances needed
for profiles f_i/F and physical rates F/heating are absent from final statepoints.
Ratio uncertainties are first-order estimates with covariance set to zero;
additional bounds use |Cov(X,Y)| <= sigma_X*sigma_Y. These bounds are not
confidence intervals. Distinct seeds justify independent-run error propagation,
but source convergence and batch autocorrelations are not independently tested.
Results are conditional on fixed temperatures, this single target and the
near-vacuum He assumption; they do not assess simultaneous full-core deformation
or temperature/composition feedback.
