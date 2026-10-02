# FE101 B1 axial and radial-axial static verification — completed

Static eigenvalue transport only; no depletion and no new fuel materials.
The physical B initial baseline is reconstructed from the radial verification
sources, with its exact initial composition, gap/clad and unchanged ordinary
core. This needs ../fe101_radial_depletion sources, not generated XMLs or an
external notebook. Preliminary CLEAN: 5000 particles, 30 batches, 10 inactive.

Fuel R=1.791 cm, He4 gap to 1.804 cm (130 um), Al clad to 1.880 cm.
Axial tally boundaries: 10.200, 17.312, 24.424, 31.536, 38.648, 45.760 cm.
A second cylindrical tally mesh crosses these with five equal-area radial bins.
The discretization is only in the tallies, not in physical fuel geometry.
Each axial volume is 71.669238883859 cm3; each 2D bin 14.333847776772 cm3.

## Results

Exit code 0, transport runtime 19.098 s, keff=1.020778 ±0.003241 (1 sigma).
Densities below use whole-reactor heating-local normalization at 250 kW.

| Axial zone | z [cm] | Fission density [fissions/s/cm3] | Power density [W/cm3] | Fission/FE mean | Outer/inner ring |
|---:|---|---:|---:|---:|---:|
| 1 | 10.200–17.312 | 3.263790e11 | 10.1495 | 0.664793 | 1.579276 |
| 2 | 17.312–24.424 | 5.695057e11 | 17.7040 | 1.160011 | 1.240038 |
| 3 | 24.424–31.536 | 6.274926e11 | 19.4956 | 1.278123 | 1.380307 |
| 4 | 31.536–38.648 | 5.469351e11 | 16.9951 | 1.114038 | 1.330184 |
| 5 | 38.648–45.760 | 3.844295e11 | 11.9405 | 0.783034 | 1.396092 |

Axial maximum/minimum fission density=1.922589. The center exceeds the mean
of the two end zones by 76.557%; it exceeds the lower/upper ends by 92.259%
and 63.227%. Both the five axial tally sums and 25-bin sums reconstruct the
integrated FE fission and heating tallies within 3e-16 relative error.
Volumes and baseline geometry checks passed.

The strong axial amplitude, combined with the radial gradient, justifies a
future 5x5 depletion (conclusion C). Radial shape varies with height in the
observed results, but preliminary statistics and absent covariance propagation
do not establish the significance of that variation or mesh convergence.

## Shared reference and reproduction

reference/static_5x5.json contains the axial table, keff, checks and BOTH 5x5
fission and power matrices normalized to FE mean. This is the complete small
reference required by future benchmarks; no statepoint is needed to read it.
The extraction used existing completed results, without rerunning transport.

Set OPENMC_CROSS_SECTIONS and OPENMC_CHAIN_FILE to the same ENDF/B-VIII.0
libraries used by the radial case. `python run_static.py` prepares and runs a
new static calculation; it rejects nonempty run_preliminary/. Do not invoke
it merely to inspect the reference. `python postprocess.py` processes an
existing local statepoint. Large outputs and generated XML/log/summary files
are ignored. The initial failed launch was a missing cross-sections environment
setting, resolved before the successful transport; it did not change physics.
