# FE101 B1 5 radial x 5 axial depletion — prepared, NOT run

triga_single_FE101_B1_1965_5x5_CLEAN.ipynb is a separate notebook, constructed
from the verified initial radial-case sources. It requires no previously
exported XMLs and no external CLEAN notebook. Only active B1 fuel is partitioned;
all other cells/materials, radial surfaces, active end planes, gap, clad,
reflectors, samarium disks and rods are preserved.

25 independent materials Fuel101_B1_r1_z1 through Fuel101_B1_r5_z5 retain the
exact initial FE101 composition, 6.3 g/cm3, 294 K, original S(a,b), and
H/Zr=0.994903480489389. Distinct IDs and one physical instance per material.
All 60 ordinary FE101 and all non-target materials remain nondepletable.

Radial edges use unrounded 1.791*sqrt(i/5). Axial boundaries are 10.200,
17.312, 24.424, 31.536, 38.648, 45.760 cm. Domain volume
14.333847776772 cm3; total 358.346194419294 cm3; fuel mass
2257.581024841551 g; initial U mass 0.180642610509426 kg.
Fuel/U mass and every initial nuclide inventory are conserved; maximum recorded
relative isotope discrepancy is below 3e-16. Five-domain aggregates reconstruct
each radial ring and axial section. Point-location and non-target snapshots
verify the 130 um He gap and 1.804–1.880 cm Al clad.

## Execution state and configuration

No OpenMC transport, operator initialization or depletion has been executed.
RUN_DEPLETION, RUN_POSTPROCESS and RUN_BENCHMARKS are all False. Preparation
and synthetic postprocessing validations have passed. The prior radial
verification completed 3500 h and the static check found center/end-average
fission density +76.557%, motivating this two-dimensional preparation.

Future Predictor: whole-reactor 250000 W, 20 x 175 h, energy-deposition,
diff_burnable_mats=False, write_rates=True, final_step=True. Preliminary CLEAN
5000 particles, 30 batches, 10 inactive; no high-statistics substitution.
A future authorized run requires a fresh run_3500h_ED/ and saves 21 time points.

## Portable use and validations

Set OPENMC_CROSS_SECTIONS to the ENDF/B-VIII.0 cross_sections.xml and
OPENMC_CHAIN_FILE to chain-endf-b8.0.xml. These are environment-specific data
locations, not shared source paths. Notebook discovery supports the repository
root or this directory; FE101_CASE_DIR/FE101_OUTPUT_ROOT can override them.

`python build_notebook.py` rebuilds the notebook only.
`python validate_preparation.py` constructs/checks the model in a temporary
output directory. `python validate_postprocessing_synthetic.py` checks outputs
and aggregation on synthetic data in a removed temporary directory. Neither
validation runs transport/depletion or overwrites historical run artifacts.
`python validate_reference_data.py` checks the small versioned benchmark datasets.
Generated XMLs/audits, logs, cache and run directories are ignored.

## Future postprocessing and benchmark interfaces

525 domain/time records will contain r/z bounds, volume, U235/U238/Pu239/Cs137,
all fissioning-chain-nuclide rates, summed fissions/s, cumulative fissions,
heating-local, local deposited power/energy, FIMA, %FIMA and independent BU.
FIMA uses fixed initial U+Zr atoms, H excluded. BU uses deposited energy and
fixed initial U mass; no conversion of BU into the transferred FIMA.
Both cumulative integrals use first-order BOS left endpoints, consistent with
the prior Predictor verification; final transport adds no extra interval.
Statistical/time-integration uncertainties are not propagated.

The future zest_FIMA_2D.csv contains time_h, radial_index, axial_index,
r_in_cm, r_out_cm, z_low_cm, z_high_cm and FIMA fraction. Its schema is versioned;
no real transfer field has been fabricated and no ZEST coupling performed.

benchmarks.py now reads only versioned small references:
- ../fe101_radial_depletion/reference/radial_history.json: 105 time/ring records;
- ../fe101_axial_static/reference/static_5x5.json: initial fission/power matrices.

It aggregates five axial domains per ring, reports differences from the prior
radial run, and verifies weighted FE FIMA/BU means against total fissions/initial
atoms and deposited energy/initial U mass. Initial FIMA is zero; pattern
comparisons use dFIMA/dt=initial fission rate/initial U+Zr atoms.
Reference-zero relative differences are undefined, stored as NaN in diagnostic
CSV outputs. Agreement is not required to be exact and is not a formal
statistical test. No benchmark depends on ignored statepoints/H5/CSV files.
