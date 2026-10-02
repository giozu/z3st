# FE101 B1 neutronic geometry verification — completed

The frozen `source_model.py` contains the original CLEAN full 1965 core
(61 FE101) definitions. The source notebook was not changed and is not needed
at runtime; its name/hash and snapshot provenance are in source_provenance.json.
Only B1 active fuel/gap/clad changes in the sensitivity cases. Rod positions,
source, axial reflectors, samarium disks and all other elements remain unchanged.

| Case | Fuel radius [cm] | He gap [cm] | Al clad [cm] |
|---|---:|---|---|
| A | 1.791 | none | 1.791–1.880 |
| B (baseline) | 1.791 | 1.791–1.804 | 1.804–1.880 |
| C | 1.7975 | 1.7975–1.804 | 1.804–1.880 |
| D | 1.804 | none | 1.804–1.880 |

C/D preserve the original fuel mass by adjusting density. B preserves the
initial FE101 composition and 6.3 g/cm3 density. Gap He4 is near-vacuum,
1e-10 g/cm3 at 294 K, a sensitivity assumption rather than a measured inventory.
Fuel active z is 10.20–45.76 cm. Five equal-area tally rings use R*sqrt(i/5).
No depletion was run in this static study.

## Completed results

Preliminary A/B runs used 5000 particles, 30 batches, 10 inactive. A–D
high-statistics runs used 5000 particles, 220 batches, 20 inactive and distinct
seeds. Those high-statistics settings are confined to the sensitivity study.

| Case | Preliminary keff ± 1 sigma | High-statistics keff ± 1 sigma |
|---|---|---|
| A | 1.024921 ± 0.002375 | 1.018665 ± 0.000944 |
| B | 1.019921 ± 0.002765 | 1.017535 ± 0.000978 |
| C | — | 1.019200 ± 0.001025 |
| D | — | 1.019629 ± 0.001060 |

High-statistics delta keff relative to B: A +113 ±136 pcm, C +167 ±142 pcm,
D +209 ±144 pcm. These differences are below 1.5 sigma. The B radial fission
profile, normalized to FE mean, is [0.876446, 0.921888, 0.996179, 1.057541,
1.147946]. All saved ring sums reconstruct the integrated FE tally.
Compact reference/geometry_results.json contains only the relevant results.
Raw statepoints, generated XMLs and logs are deliberately excluded from Git.

## Reproduction

Use OpenMC with the same ENDF/B-VIII.0 HDF5 data. Set OPENMC_CROSS_SECTIONS to
its cross_sections.xml; set OPENMC_CHAIN_FILE to chain-endf-b8.0.xml when using
these definitions for depletion preparation. Paths are supplied by the user,
not embedded in the shared source. Python dependencies include numpy, pandas,
matplotlib and OpenMC; high_statistics/analyze.py also uses h5py.

From this directory, `python prepare.py --case B` constructs inputs and checks
geometry only. `--high-statistics` explicitly selects the sensitivity settings.
Running transport is a separate action. `python read_results.py case_B` reads
an existing result. `high_statistics/analyze.py` requires locally regenerated
A–D outputs; the versioned compact reference needs no HDF5 to inspect.

Tally uncertainties are OpenMC standard errors. Normalized ratios require
covariance; estimates assuming zero covariance are approximate. Independent
seeds do not establish source convergence or remove batch correlations.
The axial verification in ../fe101_axial_static uses the same physical B
baseline and preliminary statistics.
