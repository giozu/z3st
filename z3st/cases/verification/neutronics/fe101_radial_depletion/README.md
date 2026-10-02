# FE101 B1 five-ring depletion verification — completed

triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb is a separate, shareable notebook.
It embeds selected definitions from ../fe101_geometry_static/source_model.py,
so no external CLEAN notebook is required or modified. The physical B target
has fuel radius 1.791 cm, near-vacuum He4 gap to 1.804 cm (130 um), and Al
cladding to 1.880 cm. Active z is 10.20–45.76 cm. Only B1 depletes; all other
60 FE101 remain unchanged and nondepletable.

Five equal-area rings use unrounded R*sqrt(i/5). Each initial material is a
fresh exact-composition clone at 6.3 g/cm3, 294 K, with original S(a,b) and
H/Zr=0.994903480489389. One physical instance per material. Total fuel volume
358.346194419294 cm3, mass 2257.581024841551 g, initial U mass
0.180642610509426 kg. Every initial isotope inventory is conserved.

## Completed run

Predictor, 20 x 175 h = 3500 h, whole-reactor power 250000 W,
energy-deposition, diff_burnable_mats=False, write_rates=True, final_step=True.
Preliminary CLEAN settings: 5000 particles, 30 batches, 10 inactive.
Exit code 0; runtime including prepared postprocessing 911.182 s.
All 21 time points/statepoints were saved. Initial/final keff are
1.024662 ±0.003183 and 1.012377 ±0.003041.

| Ring | Final FIMA | Final %FIMA | BU [MWd/kgU] | FIMA/FE mean |
|---:|---:|---:|---:|---:|
| 1 | 0.0001301083 | 0.01301083 | 3.644413 | 0.853181 |
| 2 | 0.0001405951 | 0.01405951 | 3.933373 | 0.921948 |
| 3 | 0.0001509545 | 0.01509545 | 4.217579 | 0.989880 |
| 4 | 0.0001624079 | 0.01624079 | 4.532258 | 1.064985 |
| 5 | 0.0001784231 | 0.01784231 | 4.975189 | 1.170005 |

FE mean FIMA 0.00015249778 (0.01524978%), BU 4.260562 MWd/kg initial U.
Weighted means reconstruct the FE values within 2.5e-13 relative error.
Outer/inner FIMA is about 1.371: a marked outward gradient. Its formal
statistical significance is not computed. zest_FIMA.csv was generated and
validated with 105 time/ring rows; coupling to ZEST was not performed.

## Portable preparation and references

Set OPENMC_CROSS_SECTIONS to the same ENDF/B-VIII.0 cross_sections.xml and
OPENMC_CHAIN_FILE to chain-endf-b8.0.xml. No nuclear-data location is hardcoded.
From this directory, `python build_notebook.py` rebuilds the notebook without
OpenMC execution. `python validate_preparation.py` constructs/checks the model
in a temporary directory, without an operator, transport or depletion.
`python prepare_initial.py` exports fresh initial XMLs under preparation/;
its reusable prepare_initial(export=False) API constructs objects only.

The notebook works from the repository root or its case directory; FE101_CASE_DIR
and FE101_OUTPUT_ROOT can override discovery/output when executing programmatically.
RUN_DEPLETION and RUN_POSTPROCESS default to False. A future authorized rerun
requires explicitly enabling them and using a fresh run directory. Existing
historical outputs are not overwritten by validation.

reference/initial_inventory.json preserves the initial physical/inventory
baseline. reference/radial_history.json has only the 105 records needed by
future 5x5 aggregation comparisons. reference/radial_summary.json records the
run outcome/FE means. These small JSONs require no ignored outputs or HDF5.
Generated XMLs, logs, CSV outputs, caches and large run directories stay local.

## FIMA, energetic BU and limits

FIMA = cumulative fissions / fixed initial U+Zr atoms (H excluded); %FIMA=100*FIMA.
Fission rates sum every chain nuclide with a fission reaction, including U238
and actinides produced later. Results.get_reaction_rate returns fissions/s.
BU = ring deposited energy / (8.64e10 * initial U mass in kg), kept separate
from FIMA. heating-local is normalized to whole-reactor 250 kW; no fixed
200 MeV/fission conversion and no assignment of 250 kW to one FE or ring.
First-order BOS/left-endpoint integration matches this Predictor verification.
Final transport supplies final rates without another irradiation interval.
Statistical and time-integration uncertainties are not propagated.
