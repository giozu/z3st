"""Create a self-contained separate notebook; no OpenMC imports or runs."""
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROVENANCE = json.loads((ROOT.parent/"fe101_geometry_static/source_provenance.json").read_text())
SOURCE = ROOT.parent/"fe101_geometry_static/source_model.py"
DESTINATION = ROOT/"triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb"
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == PROVENANCE["source_model_sha256"]
parts = re.split(r"# Notebook cell (\d+) \(zero-based\)\n", SOURCE.read_text())
source_cells = {int(parts[i]): parts[i+1] for i in range(1,len(parts),2)}
source = {"cells": {i:{"source":text} for i,text in source_cells.items()},
          "metadata": {"kernelspec":{"display_name":"Python 3", "language":"python", "name":"python3"},
                       "language_info":{"name":"python"}}}
cells = []


def add(kind, text):
    cell = {"cell_type": kind, "metadata": {}, "source": text.splitlines(keepends=True),
            "id": hashlib.sha256((str(len(cells))+text).encode()).hexdigest()[:12]}
    if kind == "code":
        cell.update(execution_count=None, outputs=[])
    cells.append(cell)


add("markdown", """# TRIGA FE101 B1 — five-ring radial depletion preparation

Separate copy derived from `triga_single_FE101_B1_1965_2026_CLEAN.ipynb`.
The original notebook is not modified. Full 1965 core: 61 FE101; only five
distinct equal-area B1 materials deplete, other 60 FE101 remain frozen.

**Reference geometry is physical B:** fuel R=1.791 cm; near-vacuum He4 gap
1.791–1.804 cm; Al clad 1.804–1.880 cm in the active segment. Legacy A is
retained only in ordinary, unchanged FE101, not as the target baseline.

Transport and depletion are disabled by default. Preparation performs only
Python object construction, deterministic geometry checks and initial XML export.
""")
add("code", """import json
import re
from pathlib import Path

RUN_DEPLETION = False
RUN_POSTPROCESS = False
EXPORT_PREPARED_INPUTS = True
CASE_DIR = Path(globals().get("FE101_CASE_DIR", Path.cwd()))
if not (CASE_DIR/"radial_geometry.py").is_file():
    CASE_DIR = next(parent/"z3st/cases/verification/neutronics/fe101_radial_depletion"
                    for parent in [Path.cwd(), *Path.cwd().parents]
                    if (parent/"z3st/cases/verification/neutronics/fe101_radial_depletion/radial_geometry.py").is_file())
SOURCE_MODEL = CASE_DIR.parent/"fe101_geometry_static/source_model.py"
OUTPUT_ROOT = Path(globals().get("FE101_OUTPUT_ROOT", CASE_DIR/"preparation"))
PREPARED_DIR = OUTPUT_ROOT/"prepared_inputs"
RUN_DIR = OUTPUT_ROOT/"run_3500h_ED"
""")
add("markdown", "## Original CLEAN transport definitions — unchanged core, materials and rods\n")
# Copy model definitions, omitting the old one-material B1/depletion/results
# workflow. The replacement below creates the dedicated physical B1 universe.
for i in [1, 3, 5, 6, 8, 9, *range(13, 23), 26]:
    text = "".join(source["cells"][i]["source"])
    if i == 8:
        text = text.replace("inline_plot = True", "inline_plot = False")
    add("code", f"# Original CLEAN cell {i} (zero-based)\n"+text)
add("markdown", "## Dedicated physical B1: five rings and preparation audit\n")
add("code", (ROOT/"radial_geometry.py").read_text())
add("markdown", """## Settings and nuclear data

Settings remain those of the CLEAN single-FE101 preliminary verification:
30 batches, 10 inactive, 5000 particles/batch. No high-statistics sensitivity
settings are silently substituted. Original temperatures, rods, source and
nuclear data are retained. The following tally definitions preserve the original
flux and global heating diagnostics and add material-resolved ring diagnostics.
""")
for i in (31, 33, 34):
    add("code", f"# Original CLEAN cell {i} (zero-based)\n"+"".join(source["cells"][i]["source"]))
add("code", """ring_filter = openmc.MaterialFilter(ring_materials)
ring_fission_tally = openmc.Tally(name="FE101 B1 ring total fission")
ring_fission_tally.filters = [ring_filter]
ring_fission_tally.scores = ["fission"]
ring_fission_tally.estimator = "tracklength"
ring_energy_tally = openmc.Tally(name="FE101 B1 ring energy deposition")
ring_energy_tally.filters = [ring_filter]
ring_energy_tally.scores = ["heating-local"]
tallies_file.extend([ring_fission_tally, ring_energy_tally])
radial_model = openmc.Model(geometry=geometry, materials=materials_file,
                            settings=settings_file, tallies=tallies_file)
""")
add("code", f"""audit.update({{
    "source_notebook": {PROVENANCE["notebook"]!r},
    "source_model": "../fe101_geometry_static/source_model.py",
    "source_model_sha256": {PROVENANCE["source_model_sha256"]!r},
    "source_sha256": {PROVENANCE['sha256']!r},
    "cross_sections": xs, "chain_file": chain,
    "settings": {{"batches": settings_file.batches, "inactive": settings_file.inactive,
                 "particles": settings_file.particles}},
    "normalization_mode": "energy-deposition", "diff_burnable_mats": False,
    "power_W": 250000.0, "steps": 20, "step_hours": 175.0,
    "write_rates": True, "final_step": True,
    "FIMA_denominator": "fixed initial U+Zr atoms from exact input composition; H excluded",
    "FIMA_primary_for_ZEST": True,
}})
if not chain:
    raise ValueError("Set OPENMC_CHAIN_FILE to chain-endf-b8.0.xml for depletion preparation")
if EXPORT_PREPARED_INPUTS:
    PREPARED_DIR.mkdir(parents=True, exist_ok=True)
    geometry.export_to_xml(PREPARED_DIR/"geometry.xml")
    materials_file.export_to_xml(PREPARED_DIR/"materials.xml")
    settings_file.export_to_xml(PREPARED_DIR/"settings.xml")
    tallies_file.export_to_xml(PREPARED_DIR/"tallies.xml")
    (PREPARED_DIR/"preparation_audit.json").write_text(json.dumps(audit, indent=2)+"\\n")
    print("Prepared initial XMLs and fixed initial inventories:", PREPARED_DIR)
""")
add("markdown", """## Depletion — disabled until explicitly authorized

Predictor: 20 intervals × 175 h = 3500 h, at **250 kW whole-reactor power**.
Use `energy-deposition`, `diff_burnable_mats=False`, `write_rates=True`,
`final_step=True`. There will be 21 saved time points and 21 transport
evaluations, including the genuine final-rate calculation. Five ring materials
are already distinct; do not enable automatic material differentiation.

Changing `RUN_DEPLETION` is an explicit execution switch. Preparation does not
construct or initialize a depletion operator. A fresh run directory is required
to prevent accidental mixing with old one-material or prior five-ring results.
""")
add("code", """if RUN_DEPLETION:
    import openmc.deplete.pool
    openmc.deplete.pool.USE_MULTIPROCESSING = False
    if RUN_DIR.exists() and any(RUN_DIR.iterdir()):
        raise FileExistsError(f"Use a fresh depletion directory: {RUN_DIR}")
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    if not (PREPARED_DIR/"preparation_audit.json").exists():
        raise FileNotFoundError("Export preparation_audit.json before running")
    operator = openmc.deplete.CoupledOperator(
        radial_model, chain_file=openmc.config["chain_file"],
        normalization_mode="energy-deposition", diff_burnable_mats=False)
    operator.output_dir = RUN_DIR
    assert len(operator.burnable_mats) == 5
    assert set(operator.burnable_mats) == {str(m.id) for m in ring_materials}
    integrator = openmc.deplete.PredictorIntegrator(
        operator, timesteps=[175.0]*20, power=250000.0, timestep_units="h")
    integrator.integrate(
        path="depletion_results_FE101_B1_1965_5RINGS_3500h_ED.h5",
        write_rates=True, final_step=True)
else:
    print("Depletion/transport disabled: preparation only.")
""")
add("markdown", """## Postprocessing plan — inventories, direct FIMA and independent energetic BU

For each ring and each of 21 times, save U235/U238/Pu239/Cs137 inventories
in atoms; total physical fission rate summed over all fissioning chain nuclides;
cumulative fissions; FIMA; deposited power and energy; BU in MWd/kg initial U.

Rates read through `Results.get_reaction_rate()` are already fissions/s, not
per-atom rates. FIMA uses fixed initial U+Zr atoms, excludes H and is a fraction
(0.01 means 1% FIMA). No 200 MeV/fission conversion is used.

For this first Predictor verification, cumulative fissions and deposited energy
use first-order BOS/left-endpoint integration: sum(rate[k]*dt[k]). This is a
time-discretization approximation, not an exact integrated Bateman counter.
The final rate is real because `final_step=True`, but adds no further interval.

BU integrates the material-resolved `heating-local` tally normalized by global
heating at 250 kW and divides by fixed initial U mass. It means **energy
deposited in that ring**, with local secondary-photon deposition assumption;
it is distinct from a fission-energy-production burnup. Capture heating is
included. Neither the FE nor an individual ring is assigned the whole 250 kW.

Outputs in `run_3500h_ED/postprocessing/`:
- `ring_history.csv` and `ring_history.json`: 105 ring/time records.
- `fission_rates_by_nuclide.csv`: isotope-resolved rates used in the sum.
- `zest_FIMA.csv`: primary transfer file, with time, ring radii/r/R and FIMA.

Initial geometry, volumes, U mass and atom-count denominators are preserved in
`prepared_inputs/preparation_audit.json`. All 21 BOS/final statepoints and
the depletion H5 are retained. Statistical/time-integration uncertainties are
not silently claimed to be propagated into cumulative FIMA/BU.
""")
add("code", (ROOT/"postprocess.py").read_text())
add("code", """if RUN_POSTPROCESS:
    ring_history = postprocess_radial(RUN_DIR, PREPARED_DIR/"preparation_audit.json")
    display(ring_history.tail(5))
else:
    print("Postprocessing disabled until the five-ring result exists.")
""")
notebook = {"cells": cells, "metadata": source["metadata"],
            "nbformat": 4, "nbformat_minor": 5}
notebook["metadata"]["fe101_radial_preparation"] = {
    "source": "../fe101_geometry_static/source_model.py", "source_model_sha256": PROVENANCE["source_model_sha256"], "source_sha256": PROVENANCE["sha256"],
    "geometry_baseline": "physical B", "transport_executed": False,
    "depletion_executed": False}
DESTINATION.write_text(json.dumps(notebook, ensure_ascii=False, indent=1)+"\n")
print(DESTINATION)
