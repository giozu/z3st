"""Build separate notebook with execution switches disabled."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent/'fe101_radial_depletion'
reference=json.loads((BASE/'triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb').read_text())
cells=[]
def add(kind,text):
    c={'cell_type':kind,'metadata':{},'source':text.splitlines(keepends=True),'id':hashlib.sha256((str(len(cells))+text).encode()).hexdigest()[:12]}
    if kind=='code': c.update(execution_count=None,outputs=[])
    cells.append(c)
add('markdown','''# FE101 B1 1965 — 5 radial × 5 axial CLEAN depletion preparation

Separate notebook based on verified initial radial-case sources. Physical B,
25 independent fuel materials, original composition/density/temperature/S(a,b).
All other cells and materials preserved. No OpenMC transport/depletion during
preparation. Preliminary CLEAN: 5000 particles, 30 batches, 10 inactive.

Radial edges use 1.791*sqrt(i/5) without rounding. Axial edges are 10.200,
17.312, 24.424, 31.536, 38.648, 45.760 cm. Each material occurs once.
''')
add('code','''from pathlib import Path
RUN_DEPLETION = False
RUN_POSTPROCESS = False
RUN_BENCHMARKS = False
EXPORT_PREPARED_INPUTS = True
CASE_DIR = Path(globals().get('FE101_CASE_DIR', Path.cwd()))
if not (CASE_DIR/'prepare_geometry.py').is_file():
    CASE_DIR = next(parent/'z3st/cases/verification/neutronics/fe101_5x5_depletion'
                    for parent in [Path.cwd(), *Path.cwd().parents]
                    if (parent/'z3st/cases/verification/neutronics/fe101_5x5_depletion/prepare_geometry.py').is_file())
RADIAL_CASE_DIR = CASE_DIR.parent/'fe101_radial_depletion'
RADIAL_NOTEBOOK = RADIAL_CASE_DIR/'triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb'
OUTPUT_ROOT = Path(globals().get('FE101_OUTPUT_ROOT', CASE_DIR/'preparation'))
PREPARED_DIR = OUTPUT_ROOT/'prepared_inputs'
RUN_DIR = OUTPUT_ROOT/'run_3500h_ED'
''')
add('code',(ROOT/'prepare_geometry.py').read_text())
add('markdown','''## Depletion — disabled

Predictor, 20 × 175 h at 250000 W whole-reactor power; energy-deposition,
diff_burnable_mats=False, write_rates=True, final_step=True. This creates
21 time points, including genuine final rates. Do not differentiate materials:
the 25 physical domains already have distinct material IDs. Require a fresh
run directory. No operator is constructed unless RUN_DEPLETION is enabled.
''')
add('code','''if RUN_DEPLETION:
    import openmc.deplete
    import openmc.deplete.pool
    openmc.deplete.pool.USE_MULTIPROCESSING = False
    if RUN_DIR.exists() and any(RUN_DIR.iterdir()):
        raise FileExistsError(f'Use a fresh depletion directory: {RUN_DIR}')
    if not (PREPARED_DIR/'preparation_audit.json').is_file():
        raise FileNotFoundError('Export preparation audit before running')
    RUN_DIR.mkdir(parents=True,exist_ok=True)
    operator = openmc.deplete.CoupledOperator(radial_model,chain_file=audit['chain_file'],
                 normalization_mode='energy-deposition',diff_burnable_mats=False)
    operator.output_dir=RUN_DIR
    assert len(operator.burnable_mats)==25
    assert set(operator.burnable_mats)=={str(m.id) for m in domain_materials}
    integrator=openmc.deplete.PredictorIntegrator(operator,timesteps=[175.]*20,
                     power=250000.,timestep_units='h')
    integrator.integrate(path='depletion_results_FE101_B1_1965_5x5_3500h_ED.h5',
                         write_rates=True,final_step=True)
else:
    print('Transport/depletion disabled: preparation only.')
''')
add('markdown','''## Future postprocessing and benchmarks — disabled

525 domain/time records. FIMA = cumulative fissions / fixed initial U+Zr
atoms, excluding H; %FIMA = 100*FIMA. Sum all chain nuclides with fission,
including U238 and actinides formed later. Rates are physical fissions/s.
BU = heating-local deposited energy / (8.64e10 * initial U mass in kg),
kept separate from FIMA. Both integrations use first-order BOS left endpoints,
consistent with Predictor; final transport supplies rates, not another interval.
Statistical and time-integration errors are not propagated.

Future zest_FIMA_2D.csv has times, radial/axial indices, radii, z bounds and
FIMA fraction. It is not fabricated during preparation.

Benchmarks aggregate along z to compare with the existing radial depletion,
compare initial fission/power patterns with the static 5×5 tally, and check
weighted FE FIMA/BU reconstruction. FIMA(0)=0, so the initial FIMA pattern
uses dFIMA/dt=fission_rate/N_initial, not normalization of zero inventories.
''')
add('code',(ROOT/'postprocess.py').read_text())
add('code',(ROOT/'benchmarks.py').read_text())
add('code','''if RUN_POSTPROCESS:
    domain_history=postprocess_2d(RUN_DIR,PREPARED_DIR/'preparation_audit.json')
else:
    print('Postprocessing disabled until actual 5x5 results exist.')
if RUN_BENCHMARKS:
    benchmark_summary=compare_benchmarks(RUN_DIR,PREPARED_DIR/'preparation_audit.json')
else:
    print('Post-run benchmarks disabled.')
''')
nb={'cells':cells,'metadata':reference['metadata'].copy(),'nbformat':4,'nbformat_minor':5}
nb['metadata'].pop('fe101_radial_preparation',None)
nb['metadata']['fe101_5x5_preparation']={'baseline':'../fe101_radial_depletion','geometry':'physical B','transport_executed':False,'depletion_executed':False}
(ROOT/'triga_single_FE101_B1_1965_5x5_CLEAN.ipynb').write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n')
