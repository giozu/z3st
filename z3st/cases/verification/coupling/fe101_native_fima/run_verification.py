"""Noninteractive A-D gates, then native phase-1 and equal-mean controls.

Each Z3ST invocation has a recorded resource gate. No OpenMC/depletion calls.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[4]
STATUS = CASE/'verification_status.json'


def preflight():
    memory = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        name,value,*_ = line.split()
        memory[name.rstrip(':')] = int(value)*1024
    available = memory['MemAvailable']
    swap_used = memory['SwapTotal']-memory['SwapFree']
    disk = shutil.disk_usage(CASE).free
    stable = available >= 2*1024**3 and disk >= 2*1024**3 and swap_used <= max(1024**3,.5*memory['SwapTotal'])
    return {'RAM_available_bytes':available,'RAM_available_GiB':available/1024**3,
            'swap_total_bytes':memory['SwapTotal'],'swap_used_bytes':swap_used,
            'disk_available_bytes':disk,'disk_available_GiB':disk/1024**3,
            'logical_CPUs':os.cpu_count(),'gate_pass':stable,
            'thresholds':'RAM>=2GiB; disk>=2GiB; swap_used<=max(1GiB,50% total)',
            'OMP_NUM_THREADS':2,'MPI_processes':1,'BLAS_NumExpr_threads':1}


def environment():
    env = os.environ.copy()
    env.update({'PYTHONPATH':str(ROOT),'UCX_TLS':'self','OMP_NUM_THREADS':'2',
                'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1',
                'BLIS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','Z3ST_PLAIN_LOG':'1',
                'MPLCONFIGDIR':str(CASE/'cache/matplotlib'),'XDG_CACHE_HOME':str(CASE/'cache'),
                'OMPI_MCA_btl':'self'})
    return env


def load_fields(mode):
    return np.load(CASE/'runs'/mode/'output/final_fields.npz')


def run():
    state = {'phase':'A-B interface','completed':False,'solves':{},'thermal_power_coupled':False}
    def save():STATUS.write_text(json.dumps(state,indent=2)+'\n')
    save()
    with (CASE/'interface_tests.log').open('w') as log:
        subprocess.run([sys.executable,str(CASE/'test_interface.py')],cwd=CASE,env=environment(),stdout=log,stderr=subprocess.STDOUT,check=True)
    state['A_B'] = 'PASS';save()
    with (CASE/'binding_tests.log').open('w') as log:
        subprocess.run([sys.executable,str(CASE/'test_binding.py')],cwd=CASE,env=environment(),stdout=log,stderr=subprocess.STDOUT,check=True)
    state['FE_binding'] = 'PASS';save()
    def solve(mode):
        state['phase'] = mode
        resources = preflight()
        run_dir = CASE/'runs'/mode
        (run_dir/'resource_preflight.json').write_text(json.dumps(resources,indent=2)+'\n')
        state['solves'][mode] = {'resources_before_solve':resources};save()
        if not resources['gate_pass']:
            raise RuntimeError(f'Resource gate failed before {mode}; no solve launched')
        if (run_dir/'output').exists():
            raise RuntimeError(f'Refusing to overwrite outputs in {mode}')
        started = time.monotonic()
        with (run_dir/'log_z3st.md').open('w') as log:
            result = subprocess.run([sys.executable,'-m','z3st'],cwd=run_dir,env=environment(),stdout=log,stderr=subprocess.STDOUT)
        state['solves'][mode].update({'exit_code':result.returncode,'runtime_s':time.monotonic()-started});save()
        text = (run_dir/'log_z3st.md').read_text()
        if result.returncode or 'proceeding with last-iteration state' in text or '[WARNING] diagnostics.per_step failed' in text:
            raise RuntimeError(f'{mode} failed solve/diagnostics validity gate; inspect log')
        history = [json.loads(row) for row in (run_dir/'output/diagnostic_history.jsonl').read_text().splitlines()]
        target = 1.0 if mode.startswith('uniform_a') or mode == 'uniform_constant_reference' else 12600000.0
        if not history or history[-1]['time_s'] != target:
            raise RuntimeError(f'{mode}: incomplete saved history')
        state['solves'][mode]['final'] = history[-1];save()

    # C: native uniform field vs analytic free swelling and existing constant
    # volumetric-eigenstrain channel, in an isothermal zero-power control.
    solve('uniform_analytic');solve('uniform_constant_reference')
    a,b = load_fields('uniform_analytic'),load_fields('uniform_constant_reference')
    np.testing.assert_allclose(a['vector_coords'],b['vector_coords'],rtol=0,atol=0)
    xyz=a['vector_coords'];fuel=xyz[:,0] <= .01791+1e-12
    expected=2e-4*xyz[fuel,:2]
    analytic_error=float(np.max(abs(a['u'][fuel]-expected)))
    constant_error=float(np.max(abs(a['u']-b['u'])))
    if analytic_error > 1e-9 or constant_error > 1e-9:
        raise RuntimeError('C uniform swelling regression failed (absolute displacement tolerance 1 nm)')
    state['C'] = {'result':'PASS','analytic_displacement_max_abs_error_m':analytic_error,'constant_eigenstrain_control_max_abs_difference_m':constant_error,'absolute_tolerance_m':1e-9};save()

    # D: imported zero excludes legacy BU-driven swelling, despite nonzero BU.
    solve('baseline_no_swelling');solve('zero_import')
    a,b=load_fields('baseline_no_swelling'),load_fields('zero_import')
    np.testing.assert_allclose(a['scalar_coords'],b['scalar_coords'],rtol=0,atol=0)
    terr=float(np.max(abs(a['T']-b['T'])));uerr=float(np.max(abs(a['u']-b['u'])))
    if terr > 1e-6 or uerr > 1e-10 or not np.all(b['FIMA_native']==0):
        raise RuntimeError('D imported-zero thermoelastic baseline regression failed')
    state['D']={'result':'PASS','Tmax_abs_difference_K':terr,'umax_abs_difference_m':uerr,'T_tolerance_K':1e-6,'u_tolerance_m':1e-10};save()
    assert state['A_B']=='PASS' and state['C']['result']=='PASS' and state['D']['result']=='PASS'
    state['A_D_gate_before_native']='PASS';save()
    solve('uniform_mean');solve('native')
    state['phase']='postprocessing';save()
    subprocess.run([sys.executable,str(CASE/'postprocess.py')],cwd=CASE,env=environment(),check=True)
    state['completed']=True;state['phase']='complete';save()


if __name__=='__main__':
    try:
        if '--postprocess-only' in sys.argv:
            state=json.loads(STATUS.read_text())
            assert state['A_D_gate_before_native']=='PASS'
            assert all(v.get('exit_code')==0 and 'final' in v for v in state['solves'].values())
            with (CASE/'binding_tests.log').open('w') as log:
                subprocess.run([sys.executable,str(CASE/'test_binding.py')],cwd=CASE,env=environment(),stdout=log,stderr=subprocess.STDOUT,check=True)
            state['FE_binding']='PASS'
            if 'error' in state:
                state.setdefault('technical_recoveries',[]).append(state.pop('error'))
            STATUS.write_text(json.dumps(state,indent=2)+'\n')
            subprocess.run([sys.executable,str(CASE/'postprocess.py')],cwd=CASE,env=environment(),check=True)
            state['completed']=True;state['phase']='complete'
            STATUS.write_text(json.dumps(state,indent=2)+'\n')
        else:
            run()
    except Exception as error:
        state=json.loads(STATUS.read_text()) if STATUS.exists() else {}
        state.update({'completed':False,'error':repr(error)})
        STATUS.write_text(json.dumps(state,indent=2)+'\n')
        raise
