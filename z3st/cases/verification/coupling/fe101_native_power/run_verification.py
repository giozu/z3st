"""Resource-gated A-D tests before the fully coupled one-way verification."""
import json,os,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
from prepare import CASE,ROOT,PHASE1

STATUS=CASE/'verification_status.json'

def environment():
    e=os.environ.copy();e.update({'PYTHONPATH':str(ROOT),'PYTHONDONTWRITEBYTECODE':'1','UCX_TLS':'self','OMPI_MCA_btl':'self','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','BLIS_NUM_THREADS':'1','Z3ST_PLAIN_LOG':'1','MPLCONFIGDIR':str(CASE/'cache/matplotlib'),'XDG_CACHE_HOME':str(CASE/'cache')});return e

def preflight(run):
    commands=[['free','-h'],['df','-h','/'],['df','-h','/mnt/c'],['nproc'],['vmstat','1','3']]
    results=[]
    for cmd in commands:
        r=subprocess.run(cmd,capture_output=True,text=True);results.append({'command':cmd,'exit_code':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
    mem={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines()}
    avail=mem['MemAvailable'];swap=mem['SwapTotal']-mem['SwapFree'];disk=shutil.disk_usage(run).free;windows=shutil.disk_usage('/mnt/c').free
    vm=results[-1]['stdout'].splitlines();rates=[(int(row.split()[6]),int(row.split()[7])) for row in vm[3:] if len(row.split())>=17]
    active_swap=any(si>1024 or so>1024 for si,so in rates)
    gate=all(r['exit_code']==0 for r in results) and avail>=2*1024**3 and disk>=2*1024**3 and swap<=max(1024**3,.5*mem['SwapTotal']) and not active_swap
    data={'commands':results,'RAM_available_bytes':avail,'RAM_available_GiB':avail/1024**3,'swap_used_bytes':swap,'swap_total_bytes':mem['SwapTotal'],'disk_free_bytes':disk,'disk_free_GiB':disk/1024**3,'mnt_c_free_bytes':windows,'logical_CPUs':os.cpu_count(),'gate_pass':gate,'OMP_NUM_THREADS':2,'MPI_processes':1,'BLAS_NumExpr_threads':1,'gate':'RAM and disk >=2GiB; swap_used<=max(1GiB,50% total); no sustained vmstat si/so >1024KiB/s'}
    (run/'resource_preflight.json').write_text(json.dumps(data,indent=2)+'\n');return data

def fields(mode):return np.load(CASE/'runs'/mode/'output/final_fields.npz')

def run():
    s=json.loads(STATUS.read_text()) if '--resume' in sys.argv and STATUS.exists() else {'completed':False,'solves':{}}
    def save():STATUS.write_text(json.dumps(s,indent=2)+'\n')
    with (CASE/'interface_tests.log').open('w') as f:subprocess.run([sys.executable,str(CASE/'test_interface.py')],env=environment(),stdout=f,stderr=subprocess.STDOUT,check=True)
    s['A']='PASS';save()
    with (CASE/'binding_tests.log').open('w') as f:subprocess.run([sys.executable,str(CASE/'test_binding.py')],cwd=CASE,env=environment(),stdout=f,stderr=subprocess.STDOUT,check=True)
    s['binding']='PASS';save()
    def solve(mode):
        if s['solves'].get(mode,{}).get('PASS'):return
        path=CASE/'runs'/mode
        if (path/'output').exists():raise RuntimeError('Refusing overwrite of incomplete outputs: '+str(path))
        s['phase']=mode;resources=preflight(path);s['solves'][mode]={'resources':resources};save()
        if not resources['gate_pass']:raise RuntimeError('Resource gate failed '+mode)
        start=time.monotonic()
        with (path/'log_z3st.md').open('w') as f:r=subprocess.run([sys.executable,'-m','z3st'],cwd=path,env=environment(),stdout=f,stderr=subprocess.STDOUT)
        text=(path/'log_z3st.md').read_text();s['solves'][mode].update({'runtime_s':time.monotonic()-start,'exit_code':r.returncode});save()
        if r.returncode or 'proceeding with last-iteration state' in text or 'diagnostics.per_step failed' in text:raise RuntimeError('Solve/diagnostic failure '+mode)
        rows=[json.loads(line) for line in (path/'output/diagnostic_history.jsonl').read_text().splitlines()]
        assert len(rows)==21 and rows[-1]['time_s']==12600000
        assert all(np.isfinite([v for v in row.values() if isinstance(v,(float,int))]).all() for row in rows)
        s['solves'][mode].update({'PASS':True,'final':rows[-1]});save()
    solve('zero_native')
    a=fields('zero_native');assert np.max(abs(a['T']-300))<1e-5
    s['B']={'PASS':True,'max_abs_T_minus_BC_K':float(np.max(abs(a['T']-300)))};save()
    solve('uniform_native');solve('uniform_baseline')
    a,b=fields('uniform_native'),fields('uniform_baseline')
    np.testing.assert_array_equal(a['scalar_coords'],b['scalar_coords']);err=float(np.max(abs(a['T']-b['T'])))
    assert err<1e-5
    s['C']={'PASS':True,'max_abs_T_difference_K':err,'tolerance_K':1e-5,'known_power_W':5000};save()
    solve('native_thermal');s['D']='PASS';s['A_D_before_fully_coupled']='PASS';save()
    solve('phase1_replay')
    a=fields('phase1_replay');b=np.load(PHASE1/'runs/native/output/final_fields.npz')
    terr=float(np.max(abs(a['T']-b['T'])));uerr=float(np.max(abs(a['u']-b['u'])))
    assert terr<1e-6 and uerr<1e-10
    s['phase1_compatibility']={'PASS':True,'T_error_K':terr,'u_error_m':uerr};save()
    solve('uniform_power_native_fima');solve('native_uniform_fima');solve('full_native')
    subprocess.run([sys.executable,str(CASE/'postprocess.py')],cwd=CASE,env=environment(),check=True)
    s['completed']=True;s['phase']='complete';s.pop('error',None);save()

if __name__=='__main__':
    try:run()
    except Exception as e:
        s=json.loads(STATUS.read_text()) if STATUS.exists() else {};s['error']=repr(e);s['completed']=False;STATUS.write_text(json.dumps(s,indent=2)+'\n');raise
