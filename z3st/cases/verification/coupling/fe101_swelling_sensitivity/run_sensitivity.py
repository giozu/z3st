"""Sequential resource-gated five-case verification; no output overwrites."""
import importlib.util,json,os,subprocess,sys,time
from pathlib import Path
CASE=Path(__file__).resolve().parent;ROOT=CASE.parents[4]
os.environ['MPLCONFIGDIR']=str(CASE/'cache/matplotlib');os.environ['XDG_CACHE_HOME']=str(CASE/'cache')
def main():
    path=CASE.parent/'fe101_native_power/run_verification.py';sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location('readonly_phase2_runner',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    status={'completed':False,'solves':{}}
    for s in [0.,.5,1.,1.5,2.]:
        run=CASE/'runs'/f's_{s:.1f}'
        if (run/'output').exists():raise RuntimeError('Refusing existing output '+str(run))
        resources=mod.preflight(run)
        if not resources['gate_pass']:raise RuntimeError('Resource gate failed; no solve started')
        env=mod.environment();env.update({'PYTHONPATH':str(CASE)+os.pathsep+str(ROOT),'MPLCONFIGDIR':str(CASE/'cache/matplotlib'),'XDG_CACHE_HOME':str(CASE/'cache')})
        start=time.monotonic()
        with (run/'log_z3st.md').open('w') as f:r=subprocess.run([sys.executable,'-m','z3st'],cwd=run,env=env,stdout=f,stderr=subprocess.STDOUT)
        log=(run/'log_z3st.md').read_text();passed=r.returncode==0 and 'diagnostics.per_step failed' not in log and 'proceeding with last-iteration state' not in log
        status['solves'][str(s)]={'runtime_s':time.monotonic()-start,'exit_code':r.returncode,'resources':resources,'PASS':passed}
        (CASE/'solve_status.json').write_text(json.dumps(status,indent=2))
        print(f's={s}: PASS={passed}',flush=True)
        if not passed:raise RuntimeError('Solve/diagnostics failure '+str(run))
    status['completed']=True;(CASE/'solve_status.json').write_text(json.dumps(status,indent=2))
if __name__=='__main__':main()
