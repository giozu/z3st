"""One isolated resource-gated solve; never overwrite existing results."""
import importlib.util,json,os,subprocess,sys,time
from pathlib import Path
CASE=Path(__file__).resolve().parent
ROOT=CASE.parents[4]
def main():
    run=CASE/'runs/swelling_off'
    if (run/'output').exists():raise RuntimeError('Refusing existing output')
    path=CASE.parent/'fe101_native_power/run_verification.py'
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location('phase2_runner_readonly',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    resources=mod.preflight(run)
    if not resources['gate_pass']:raise RuntimeError('Resource gate failed; no solve started')
    env=mod.environment();env.update({'MPLCONFIGDIR':str(CASE/'cache/matplotlib'),'XDG_CACHE_HOME':str(CASE/'cache')})
    start=time.monotonic()
    with (run/'log_z3st.md').open('w') as f:r=subprocess.run([sys.executable,'-m','z3st'],cwd=run,env=env,stdout=f,stderr=subprocess.STDOUT)
    text=(run/'log_z3st.md').read_text()
    status={'runtime_s':time.monotonic()-start,'exit_code':r.returncode,'resources':resources,'PASS':r.returncode==0 and 'diagnostics.per_step failed' not in text and 'proceeding with last-iteration state' not in text}
    (CASE/'solve_status.json').write_text(json.dumps(status,indent=2))
    if not status['PASS']:raise RuntimeError('Solve failed; see log/status')
if __name__=='__main__':main()
