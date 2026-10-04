"""Isolated phase-2 cases; phase-1 inputs and outputs remain immutable."""
import copy,csv,json,os,shutil
from pathlib import Path
import numpy as np
import yaml
from z3st.coupling.openmc.depletion_fields import HeatingHistory

CASE=Path(__file__).resolve().parent
ROOT=CASE.parents[4]
PHASE1=CASE.parent/'fe101_native_fima'
HIGH=ROOT/'z3st/cases/verification/neutronics/fe101_10x10_depletion/preparation/high_statistics_40000_40_15'
SOURCE=HIGH/'run_3500h_ED/postprocessing/zest_FIMA_2D_10x10.csv'


def prepare():
    h=HeatingHistory(SOURCE)
    tables=CASE/'synthetic_histories';tables.mkdir(exist_ok=True)
    with SOURCE.open() as stream:original=list(csv.DictReader(stream))
    for mode in ['zero','known_uniform','same_power_uniform']:
        rows=copy.deepcopy(original)
        for row in rows:
            t=float(row['time_s']);v=float(row['volume_cm3'])
            total=0 if mode=='zero' else (5000 if mode=='known_uniform' else h.power_at(t).sum())
            row['deposited_power_W']=str(total*v/h.volumes_cm3.sum())
        with (tables/(mode+'.csv')).open('w',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    shutil.copy2(PHASE1/'synthetic_histories/mean.csv',tables/'mean_fima.csv')
    for name in ['mesh.geo','mesh.msh','geometry.yaml','clad.yaml','boundary_conditions.yaml']:
        shutil.copy2(PHASE1/name,CASE/name)
    for mode in ['zero_native','uniform_native','uniform_baseline','native_thermal','native_uniform_fima','uniform_power_native_fima','full_native','phase1_replay']:
        run=CASE/'runs'/mode;run.mkdir(parents=True,exist_ok=True)
        if (run/'output').exists():raise RuntimeError('Refusing to overwrite '+str(run))
        inp=yaml.safe_load((PHASE1/'input.yaml').read_text());fuel=yaml.safe_load((PHASE1/'fuel.yaml').read_text())
        inp['coupling']={}
        fima_mode=mode in ['native_uniform_fima','uniform_power_native_fima','full_native','phase1_replay']
        if fima_mode:
            fi=tables/'mean_fima.csv' if mode=='native_uniform_fima' else SOURCE
            inp['coupling']['native_fima']={'enabled':True,'material':'fuel','axial_origin_cm':10.20,'history_path':os.path.relpath(fi,run)}
        else:
            fuel.pop('eigenstrain');fuel.pop('fima_source')
        if mode not in ['uniform_baseline','phase1_replay']:
            power={'zero_native':tables/'zero.csv','uniform_native':tables/'known_uniform.csv','uniform_power_native_fima':tables/'same_power_uniform.csv'}.get(mode,SOURCE)
            inp['coupling']['native_power']={'enabled':True,'material':'fuel','axial_origin_cm':10.20,'history_path':os.path.relpath(power,run)}
        if mode=='uniform_baseline':inp['lhr']=[5000/.3556]*2
        for name,data in [('input',inp),('fuel',fuel)]:
            (run/(name+'.yaml')).write_text(yaml.safe_dump(data,sort_keys=False))
        for name in ['mesh.msh','geometry.yaml','clad.yaml','boundary_conditions.yaml']:
            shutil.copy2(CASE/name,run/name)
        shutil.copy2(CASE/'diagnostics.py',run/'diagnostics.py')
    (CASE/'source_metadata.json').write_text(json.dumps({'source':str(SOURCE),'sha256':h.sha256,'tally':'FE101 B1 domain energy deposition','score':'heating-local','native_unit':'eV/source particle','normalization_tally':'factor-for-normalization','normalization':'P_i[W]=250000*H_i/H_core; saved deposited_power_W','SI_conversion':'q_i[W/m3]=deposited_power_W/(volume_cm3*1e-6)','interpolation':'linear instantaneous power, same physical time as linear cumulative FIMA; no extrapolation; different from BOS cumulative-energy integration','start_fuel_power_W':float(h.power_W[0].sum()),'final_fuel_power_W':float(h.power_W[-1].sum())},indent=2)+'\n')


if __name__=='__main__':prepare()
