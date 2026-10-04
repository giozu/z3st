"""Prepare isolated phase-1 verification cases; never run neutronics."""
import copy
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import yaml

from z3st.coupling.openmc.depletion_fields import DepletionHistory

CASE = Path(__file__).resolve().parent
PROJECT = CASE.parents[3]
REFERENCE = PROJECT/'cases/regression/101_rod_2D'
SOURCE = PROJECT/'cases/verification/neutronics/fe101_10x10_depletion/preparation/high_statistics_40000_40_15/run_3500h_ED/postprocessing/zest_FIMA_2D_10x10.csv'


def synthetic_table(history, path, kind):
    columns = ['time_s','radial_index','axial_index','r_in_cm','r_out_cm','z_low_cm','z_high_cm','volume_cm3','initial_U_Zr_atoms','cumulative_fissions','FIMA']
    times = [0.0, 1.0] if kind == 'analytic' else history.times
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for t in times:
            if kind == 'zero':
                value = 0.0
            elif kind == 'analytic':
                value = 0.0 if t == 0 else 2.0e-4
            else:
                value = np.sum(history.fissions_at(t))/np.sum(history.initial_atoms)
            for key, bounds, volume, atoms in zip(history.keys,history.bounds_cm,history.volumes_cm3,history.initial_atoms):
                writer.writerow([t,*key,*bounds,volume,atoms,value*atoms,value])


def prepare():
    history = DepletionHistory(SOURCE)
    assert history.shape == (10,10) and len(history.times) == 21 and history.times[-1] == 3500*3600
    tables = CASE/'synthetic_histories'
    tables.mkdir(exist_ok=True)
    for kind in ['zero','mean','analytic']:
        synthetic_table(history,tables/(kind+'.csv'),kind)
    inp = yaml.safe_load((REFERENCE/'input.yaml').read_text())
    geom = yaml.safe_load((REFERENCE/'geometry.yaml').read_text())
    fuel = yaml.safe_load((REFERENCE/'fuel.yaml').read_text())
    clad = yaml.safe_load((REFERENCE/'clad.yaml').read_text())
    bc = yaml.safe_load((REFERENCE/'boundary_conditions.yaml').read_text())
    geom['Lz'] = 0.3556
    inp['time'] = [0.0,12600000.0]
    inp['n_steps'] = [20]
    inp['output'] = {'format':'vtu'}
    inp['coupling'] = {'native_fima': {'enabled':True,'material':'fuel','axial_origin_cm':10.20,'history_path':os.path.relpath(SOURCE,CASE)}}
    fuel['fima_source'] = 'native_openmc'
    for name,data in [('input',inp),('geometry',geom),('fuel',fuel),('clad',clad),('boundary_conditions',bc)]:
        (CASE/(name+'.yaml')).write_text(yaml.safe_dump(data,sort_keys=False))
    geo = (REFERENCE/'mesh.geo').read_text().replace('h1    = 0.356;', 'h1    = 0.3556;').replace('h2    = 0.356;', 'h2    = 0.3556;')
    (CASE/'mesh.geo').write_text(geo)
    with (CASE/'log_mesh.md').open('w') as log:
        subprocess.run(['/usr/bin/gmsh','mesh.geo','-2','-nt','1'],cwd=CASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    for mode in ['baseline_no_swelling','zero_import','uniform_analytic','uniform_constant_reference','uniform_mean','native']:
        run = CASE/'runs'/mode
        run.mkdir(parents=True,exist_ok=True)
        if (run/'output').exists():
            raise RuntimeError(f'Refusing to overwrite previous solve: {run}')
        config,material,al = copy.deepcopy(inp),copy.deepcopy(fuel),copy.deepcopy(clad)
        if mode == 'baseline_no_swelling':
            config.pop('coupling');material.pop('eigenstrain');material.pop('fima_source')
        elif mode == 'zero_import':
            config['coupling']['native_fima']['history_path'] = os.path.relpath(tables/'zero.csv',run)
        elif mode == 'uniform_mean':
            config['coupling']['native_fima']['history_path'] = os.path.relpath(tables/'mean.csv',run)
        elif mode in ['uniform_analytic','uniform_constant_reference']:
            config['time'] = [0.0,1.0];config['n_steps'] = [1];config['lhr'] = [0.0,0.0]
            config['time_adaptivity']['enabled'] = False
            material['T_ref'] = al['T_ref'] = 300.0
            material['T_initial'] = al['T_initial'] = 300.0
            config['coupling']['native_fima']['history_path'] = os.path.relpath(tables/'analytic.csv',run)
            if mode == 'uniform_constant_reference':
                config.pop('coupling');material.pop('eigenstrain');material.pop('fima_source');material['swelling'] = 6.0e-4
        else:
            config['coupling']['native_fima']['history_path'] = os.path.relpath(SOURCE,run)
        for name,data in [('input',config),('geometry',geom),('fuel',material),('clad',al),('boundary_conditions',bc)]:
            (run/(name+'.yaml')).write_text(yaml.safe_dump(data,sort_keys=False))
        shutil.copy2(CASE/'mesh.msh',run/'mesh.msh')
        shutil.copy2(CASE/'diagnostics.py',run/'diagnostics.py')
    (CASE/'source_metadata.json').write_text(json.dumps({'source':str(SOURCE),'sha256':history.sha256,'records':len(history.times)*len(history.keys),'coordinate_transform':'r_FE_m=r_OpenMC_cm/100; z_FE_m=(z_OpenMC_cm-10.20)/100','old_regression_length_m':0.356,'new_verification_length_m':0.3556,'axial_stretch':False,'baseline_LHR_W_m':8780,'power_coupled':False,'history_end_h':3500},indent=2)+'\n')


if __name__ == '__main__':
    prepare()
