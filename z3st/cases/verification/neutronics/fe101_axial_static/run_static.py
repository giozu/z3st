"""Static tally-only axial check from unchanged initial radial-case sources."""
import hashlib
import json
import os
import runpy
import subprocess
import time
from pathlib import Path
import numpy as np
import openmc

ROOT = Path(__file__).resolve().parent
RADIAL = ROOT.parent/'fe101_radial_depletion'
RUN = ROOT/'run_preliminary'
if RUN.exists() and any(RUN.iterdir()):
    raise FileExistsError(f'Nonempty run directory: {RUN}')
RUN.mkdir(exist_ok=True)
initial = runpy.run_path(str(RADIAL/'prepare_initial.py'))['prepare_initial'](export=False)
audit = initial['audit']
model = initial['radial_model']
protected_paths = [RADIAL/'triga_single_FE101_B1_1965_5RINGS_CLEAN.ipynb', ROOT.parent/'fe101_geometry_static/source_model.py']
hash_file = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
protected = {str(p.relative_to(ROOT.parent)):hash_file(p) for p in protected_paths}
for name, obj in [('geometry',model.geometry),('materials',model.materials),('settings',model.settings)]:
    obj.export_to_xml(RUN/(name+'.xml'))
baseline_hashes = {name:hash_file(RUN/(name+'.xml')) for name in ['geometry','materials','settings']}
assert (model.settings.particles, model.settings.batches, model.settings.inactive)==(5000,30,10)
z = np.array([10.200,17.312,24.424,31.536,38.648,45.760])
r = 1.791*np.sqrt(np.arange(6)/5)
origin = audit['B1_origin_cm']
assert origin[2] == 0.0
assert np.allclose(np.diff(z),7.112,rtol=0,atol=1e-13)
cell_ids = [a['cell_id'] for a in audit['rings']]
material_ids = [a['material_id'] for a in audit['rings']]
for ring in audit['rings']:
    mat = model.geometry.get_all_materials()[ring['material_id']]
    atoms = mat.get_nuclide_atom_densities()
    for nuc, n in ring['initial_atoms_by_nuclide'].items():
        assert np.isclose(atoms[nuc]*1e24*ring['volume_cm3'], n, rtol=1e-12)
filter_target = openmc.CellFilter(cell_ids)
ax_mesh = openmc.CylindricalMesh(r_grid=[0.,1.791], z_grid=z, origin=origin, name='B1 axial five zones')
rz_mesh = openmc.CylindricalMesh(r_grid=r, z_grid=z, origin=origin, name='B1 radial axial 5x5 tally mesh')
for name, mesh in [('B1 static integrated',None),('B1 static axial',ax_mesh),('B1 static radial axial',rz_mesh)]:
    tally = openmc.Tally(name=name)
    tally.filters = [filter_target] + ([] if mesh is None else [openmc.MeshFilter(mesh)])
    tally.scores = ['fission','heating-local']
    tally.estimator = 'tracklength'
    model.tallies.append(tally)
model.tallies.export_to_xml(RUN/'tallies.xml')
volumes = np.pi*np.diff(r*r)[:,None]*np.diff(z)[None,:]
assert np.allclose(volumes, audit['original_volume_cm3']/25, rtol=1e-14)
assert np.isclose(volumes.sum(),audit['original_volume_cm3'],rtol=1e-14)
for name in ['geometry.xml','materials.xml','settings.xml']:
    assert baseline_hashes[Path(name).stem]==hash_file(RUN/name)
manifest = {'baseline':'physical B, unchanged initial radial-case inputs','radial_edges_cm':r.tolist(), 'axial_edges_cm':z.tolist(), 'origin_cm':origin,'bin_volumes_cm3':volumes.tolist(),'FE_volume_cm3':float(volumes.sum()),'cell_ids':cell_ids,'material_ids':material_ids,'power_W_for_density':250000.,'settings':{'particles':5000,'batches':30,'inactive':10},'protected_sha256':protected,'no_depletion_operator':True,'new_materials':0}
(ROOT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
start=time.perf_counter()
with (ROOT/'transport.log').open('w') as log:
    env = dict(os.environ, OPENMC_CROSS_SECTIONS=audit['cross_sections'])
    result=subprocess.run(['openmc','-s','8'],cwd=RUN,stdout=log,stderr=subprocess.STDOUT,env=env)
status={'exit_code':result.returncode,'transport_runtime_s':time.perf_counter()-start}
for p,h in protected.items():
    assert hash_file(ROOT.parent/p)==h
status['protected_inputs_unchanged']=True
(ROOT/'execution_status.json').write_text(json.dumps(status,indent=2)+'\n')
print(status,flush=True)
raise SystemExit(result.returncode)
