"""Deterministic object construction and checks only; no OpenMC execution."""
import hashlib
import json
import os
import runpy
from pathlib import Path
import numpy as np
import openmc

initial = runpy.run_path(str(RADIAL_CASE_DIR/'prepare_initial.py'))['prepare_initial'](export=False)
baseline_audit = initial['audit']
protected_paths = [RADIAL_NOTEBOOK, RADIAL_CASE_DIR/'radial_geometry.py', RADIAL_CASE_DIR.parent/'fe101_geometry_static/source_model.py']
protected_hashes = {str(p.relative_to(RADIAL_CASE_DIR.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected_paths}
radial_model = initial['radial_model']
geometry, materials_file = radial_model.geometry, radial_model.materials
settings_file = radial_model.settings
tallies_file = radial_model.tallies
all_cells_before = geometry.get_all_cells()
old_cell_ids = {r['cell_id'] for r in baseline_audit['rings']}
old_material_ids = {r['material_id'] for r in baseline_audit['rings']}
unchanged_cells = {cid:(str(c.region),c.fill_type,getattr(c.fill,'id',None),c.temperature)
                   for cid,c in all_cells_before.items() if cid not in old_cell_ids}
unchanged_materials = {m.id:(m.nuclides.copy(),m.density,m.density_units,m.temperature,list(m._sab),m.depletable)
                       for m in materials_file if m.id not in old_material_ids}
target = next(u for u in geometry.get_all_universes().values() if old_cell_ids <= set(u.cells))
R_FUEL_CM=1.791
radial_edges_cm=R_FUEL_CM*np.sqrt(np.arange(6)/5.0)
axial_edges_cm=np.array([10.200,17.312,24.424,31.536,38.648,45.760])
assert np.allclose(np.diff(axial_edges_cm),7.112,rtol=0,atol=1e-13)
assert np.allclose(radial_edges_cm,[baseline_audit['rings'][0]['r_in_cm']]+[r['r_out_cm'] for r in baseline_audit['rings']],rtol=0,atol=0)
# Reuse exact existing end planes; only add internal axial partition planes.
first_cell=all_cells_before[baseline_audit['rings'][0]['cell_id']]
planes={s.z0:s for s in first_cell.region.get_surfaces().values() if isinstance(s,openmc.ZPlane)}
axial_planes=[planes[axial_edges_cm[0]]]+[openmc.ZPlane(z0=float(z)) for z in axial_edges_cm[1:-1]]+[planes[axial_edges_cm[-1]]]
domain_materials=[]
domain_cells=[]
domains=[]
for ir,ring in enumerate(baseline_audit['rings'],1):
    old_cell=all_cells_before[ring['cell_id']]
    source_material=old_cell.fill
    surfaces=list(old_cell.region.get_surfaces().values())
    cylinders=[s for s in surfaces if isinstance(s,openmc.ZCylinder)]
    outer=next(s for s in cylinders if s.r==radial_edges_cm[ir])
    radial_region=-outer
    if ir>1:
        inner=next(s for s in cylinders if s.r==radial_edges_cm[ir-1])
        radial_region=radial_region & +inner
    target.remove_cell(old_cell)
    materials_file.remove(source_material)
    for iz in range(1,6):
        mat=source_material.clone()
        mat.name=f'Fuel101_B1_r{ir}_z{iz}'
        mat.depletable=True
        volume=float(np.pi*(radial_edges_cm[ir]**2-radial_edges_cm[ir-1]**2)*(axial_edges_cm[iz]-axial_edges_cm[iz-1]))
        mat.volume=volume
        cell=openmc.Cell(name=f'FE101 B1 active r{ir} z{iz}',fill=mat,
                         region=radial_region & +axial_planes[iz-1] & -axial_planes[iz])
        target.add_cell(cell)
        materials_file.append(mat)
        domain_materials.append(mat)
        domain_cells.append(cell)
        assert mat.nuclides==source_material.nuclides
        assert mat.density==6.3 and mat.density_units=='g/cm3'
        assert mat.temperature==source_material.temperature and mat._sab==source_material._sab
        densities=mat.get_nuclide_atom_densities()
        initial_atoms={n:float(d*1e24*volume) for n,d in densities.items()}
        element=lambda n: ''.join(c for c in n if c.isalpha())
        metal_atoms=sum(a for n,a in initial_atoms.items() if element(n) in ('U','Zr'))
        u_mass=sum(mat.get_mass(n) for n in densities if element(n)=='U')/1000
        hzr=sum(d for n,d in densities.items() if element(n)=='H')/sum(d for n,d in densities.items() if element(n)=='Zr')
        domains.append({'radial_index':ir,'axial_index':iz,'material_id':mat.id,'material_name':mat.name,'cell_id':cell.id,
                        'r_in_cm':float(radial_edges_cm[ir-1]),'r_out_cm':float(radial_edges_cm[ir]),
                        'r_in_over_R':float(radial_edges_cm[ir-1]/R_FUEL_CM),'r_out_over_R':float(radial_edges_cm[ir]/R_FUEL_CM),
                        'z_low_cm':float(axial_edges_cm[iz-1]),'z_high_cm':float(axial_edges_cm[iz]),'volume_cm3':volume,
                        'mass_g':mat.get_mass(),'initial_U_mass_kg':u_mass,'initial_U_Zr_atoms':metal_atoms,
                        'initial_atoms_by_nuclide':initial_atoms,'H_Zr_atom_ratio':hzr,'temperature_K':mat.temperature})
geometry.determine_paths(instances_only=True)
assert len({m.id for m in domain_materials})==25
assert len([m for m in materials_file if m.depletable])==25
assert {m.id for m in geometry.get_all_materials().values() if m.depletable}=={m.id for m in domain_materials}
assert all(m.num_instances==c.num_instances==1 for m,c in zip(domain_materials,domain_cells))
for d,m in zip(domains,domain_materials): d['instances']=m.num_instances
for cid,snapshot in unchanged_cells.items():
    c=geometry.get_all_cells()[cid]
    assert (str(c.region),c.fill_type,getattr(c.fill,'id',None),c.temperature)==snapshot
for mid,snapshot in unchanged_materials.items():
    m=next(m for m in materials_file if m.id==mid)
    assert (m.nuclides,m.density,m.density_units,m.temperature,list(m._sab),m.depletable)==snapshot
    assert not m.depletable
original_volume=baseline_audit['original_volume_cm3']
assert np.allclose([d['volume_cm3'] for d in domains],original_volume/25,rtol=1e-14)
assert np.isclose(sum(d['volume_cm3'] for d in domains),original_volume,rtol=1e-14)
assert np.isclose(original_volume,358.346194419,rtol=0,atol=5e-10)
original_mass=baseline_audit['original_mass_g']
original_u_mass=sum(r['initial_U_mass_kg'] for r in baseline_audit['rings'])
assert np.isclose(sum(d['mass_g'] for d in domains),original_mass,rtol=1e-14)
assert np.isclose(sum(d['initial_U_mass_kg'] for d in domains),original_u_mass,rtol=1e-14)
nuclide_checks={}
for nuc in baseline_audit['rings'][0]['initial_atoms_by_nuclide']:
    expected=sum(r['initial_atoms_by_nuclide'][nuc] for r in baseline_audit['rings'])
    actual=sum(d['initial_atoms_by_nuclide'][nuc] for d in domains)
    assert np.isclose(actual,expected,rtol=1e-14)
    nuclide_checks[nuc]={'original_atoms':expected,'reconstructed_atoms':actual,'relative_error':(actual-expected)/expected}
radial_volumes=[sum(d['volume_cm3'] for d in domains if d['radial_index']==ir) for ir in range(1,6)]
axial_volumes=[sum(d['volume_cm3'] for d in domains if d['axial_index']==iz) for iz in range(1,6)]
assert np.allclose(radial_volumes,[r['volume_cm3'] for r in baseline_audit['rings']],rtol=1e-14)
assert np.allclose(axial_volumes,original_volume/5,rtol=1e-14)
origin=np.asarray(baseline_audit['B1_origin_cm'])
for d,mat in zip(domains,domain_materials):
    rad=np.sqrt((d['r_in_cm']**2+d['r_out_cm']**2)/2)
    for phi in (0.,0.7,2.,4.):
        for zz in (d['z_low_cm']+1e-6,(d['z_low_cm']+d['z_high_cm'])/2,d['z_high_cm']-1e-6):
            assert geometry.find(tuple(origin+[rad*np.cos(phi),rad*np.sin(phi),zz]))[-1].fill is mat
for zz in (10.200001,17.312,27.98,38.648,45.759999):
    for rad,expected_name in [(1.791001,'He4'),(1.803999,'He4'),(1.804001,'Al'),(1.879999,'Al')]:
        mat=geometry.find(tuple(origin+[rad,0,zz]))[-1].fill
        if expected_name=='He4': assert mat.nuclides[0].name=='He4' and mat.density==1e-10
        else: assert mat.id==geometry.find(tuple(origin+[1.85,0,27.98]))[-1].fill.id
ordinary_fuel=next(m for m in materials_file if m.name=='Fuel 101')
assert ordinary_fuel.num_instances==60 and not ordinary_fuel.depletable
# Replace only radial material diagnostics; other original diagnostics remain.
tallies_file=openmc.Tallies([t for t in tallies_file if not t.name.startswith('FE101 B1 ring')])
domain_filter=openmc.MaterialFilter(domain_materials, filter_id=max(f.id for t in radial_model.tallies for f in t.filters)+1)
for name,score in [('FE101 B1 domain total fission','fission'),('FE101 B1 domain energy deposition','heating-local')]:
    t=openmc.Tally(name=name)
    t.filters=[domain_filter]
    t.scores=[score]
    if score=='fission': t.estimator='tracklength'
    tallies_file.append(t)
radial_model.tallies=tallies_file
assert (settings_file.particles,settings_file.batches,settings_file.inactive)==(5000,30,10)
openmc.config['cross_sections']=baseline_audit['cross_sections']
openmc.config['chain_file']=baseline_audit['chain_file']
audit={'geometry_baseline':'physical B','fuel_outer_cm':1.791,'clad_inner_cm':1.804,'clad_outer_cm':1.880,
       'gap_radial_um':130.,'radial_edges_cm':radial_edges_cm.tolist(),'axial_edges_cm':axial_edges_cm.tolist(),
       'B1_origin_cm':origin.tolist(),'domains':domains,'depletable_materials':25,'ordinary_FE101_instances':60,
       'original_volume_cm3':original_volume,'reconstructed_volume_cm3':sum(d['volume_cm3'] for d in domains),
       'original_mass_g':original_mass,'reconstructed_mass_g':sum(d['mass_g'] for d in domains),
       'original_U_mass_kg':original_u_mass,'reconstructed_U_mass_kg':sum(d['initial_U_mass_kg'] for d in domains),
       'nuclide_conservation':nuclide_checks,'radial_group_volumes_cm3':radial_volumes,'axial_group_volumes_cm3':axial_volumes,
       'cross_sections':baseline_audit['cross_sections'],'chain_file':baseline_audit['chain_file'],
       'settings':baseline_audit['settings'],'normalization_mode':'energy-deposition','diff_burnable_mats':False,
       'write_rates':True,'final_step':True,'power_W':250000.,'steps':20,'step_hours':175.,
       'protected_sha256':protected_hashes,'checks_passed':True,'transport_executed':False,'depletion_executed':False,
       'radial_benchmark':os.path.relpath(RADIAL_CASE_DIR/'reference/radial_history.json',PREPARED_DIR),
       'static_benchmark':os.path.relpath(RADIAL_CASE_DIR.parent/'fe101_axial_static/reference/static_5x5.json',PREPARED_DIR)}
assert (PREPARED_DIR/audit['radial_benchmark']).resolve().is_file() and (PREPARED_DIR/audit['static_benchmark']).resolve().is_file()
if EXPORT_PREPARED_INPUTS:
    PREPARED_DIR.mkdir(parents=True,exist_ok=True)
    geometry.export_to_xml(PREPARED_DIR/'geometry.xml')
    materials_file.export_to_xml(PREPARED_DIR/'materials.xml')
    settings_file.export_to_xml(PREPARED_DIR/'settings.xml')
    tallies_file.export_to_xml(PREPARED_DIR/'tallies.xml')
    (PREPARED_DIR/'preparation_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    import pandas as pd
    pd.DataFrame(domains).drop(columns='initial_atoms_by_nuclide').to_csv(PREPARED_DIR/'domains.csv',index=False)
for p,h in protected_hashes.items(): assert hashlib.sha256((RADIAL_CASE_DIR.parent/p).read_bytes()).hexdigest()==h
print('PASS: 25 unique materials/cells, one instance each; all conservation and geometry checks passed')
print(f'Volume: {original_volume:.12f} cm3; domain: {original_volume/25:.12f} cm3')
print(f'Fuel mass: {original_mass:.12f} g; U mass: {original_u_mass:.12f} kg')
