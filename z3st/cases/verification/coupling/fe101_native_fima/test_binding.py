"""FE coefficient/state checks without any thermal/mechanical solve."""
import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import yaml
import dolfinx

from z3st.core.spine import Spine
from z3st.materials.fuel_swelling import uzrh_fission_product_swelling
from z3st.utils.writer import OutputWriter

CASE=Path(__file__).resolve().parent


def main():
    # Mirror the CLI's compatibility path for legacy material-card imports.
    sys.path.insert(0,str(CASE.parents[3]))
    os.chdir(CASE/'runs/native')
    load=lambda p:yaml.safe_load(Path(p).read_text())
    inp=load('input.yaml')
    p=Spine(inp,'mesh.msh',load('geometry.yaml'))
    p.load_materials(**{name:load(path) for name,path in inp['materials'].items()})
    p.parameters(8780)
    p.initialize_fields()
    native=p._fima_coupling
    q=p.q_third.x.array.copy()
    p.burnup.x.array[:]=1000.0
    bu=p.burnup.x.array.copy()
    swelling=uzrh_fission_product_swelling(p.T,p.materials['fuel'],model=p)
    tensor_space=dolfinx.fem.functionspace(p.mesh,('DG',0,(3,3)))
    sample=dolfinx.fem.Function(tensor_space)
    expression=dolfinx.fem.Expression(swelling,tensor_space.element.interpolation_points)
    sample.interpolate(expression)
    np.testing.assert_array_equal(sample.x.array,0)
    p.set_coupling_time(630000.0*0.37)
    np.testing.assert_allclose(p.fima_native.x.array[native.dofs],native.transfer.fima_at(630000.0*0.37),rtol=0,atol=0)
    snapshot=p.snapshot_state()
    saved=p.fima_native.x.array.copy()
    p.set_coupling_time(12600000.0)
    p.restore_state(snapshot)
    np.testing.assert_array_equal(p.fima_native.x.array,saved)
    assert native.time_s==630000.0*0.37
    np.testing.assert_array_equal(p.q_third.x.array,q)
    np.testing.assert_array_equal(p.burnup.x.array,bu)
    sample.interpolate(expression)
    tensor_dofs=np.array([tensor_space.dofmap.cell_dofs(int(c))[0] for c in native.cells])
    np.testing.assert_allclose(sample.x.array.reshape(-1,3,3)[tensor_dofs],saved[native.dofs,None,None]*np.eye(3),rtol=1e-14,atol=0)
    clad_dofs=[p.Q.dofmap.cell_dofs(int(c))[0] for c in p.cell_tags.find(p.label_map['clad'])]
    np.testing.assert_array_equal(p.fima_native.x.array[clad_dofs],0)
    with np.testing.assert_raises(ValueError):p.set_coupling_time(12600000.0+1)
    assert native.time_s==630000.0*0.37
    class MissingNative:
        burnup=p.burnup
    with np.testing.assert_raises(ValueError):
        uzrh_fission_product_swelling(p.T,{'fima_source':'native_openmc'},model=MissingNative())
    # The disabled adapter remains optional and the legacy law still operates.
    legacy=uzrh_fission_product_swelling(p.T,{k:v for k,v in p.materials['fuel'].items() if k!='fima_source'},model=p)
    sample.interpolate(dolfinx.fem.Expression(legacy,tensor_space.element.interpolation_points))
    assert (sample.x.array.reshape(-1,3,3)[tensor_dofs,0,0]>0).all()
    # Exercise the alternative DG0 XDMF export route without a solve.
    p.get_results()
    export=CASE/'interface_export_test'
    with OutputWriter(p,output_format='xdmf',output_dir=str(export),n_steps=1) as writer:
        writer.write(t=native.time_s,step=0)
    tree=ET.parse(export/'fields.xdmf')
    attributes={a.attrib['Name']:a.attrib.get('Center') for a in tree.iter('Attribute')}
    assert attributes['FIMA_native_OpenMC']=='Cell'
    assert attributes['Swelling_eigenstrain_native']=='Cell'
    report={'PASS':True,'solve_executed':False,'rollback_time_and_coefficients':True,
            'native_source_independent_of_internal_BU':True,'thermal_source_unchanged':True,
            'isotropic_diagonal_equals_imported_FIMA':True,'offdiagonal_zero':True,
            'cladding_imported_FIMA_zero':True,'missing_native_field_raises':True,
            'legacy_burnup_mode_preserved':True,'out_of_range_rejected_without_state_change':True,
            'XDMF_native_fields_exported_as_cell_data':True,
            'off_grid_time_s':native.time_s,'MPI_processes':p.mesh.comm.size}
    (CASE/'binding_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
