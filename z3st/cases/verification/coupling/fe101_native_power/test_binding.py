"""Native power ownership, projection and two-field rollback; no solve."""
import json,os,sys
from pathlib import Path
import numpy as np
import yaml
from prepare import CASE,ROOT
from z3st.core.spine import Spine
from z3st.utils.writer import OutputWriter

def main():
    sys.path.insert(0,str(ROOT/'z3st'));os.chdir(CASE/'runs/full_native')
    load=lambda path:yaml.safe_load(Path(path).read_text());inp=load('input.yaml')
    p=Spine(inp,'mesh.msh',load('geometry.yaml'))
    p.load_materials(**{name:load(path) for name,path in inp['materials'].items()})
    p.parameters(8780);p.initialize_fields()
    power=p._power_coupling
    assert p.q_third is p.qdot_native and p.q_third is not p.q_third_baseline
    p.q_third_baseline.x.array[:]=1e30
    t=630000*.37;p.set_coupling_time(t)
    q=p.qdot_native.x.array.copy();fi=p.fima_native.x.array.copy()
    p.set_power();np.testing.assert_array_equal(q,p.qdot_native.x.array)
    p0=float(power.history.power_at(t).sum())
    integral=p.material_mean(p.qdot_native,'fuel',integral=True)
    projected=p.material_mean(power.burnup_source_nodal,'fuel',integral=True)
    np.testing.assert_allclose([integral,projected],p0,rtol=3e-13)
    assert (power.burnup_source_nodal.x.array>=0).all()
    snap=p.snapshot_state();p.set_coupling_time(12600000);p.restore_state(snap)
    np.testing.assert_array_equal(p.fima_native.x.array,fi)
    np.testing.assert_array_equal(p.qdot_native.x.array,q)
    assert p._fima_coupling.time_s==power.time_s==t
    with np.testing.assert_raises(ValueError):p.set_coupling_time(12600001)
    assert p._fima_coupling.time_s==power.time_s==t
    p.get_results()
    with OutputWriter(p,output_format='xdmf',output_dir=str(CASE/'binding_export_test'),n_steps=1) as writer:writer.write(t=t,step=0)
    import xml.etree.ElementTree as ET
    attrs={a.attrib['Name']:a.attrib.get('Center') for a in ET.parse(CASE/'binding_export_test/fields.xdmf').iter('Attribute')}
    assert attrs['qdot_native_OpenMC_W_m3']=='Cell'
    (CASE/'binding_verification.json').write_text(json.dumps({'PASS':True,'no_solve':True,'absolute_power_ownership':True,'baseline_not_used':True,'simultaneous_time':True,'two_field_rollback':True,'no_extrapolation':True,'positive_conservative_nodal_BU_projection':True,'projection_power_error_W':float(projected-p0),'XDMF_DG0_export':True},indent=2)+'\n')

if __name__=='__main__':main()
