"""Synthetic data only: validate 525-row outputs, integrals and aggregations."""
import json
import runpy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import pandas as pd
import openmc
import openmc.deplete
from postprocess import postprocess_2d
from benchmarks import aggregate_history, compare_benchmarks
ROOT=Path(__file__).resolve().parent
fixture=TemporaryDirectory(prefix='fe101_5x5_fixture_')
ns=runpy.run_path(str(ROOT/'validate_preparation.py'))['validate'](fixture.name)
audit_path=Path(fixture.name)/'prepared_inputs/preparation_audit.json'
meta=json.loads(audit_path.read_text())
domains=meta['domains']
ids=[d['material_id'] for d in domains]
times=np.arange(21)*175*3600.
rates={str(d['material_id']):float((d['radial_index']+d['axial_index'])*1e12) for d in domains}
heat=np.arange(1,26,dtype=float)*1e6
global_heat=1e9
class FakeResults:
    def __init__(self,*args): pass
    def __getitem__(self,index): return SimpleNamespace(index_mat={str(i):j for j,i in enumerate(ids)},rates=SimpleNamespace(index_nuc={'U235':0,'U238':1,'Pu239':2}),index_nuc={'U235':0,'U238':1,'Pu239':2,'Cs137':3})
    def get_times(self,**kwargs): return times
    def get_atoms(self,mid,nuc,**kwargs):
        d=domains[ids.index(int(mid))]
        initial=d['initial_atoms_by_nuclide'].get(nuc,0.)
        return times,np.full(21,initial)
    def get_reaction_rate(self,mid,nuc,reaction):
        share={'U235':.8,'U238':.15,'Pu239':.05}[nuc]
        return times,np.full(21,rates[mid]*share)
class FakeStatePoint:
    def __init__(self,*args): pass
    def __enter__(self): return self
    def __exit__(self,*args): pass
    def get_tally(self,name):
        if name=='factor-for-normalization': return SimpleNamespace(mean=np.array([global_heat]))
        if name=='FE101 B1 domain energy deposition': values=heat; score='heating-local'
        else:
            source=250000./(global_heat*openmc.data.JOULE_PER_EV)
            values=np.array([rates[str(mid)]/source for mid in ids]); score='fission'
        return SimpleNamespace(mean=values[:,None,None],scores=[score],filters=[SimpleNamespace(bins=ids)])
chain=SimpleNamespace(nuclides=[SimpleNamespace(name=n,reactions=[SimpleNamespace(type='fission')]) for n in ('U235','U238','Pu239')])
with TemporaryDirectory(prefix='fe101_5x5_synthetic_') as tmp:
    with patch.object(openmc,'StatePoint',FakeStatePoint),patch.object(openmc.deplete,'Results',FakeResults),patch.object(openmc.deplete.Chain,'from_xml',return_value=chain):
        frame=postprocess_2d(tmp,audit_path)
    final=frame[frame.time_h==3500].sort_values(['radial_index','axial_index'])
    expected=np.array([rates[str(d['material_id'])]*times[-1]/d['initial_U_Zr_atoms'] for d in domains])
    assert np.allclose(final.FIMA,expected,rtol=1e-14,atol=0)
    power=250000.*heat/global_heat
    expected_bu=power*times[-1]/(8.64e10*np.array([d['initial_U_mass_kg'] for d in domains]))
    assert np.allclose(final.BU_MWd_kgU,expected_bu,rtol=1e-14,atol=0)
    assert np.allclose(frame.fission_rate_per_s,frame.tally_fission_rate_per_s,rtol=1e-14)
    whole=aggregate_history(frame,['time_h'])
    rings=aggregate_history(frame,['time_h','radial_index'])
    axial=aggregate_history(frame,['time_h','axial_index'])
    assert len(whole)==21 and len(rings)==105 and len(axial)==105
    transfer=pd.read_csv(Path(tmp)/'postprocessing/zest_FIMA_2D.csv')
    assert set(['time_h','radial_index','axial_index','r_in_cm','r_out_cm','z_low_cm','z_high_cm','FIMA'])<=set(transfer.columns)
    assert np.allclose(transfer[transfer.time_h==0].FIMA,0)
    benchmark=compare_benchmarks(tmp,audit_path)
    assert benchmark['radial_comparison_rows']==105 and benchmark['static_comparison_rows']==25
fixture.cleanup()
print('PASS: SYNTHETIC ONLY — 525 rows, all-nuclide fission sum, fixed-denominator FIMA, deposited-energy BU, transfer schema, weighted FE/radial/axial aggregation')
