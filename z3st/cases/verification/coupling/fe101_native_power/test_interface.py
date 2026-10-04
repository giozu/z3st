"""Saved tally provenance, conservative thermal transfer, invalid-data tests."""
import csv,json,runpy,tempfile,unittest
from pathlib import Path
import numpy as np
import h5py
from prepare import CASE,SOURCE,HIGH,PHASE1
from z3st.coupling.openmc.depletion_fields import HeatingHistory,DepletionHistory
from z3st.coupling.openmc.conservative_transfer import ConservativeTransfer


class PowerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Import only the pure mesh parser; never execute/write phase-1 tests.
        parser=runpy.run_path(str(PHASE1/'test_interface.py'))['mesh_fuel_bounds']
        cls.h=HeatingHistory(SOURCE);cls.b=parser(CASE/'mesh.msh')
        cls.m=ConservativeTransfer(cls.h,cls.b,10.20)

    def test_saved_heating_provenance(self):
        metadata=json.loads((HIGH/'prepared_inputs/preparation_audit.json').read_text())
        ids={(d['radial_index'],d['axial_index']):d['material_id'] for d in metadata['domains']}
        errors=[]
        for k,t in enumerate(self.h.times):
            with h5py.File(HIGH/f'run_3500h_ED/openmc_simulation_n{k}.h5') as f:
                def tally(name):
                    g=next(f['tallies'][key] for key in f['tallies'] if key.startswith('tally ') and 'name' in f['tallies'][key] and f['tallies'][key]['name'][()].decode()==name)
                    assert g['score_bins'][()].tolist()==[b'heating-local']
                    return g,g['results'][()][:,0,0]/int(g['n_realizations'][()])
                global_t,global_h=tally('factor-for-normalization')
                local,local_h=tally('FE101 B1 domain energy deposition')
                bins=list(f[f"tallies/filters/filter {int(local['filters'][()][0])}/bins"][()])
                expected=np.array([250000*local_h[bins.index(ids[key])]/global_h[0] for key in self.h.keys])
                np.testing.assert_allclose(expected,self.h.power_W[k],rtol=2e-14,atol=1e-13)
                errors.append(float(np.max(abs(expected-self.h.power_W[k]))))
        (CASE/'source_heating_verification.json').write_text(json.dumps({'PASS':True,'statepoints_verified':21,'score':'heating-local','power_column':'deposited_power_W','native_units':'eV/source','normalization':'250000 W * H_domain/H_core','max_saved_power_vs_HDF5_abs_error_W':max(errors)},indent=2)+'\n')

    def test_conservation_all_times_and_off_grid(self):
        h,m=self.h,self.m;m.validate_complete_coverage()
        rows=[]
        for t in np.sort(np.r_[h.times,.37*h.times[:-1]+.63*h.times[1:]]):
            q=m.qdot_at(t);p=h.power_at(t).sum();after=np.dot(q,m.volumes_m3)
            assert np.isfinite(q).all() and (q>=0).all()
            self.assertLess(abs(after/p-1),1e-12)
            rows.append({'time_h':float(t/3600),'source_power_W':float(p),'mapped_power_W':float(after),'absolute_error_W':float(after-p),'relative_error':float((after-p)/p)})
        with (CASE/'power_conservation_interface.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)

    def test_linear_interpolation_no_extrapolation(self):
        h=self.h
        for lo,hi in zip(h.times[:-1],h.times[1:]):
            np.testing.assert_allclose(h.power_at(lo+.37*(hi-lo)),.63*h.power_at(lo)+.37*h.power_at(hi),rtol=2e-14)
        for t in [-1,h.times[-1]+1,np.nan]:
            with self.assertRaises(ValueError):h.power_at(t)

    def test_null_uniform_and_fima_compatibility(self):
        for mode in ['zero','known_uniform','same_power_uniform']:
            h=HeatingHistory(CASE/'synthetic_histories'/(mode+'.csv'));m=ConservativeTransfer(h,self.b,10.20)
            for t in h.times:
                q=m.qdot_at(t)
                np.testing.assert_allclose(q,h.qdot_at(t)[0],rtol=3e-13,atol=1e-9)
        legacy=DepletionHistory(SOURCE)
        for t in legacy.times:
            np.testing.assert_array_equal(self.h.fima_at(t),legacy.fima_at(t))
            np.testing.assert_allclose(self.m.fima_at(t),ConservativeTransfer(legacy,self.b,10.20).fima_at(t),rtol=0,atol=0)

    def test_bad_heating_rejected(self):
        with SOURCE.open() as f:rows=list(csv.DictReader(f))
        for value in ['-1','nan','inf']:
            row=dict(rows[0]);row['deposited_power_W']=value
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/'bad.csv'
                with p.open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=row.keys());w.writeheader();w.writerows([row]+rows[1:])
                with self.assertRaises(ValueError):HeatingHistory(p)


if __name__=='__main__':unittest.main(verbosity=2)
