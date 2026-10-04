"""Interface tests run before any coupled solve; stdlib unittest + numpy."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from z3st.coupling.openmc.depletion_fields import DepletionHistory
from z3st.coupling.openmc.conservative_transfer import ConservativeTransfer, cylindrical_volumes
from prepare import CASE, SOURCE


def mesh_fuel_bounds(path):
    """Read actual ASCII Gmsh4 quads, without generating/loading a solver."""
    text = Path(path).read_text()
    tokens = iter(text.split('$Nodes\n')[1].split('$EndNodes')[0].split())
    blocks, _, _, _ = [int(next(tokens)) for _ in range(4)]
    nodes = {}
    for _ in range(blocks):
        dim,tag,param,n = [int(next(tokens)) for _ in range(4)]
        ids = [int(next(tokens)) for _ in range(n)]
        for nid in ids:
            nodes[nid] = [float(next(tokens)) for _ in range(3)]
            for j in range(dim if param else 0):
                next(tokens)
    lines = text.split('$Elements\n')[1].split('$EndElements')[0].splitlines()
    nblocks = int(lines[0].split()[0]);i = 1;bounds = []
    for _ in range(nblocks):
        dim,tag,kind,n = map(int,lines[i].split());i += 1
        for line in lines[i:i+n]:
            if dim == 2 and tag == 1:
                assert kind == 3
                xyz = np.array([nodes[j] for j in map(int,line.split()[1:])])
                bounds.append([xyz[:,0].min(),xyz[:,0].max(),xyz[:,1].min(),xyz[:,1].max()])
        i += n
    return np.array(bounds)


class InterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.h = DepletionHistory(SOURCE)
        cls.bounds = mesh_fuel_bounds(CASE/'mesh.msh')
        cls.mapping = ConservativeTransfer(cls.h,cls.bounds,10.20)

    def test_identity_by_physical_overlap(self):
        b = self.h.bounds_cm.copy()/100;b[:,2:] -= 0.102
        identity = ConservativeTransfer(self.h,b,10.20)
        self.assertLess(identity.validate_complete_coverage(),1e-12)
        np.testing.assert_allclose(identity.fractions,np.eye(100),rtol=0,atol=3e-14)
        for t in self.h.times:
            np.testing.assert_allclose(identity.fima_at(t),self.h.fima_at(t),rtol=2e-13,atol=1e-18)

    def test_conservation_and_history(self):
        m,h = self.mapping,self.h
        self.assertEqual(len(self.bounds),720)
        self.assertLess(m.validate_complete_coverage(),1e-12)
        self.assertLess(abs(m.volumes_m3.sum()/(h.volumes_cm3.sum()*1e-6)-1),1e-13)
        self.assertLess(abs(m.initial_atoms.sum()/h.initial_atoms.sum()-1),1e-13)
        np.testing.assert_array_equal(m.fima_at(0),0)
        times = np.sort(np.r_[h.times,(h.times[1:]+h.times[:-1])/2])
        values = np.array([m.fima_at(t) for t in times])
        self.assertTrue((values >= 0).all())
        self.assertTrue((np.diff(values,axis=0) >= 0).all())
        for t,fi in zip(times,values):
            target = h.fissions_at(t).sum()
            np.testing.assert_allclose(np.dot(fi,m.initial_atoms),target,rtol=3e-13,atol=0)
            mean = np.average(fi,weights=m.volumes_m3)
            expected = np.average(h.fima_at(t),weights=h.volumes_cm3)
            np.testing.assert_allclose(mean,expected,rtol=3e-13,atol=1e-18)

    def test_interpolation_and_reject_extrapolation(self):
        h = self.h
        for t0,t1 in zip(h.times[:-1],h.times[1:]):
            t = t0+0.37*(t1-t0)
            np.testing.assert_allclose(h.fima_at(t),0.63*h.fima_at(t0)+0.37*h.fima_at(t1),rtol=3e-14,atol=1e-18)
        for t in [-1,h.times[-1]+1,np.nan]:
            with self.assertRaises(ValueError):h.fima_at(t)

    def test_uniform_and_zero(self):
        for kind in ['mean','analytic','zero']:
            h = DepletionHistory(CASE/'synthetic_histories'/(kind+'.csv'))
            m = ConservativeTransfer(h,self.bounds,10.20)
            for t in h.times:
                np.testing.assert_allclose(m.fima_at(t),h.fima_at(t)[0],rtol=1e-13,atol=1e-18)

    def test_reject_mismatch_and_duplicate_coverage(self):
        bad = self.bounds.copy();bad[:,2:] *= .356/.3556
        with self.assertRaises(ValueError):ConservativeTransfer(self.h,bad,10.20)
        duplicate = ConservativeTransfer(self.h,np.vstack([self.bounds,self.bounds[:1]]),10.20)
        with self.assertRaises(ValueError):duplicate.validate_complete_coverage()

    def test_reader_keys_and_cumulative_consistency(self):
        with SOURCE.open() as stream:rows=list(csv.DictReader(stream))
        def check_bad(values):
            with tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'bad.csv'
                with path.open('w',newline='') as stream:
                    w=csv.DictWriter(stream,fieldnames=rows[0].keys());w.writeheader();w.writerows(values)
                with self.assertRaises(ValueError):DepletionHistory(path)
        check_bad(rows+[rows[0]])
        bad=[dict(row) for row in rows];bad[-1]['FIMA']='-1';check_bad(bad)
        bad=[dict(row) for row in rows];bad[-1]['cumulative_fissions']='1';check_bad(bad)
        # Row order has no physical meaning: reverse all records and recover.
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'reordered.csv'
            with path.open('w',newline='') as stream:
                w=csv.DictWriter(stream,fieldnames=rows[0].keys());w.writeheader();w.writerows(reversed(rows))
            np.testing.assert_array_equal(DepletionHistory(path).fima,self.h.fima)

    def test_native_spatial_trends(self):
        m,h=self.mapping,self.h
        native=h.fima[-1];mapped=m.fima_at(h.times[-1])
        # Rebin the mapped cell field only for comparison of spatial profiles.
        reconstructed=(m.overlap_m3.T@mapped)/(h.volumes_cm3*1e-6)
        radial=native.reshape(10,10).mean(axis=1)
        radial_mapped=reconstructed.reshape(10,10).mean(axis=1)
        axial=native.reshape(10,10).mean(axis=0)
        axial_mapped=reconstructed.reshape(10,10).mean(axis=0)
        self.assertTrue((np.diff(radial)>0).all())
        self.assertTrue((np.diff(radial_mapped)>0).all())
        np.testing.assert_allclose(axial_mapped,axial,rtol=2e-13,atol=1e-18)
        index=int(np.argmax(mapped));b=self.bounds[index]
        source_index=int(np.argmax(native));source=h.bounds_cm[source_index]/100
        source[2:] -= .102
        centre=np.array([(b[0]+b[1])/2,(b[2]+b[3])/2])
        self.assertTrue(source[0] <= centre[0] <= source[1] and source[2] <= centre[1] <= source[3])
        volume_error=abs(m.volumes_m3.sum()/(h.volumes_cm3.sum()*1e-6)-1)
        atom_error=abs(m.initial_atoms.sum()/h.initial_atoms.sum()-1)
        mean_error=abs(np.average(mapped,weights=m.volumes_m3)/np.average(native,weights=h.volumes_cm3)-1)
        metrics={'A_B_pass':True,'source_records':2100,'destination_fuel_cells':720,
                 'volume_relative_error':float(volume_error),'initial_atoms_relative_error':float(atom_error),
                 'final_FIMA_mean_relative_error':float(mean_error),
                 'max_relative_fission_sum_error':float(max(abs(np.dot(m.fima_at(t),m.initial_atoms)/h.fissions_at(t).sum()-1) for t in h.times[1:])),
                 'source_final_max_FIMA':float(native.max()),'mapped_final_max_FIMA':float(mapped.max()),
                 'peak_change_percent':float(100*(mapped.max()/native.max()-1)),
                 'native_max_domain_r_z':list(h.keys[source_index]),'mapped_max_cell_bounds_m':b.tolist(),
                 'radial_profile_native':radial.tolist(),'radial_profile_reconstructed':radial_mapped.tolist(),
                 'axial_profile_native':axial.tolist(),'axial_profile_reconstructed':axial_mapped.tolist(),
                 'max_radial_profile_remap_difference_percent':float(100*np.max(abs(radial_mapped/radial-1))),
                 'max_bin_remap_difference_percent':float(100*np.max(abs(reconstructed/native-1))),
                 'interpolation':'linear cumulative fissions; tested midpoint and off-grid 37% locations; no extrapolation'}
        (CASE/'interface_verification.json').write_text(json.dumps(metrics,indent=2)+'\n')


if __name__=='__main__':
    unittest.main(verbosity=2)
