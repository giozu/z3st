"""Cylindrical overlap transfer with independent opt-in FIMA/heating adapters.

Pure-numpy geometry/history tests do not import FEniCS. Legacy FIMA-only
behavior remains unchanged when native_power is disabled.
"""
import numpy as np

from .depletion_fields import DepletionHistory, HeatingHistory


def cylindrical_volumes(bounds):
    ri, ro, zl, zh = np.asarray(bounds, dtype=float).T
    return np.pi*(ro**2-ri**2)*(zh-zl)


class ConservativeTransfer:
    """Destination rectangular r-z cells in metres, source physical bins in cm.

    Explicit transform: r_FE=r_OpenMC/100,
    z_FE=(z_OpenMC-axial_origin_cm)/100. No axial stretching/clamping.
    Transfer fissions and initial atoms separately, then reconstruct FIMA.
    """

    def __init__(self, history, destination_bounds_m, axial_origin_cm):
        self.history = history
        self.destination_bounds_m = np.asarray(destination_bounds_m, dtype=float)
        b = history.bounds_cm.copy()/100
        b[:, 2:] -= float(axial_origin_cm)/100
        d = self.destination_bounds_m
        if d.ndim != 2 or d.shape[1] != 4 or not np.isfinite(d).all():
            raise ValueError('Destination bounds must be finite (n,4) r-z rectangles')
        self.volumes_m3 = cylindrical_volumes(d)
        if (d[:, 0] < 0).any() or (d[:, 1] <= d[:, 0]).any() or (d[:, 3] <= d[:, 2]).any():
            raise ValueError('Invalid destination cell')
        ri = np.maximum(d[:, None, 0], b[None, :, 0])
        ro = np.minimum(d[:, None, 1], b[None, :, 1])
        zl = np.maximum(d[:, None, 2], b[None, :, 2])
        zh = np.minimum(d[:, None, 3], b[None, :, 3])
        self.overlap_m3 = np.pi*np.maximum(ro**2-ri**2, 0)*np.maximum(zh-zl, 0)
        if not np.allclose(self.overlap_m3.sum(axis=1), self.volumes_m3, rtol=1e-10, atol=1e-18):
            raise ValueError('Fuel cells extend outside the saved active domain; no extrapolation allowed')
        self.fractions = self.overlap_m3/(history.volumes_cm3[None, :]*1e-6)
        self.initial_atoms = self.fractions@history.initial_atoms

    def fima_at(self, time_s):
        return (self.fractions@self.history.fissions_at(time_s))/self.initial_atoms

    def qdot_at(self, time_s):
        """Transfer extensive deposited watts, then divide by FE cell volume."""
        return (self.fractions@self.history.power_at(time_s))/self.volumes_m3

    def validate_complete_coverage(self, comm=None):
        covered = self.overlap_m3.sum(axis=0)
        if comm is not None:
            covered = comm.allreduce(covered)
        if not np.allclose(covered, self.history.volumes_cm3*1e-6, rtol=1e-10, atol=1e-18):
            raise ValueError('Destination fuel does not cover source exactly (missing or duplicate volume)')
        return float(np.max(np.abs(covered/(self.history.volumes_cm3*1e-6)-1)))


class NativeFIMACoupling:
    """Bind a validated history to persistent FEniCS cell coefficients.

    Initial atom density and cumulative fission density are retained for
    conservation diagnostics. Internal model.burnup and q_third are untouched.
    """

    def __init__(self, model, config):
        import dolfinx
        self.model = model
        self.material_name = config['material']
        material = model.materials[self.material_name]
        if model.regime != 'axisymmetric':
            raise ValueError('Native annular FIMA transfer requires axisymmetric r-z')
        if material.get('fima_source') != 'native_openmc':
            raise ValueError('Coupled fuel must explicitly select fima_source: native_openmc')
        self.history = DepletionHistory(config['history_path'])
        horizon = np.asarray(model.input_file.get('time', []), dtype=float)
        if len(horizon) and (horizon.min() < self.history.times[0] or horizon.max() > self.history.times[-1]):
            raise ValueError('Z3ST history extends beyond saved OpenMC times')
        selected = [name for name, card in model.materials.items() if card.get('fima_source') == 'native_openmc']
        if selected != [self.material_name]:
            raise ValueError('Native FIMA ownership must select exactly the configured material')
        tag = model.label_map[self.material_name]
        n_owned = model.mesh.topology.index_map(model.tdim).size_local
        self.cells = model.cell_tags.find(tag)
        self.cells = self.cells[self.cells < n_owned]
        bounds = []
        for cell in self.cells:
            xyz = model.mesh.geometry.x[model.mesh.geometry.dofmap[cell], :2]
            lo, hi = xyz.min(axis=0), xyz.max(axis=0)
            corners = np.array([[lo[0],lo[1]], [lo[0],hi[1]], [hi[0],lo[1]], [hi[0],hi[1]]])
            if len(xyz) != 4 or not all(np.any(np.all(np.isclose(xyz, point, rtol=0, atol=1e-12), axis=1)) for point in corners):
                raise ValueError('Overlap adapter currently supports axis-aligned first-order quad cells')
            bounds.append([lo[0],hi[0],lo[1],hi[1]])
        self.transfer = ConservativeTransfer(self.history, np.array(bounds).reshape(-1,4), config['axial_origin_cm'])
        self.coverage_error = self.transfer.validate_complete_coverage(model.mesh.comm)
        self.dofs = np.array([model.Q.dofmap.cell_dofs(int(c))[0] for c in self.cells], dtype=int)
        model.fima_native = dolfinx.fem.Function(model.Q, name='FIMA_native_OpenMC')
        model.swelling_eigenstrain_native = dolfinx.fem.Function(model.Q, name='Swelling_eigenstrain_native')
        model.initial_U_Zr_atom_density_native = dolfinx.fem.Function(model.Q, name='Initial_U_Zr_atom_density_native')
        model.cumulative_fission_density_native = dolfinx.fem.Function(model.Q, name='Cumulative_fission_density_native')
        model.initial_U_Zr_atom_density_native.x.array[self.dofs] = self.transfer.initial_atoms/self.transfer.volumes_m3
        model.initial_U_Zr_atom_density_native.x.scatter_forward()
        self.time_s = 0.0
        self.update(0.0)

    def update(self, time_s):
        values = self.transfer.fima_at(time_s)
        for fn in (self.model.fima_native, self.model.swelling_eigenstrain_native):
            fn.x.array[:] = 0
            fn.x.array[self.dofs] = values
            fn.x.scatter_forward()
        fn = self.model.cumulative_fission_density_native
        fn.x.array[:] = 0
        fn.x.array[self.dofs] = values*self.transfer.initial_atoms/self.transfer.volumes_m3
        fn.x.scatter_forward()
        self.time_s = float(time_s)


class NativePowerCoupling:
    """Opt-in absolute DG0 thermal source, independent of native FIMA.

    A positive, volume-conservative lumped CG1 projection feeds only the
    existing internal BU diagnostic. Thermal forms consume the DG0 source.
    """

    def __init__(self, model, config):
        import dolfinx
        import ufl
        self.model=model;self.material_name=config['material']
        if model.regime!='axisymmetric' or not model.on.get('thermal'):
            raise ValueError('Native heating requires axisymmetric thermal model')
        if model.on.get('fission_gas') or model.on.get('porosity'):
            raise ValueError('Native heating is not yet supported with fission-gas/porosity models')
        self.history=HeatingHistory(config['history_path'])
        horizon=np.asarray(model.input_file.get('time',[]),dtype=float)
        if len(horizon) and (horizon.min()<0 or horizon.max()>self.history.times[-1]):
            raise ValueError('Thermal history exceeds saved transport times')
        n_owned=model.mesh.topology.index_map(model.tdim).size_local
        cells=model.cell_tags.find(model.label_map[self.material_name])
        self.cells=cells[cells<n_owned]
        bounds=[]
        for cell in self.cells:
            xyz=model.mesh.geometry.x[model.mesh.geometry.dofmap[cell],:2]
            lo,hi=xyz.min(axis=0),xyz.max(axis=0)
            corners=np.array([[lo[0],lo[1]],[lo[0],hi[1]],[hi[0],lo[1]],[hi[0],hi[1]]])
            if len(xyz)!=4 or not all(np.any(np.all(np.isclose(xyz,p,rtol=0,atol=1e-12),axis=1)) for p in corners):
                raise ValueError('Native heating requires axis-aligned linear quad cells')
            bounds.append([lo[0],hi[0],lo[1],hi[1]])
        self.transfer=ConservativeTransfer(self.history,np.array(bounds).reshape(-1,4),config['axial_origin_cm'])
        self.transfer.validate_complete_coverage(model.mesh.comm)
        self.dofs=np.array([model.Q.dofmap.cell_dofs(int(c))[0] for c in self.cells],dtype=int)
        # Retain baseline coefficient independently; native mode bypasses it.
        model.q_third_baseline=model.q_third
        model.qdot_native=dolfinx.fem.Function(model.Q,name='qdot_native_OpenMC_W_m3')
        model.q_third=model.qdot_native
        self.burnup_source_nodal=dolfinx.fem.Function(model.V_t,name='qdot_native_BU_projection')
        v=ufl.TestFunction(model.V_t);x=ufl.SpatialCoordinate(model.mesh)
        self._load=dolfinx.fem.form(2*ufl.pi*x[0]*model.qdot_native*v*ufl.dx)
        self._mass=dolfinx.fem.assemble_vector(dolfinx.fem.form(2*ufl.pi*x[0]*v*ufl.dx))
        self._mass.scatter_reverse(dolfinx.la.InsertMode.add)
        self.time_s=0.0;self.update(0.0)

    def update(self,time_s):
        import dolfinx
        values=self.transfer.qdot_at(time_s)
        self.model.qdot_native.x.array[:]=0
        self.model.qdot_native.x.array[self.dofs]=values
        self.model.qdot_native.x.scatter_forward()
        load=dolfinx.fem.assemble_vector(self._load)
        load.scatter_reverse(dolfinx.la.InsertMode.add)
        self.burnup_source_nodal.x.array[:]=np.divide(load.array,self._mass.array,out=np.zeros_like(load.array),where=self._mass.array>0)
        self.burnup_source_nodal.x.scatter_forward()
        self.time_s=float(time_s)
