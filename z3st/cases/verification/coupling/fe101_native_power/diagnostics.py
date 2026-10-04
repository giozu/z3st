"""Single-process verification diagnostics, including local geometric gap.

The local gap is reconstructed from surface displacements; model contact
pressure/conductance still use their existing whole-interface mean.
"""
import json
from pathlib import Path
import numpy as np
import dolfinx
import ufl


def per_step(problem, step, t):
    if problem.mesh.comm.size != 1:
        raise RuntimeError('Verification diagnostics require one MPI process')
    p = problem
    xyz = p.V_t.tabulate_dof_coordinates()
    fuel = p.mgr.locate_domain_dofs(p.label_map['fuel'],p.V_t)
    sf = p.mgr.locate_facets_dofs(p.label_map['lateral_1'],p.V_t)
    sc = p.mgr.locate_facets_dofs(p.label_map['inner_2'],p.V_t)
    uf = p.mgr.locate_domain_dofs(p.label_map['fuel'],p.V_m)
    uv = p.u.x.array.reshape(-1,2)
    ux = p.V_m.tabulate_dof_coordinates()
    # Scalar/vector spaces have the same degree but do not assume dof ordering.
    def surface_displacement(dofs):
        positions = xyz[dofs]
        values = np.array([uv[np.argmin(np.linalg.norm(ux-point,axis=1)),0] for point in positions])
        order = np.argsort(positions[:,1])
        return positions[order,1],values[order]
    zf,ur_f = surface_displacement(sf)
    zc,ur_c = surface_displacement(sc)
    gap = p.geometry['inner_radius_2']-p.geometry['outer_radius_1']+np.interp(zf,zc,ur_c)-ur_f
    row = {'step':int(step),'time_s':float(t),'time_h':float(t/3600),
           'Tmax_K':float(p.T.x.array[fuel].max()),
           'fuel_surface_T_min_K':float(p.T.x.array[sf].min()),
           'fuel_surface_T_mean_K':float(np.trapezoid(p.T.x.array[sf][np.argsort(xyz[sf,1])],zf)/(zf[-1]-zf[0])),
           'fuel_surface_T_max_K':float(p.T.x.array[sf].max()),
           'gap_min_m':float(gap.min()),'gap_mean_m':float(np.trapezoid(gap,zf)/(zf[-1]-zf[0])),
           'gap_max_m':float(gap.max()),'model_mean_gap_m':float(p._last_gap),
           'model_contact_pressure_Pa':float(p._last_pressure),
           'local_gap_closed':bool((gap <= 0).any()),
           'fuel_max_radial_displacement_m':float(uv[uf,0].max()),
           'clad_inner_ur_min_m':float(ur_c.min()),'clad_inner_ur_mean_m':float(np.trapezoid(ur_c,zc)/(zc[-1]-zc[0])),
           'clad_inner_ur_max_m':float(ur_c.max()),
           'internal_BU_mean_MWd_kgU':float(p.material_mean(p.burnup,'fuel'))}
    native = getattr(p,'_fima_coupling',None)
    if native is not None:
        vals = p.fima_native.x.array[native.dofs]
        bounds = native.transfer.destination_bounds_m
        index = int(np.argmax(vals))
        row.update({'FIMA_native_mean':float(p.material_mean(p.fima_native,'fuel')),
                    'FIMA_native_max':float(vals.max()),'swelling_eigenstrain_max':float(vals.max()),
                    'FIMA_max_cell_bounds_m':bounds[index].tolist(),
                    'FIMA_max_cells_bounds_m':bounds[np.isclose(vals,vals.max(),rtol=1e-13,atol=1e-18)].tolist(),
                    'mapped_initial_atoms':float(native.transfer.initial_atoms.sum()),
                    'mapped_cumulative_fissions':float(np.dot(vals,native.transfer.initial_atoms)),
                    'coverage_max_relative_error':native.coverage_error})
    row['Tmax_position_local_m']=xyz[fuel[np.argmax(p.T.x.array[fuel])],:2].tolist()
    row['gap_min_position_local_z_m']=float(zf[np.argmin(gap)])
    if native is not None:
        row['FIMA_native_min']=float(vals.min())
    power=getattr(p,'_power_coupling',None)
    if power is not None:
        q=p.qdot_native.x.array[power.dofs]
        row.update({'qdot_min_W_m3':float(q.min()),'qdot_mean_W_m3':float(np.average(q,weights=power.transfer.volumes_m3)),
                    'qdot_max_W_m3':float(q.max()),'qdot_max_cell_bounds_m':power.transfer.destination_bounds_m[int(np.argmax(q))].tolist(),
                    'fuel_power_FE_W':float(p.material_mean(p.qdot_native,'fuel',integral=True)),
                    'fuel_power_source_W':float(power.history.power_at(t).sum()),
                    'power_time_s':power.time_s,'FIMA_time_s':native.time_s if native is not None else None})
        row['power_error_W']=row['fuel_power_FE_W']-row['fuel_power_source_W']
        row['power_relative_error']=row['power_error_W']/row['fuel_power_source_W'] if row['fuel_power_source_W'] else 0.0
        assert abs(row['power_error_W']) <= max(1e-10,1e-11*row['fuel_power_source_W'])
        if native is not None:assert native.time_s==power.time_s==t
    # Evaluate total physical stress, with thermal + swelling eigenstress.
    for name in ['fuel','clad']:
        cells = p.cell_tags.find(p.label_map[name])
        space = dolfinx.fem.functionspace(p.mesh,('DG',0,(3,3)))
        fn = dolfinx.fem.Function(space)
        expr = dolfinx.fem.Expression(p.stress[name],space.element.interpolation_points)
        fn.interpolate(expr,cells0=cells)
        indices = np.array([space.dofmap.cell_dofs(int(c))[0] for c in cells])
        tensors = fn.x.array.reshape(-1,3,3)[indices]
        dev = tensors-np.trace(tensors,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
        row[name+'_von_mises_max_Pa'] = float(np.sqrt(1.5*np.sum(dev*dev,axis=(1,2))).max())
        for component,i,j in [('rr',0,0),('hoop',1,1),('zz',2,2),('rz',0,2)]:
            row[name+'_sigma_'+component+'_min_Pa'] = float(tensors[:,i,j].min())
            row[name+'_sigma_'+component+'_max_Pa'] = float(tensors[:,i,j].max())
    out = Path('output');out.mkdir(exist_ok=True)
    with (out/'diagnostic_history.jsonl').open('w' if step==0 else 'a') as stream:
        stream.write(json.dumps(row,allow_nan=False)+'\n')
    np.savez(out/'thermal_fields.npz',qdot_values=p.qdot_native.x.array[power.dofs] if power is not None else np.array([]),bounds_m=power.transfer.destination_bounds_m if power is not None else np.empty((0,4)),surface_z_m=zf,surface_T_K=p.T.x.array[sf][np.argsort(xyz[sf,1])],gap_m=gap,fuel_ur_m=ur_f,clad_ur_m=np.interp(zf,zc,ur_c))
    np.savez(out/'final_fields.npz',scalar_coords=xyz,vector_coords=ux,T=p.T.x.array,u=uv,
             FIMA_native=p.fima_native.x.array if native is not None else np.array([]),
             native_bounds_m=native.transfer.destination_bounds_m if native is not None else np.empty((0,4)),
             native_values=p.fima_native.x.array[native.dofs] if native is not None else np.array([]),
             native_volumes_m3=native.transfer.volumes_m3 if native is not None else np.array([]))
