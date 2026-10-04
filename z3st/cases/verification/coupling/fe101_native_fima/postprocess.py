"""Phase-1 evidence, conservation, uniform comparison and spatial maps."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
from matplotlib.collections import PolyCollection

from z3st.coupling.openmc.depletion_fields import DepletionHistory
from prepare import CASE, SOURCE

ROOT = CASE.parents[4]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda:stream.read(8388608),b''):h.update(data)
    return h.hexdigest()


def history(mode):
    return [json.loads(line) for line in (CASE/'runs'/mode/'output/diagnostic_history.jsonl').read_text().splitlines()]


def main():
    state=json.loads((CASE/'verification_status.json').read_text())
    assert state['A_D_gate_before_native']=='PASS'
    src=DepletionHistory(SOURCE)
    native,uniform=history('native'),history('uniform_mean')
    assert len(native)==len(uniform)==21
    np.testing.assert_array_equal([r['time_s'] for r in native],src.times)
    mean_errors=[];fission_errors=[];comparison=[]
    for n,u in zip(native,uniform):
        source_mean=float(np.dot(src.fima_at(n['time_s']),src.volumes_cm3)/src.volumes_cm3.sum())
        np.testing.assert_allclose([n['FIMA_native_mean'],u['FIMA_native_mean']],source_mean,rtol=2e-13,atol=1e-18)
        mean_errors.append(abs(n['FIMA_native_mean']-source_mean))
        fission_errors.append(abs(n['mapped_cumulative_fissions']-src.fissions_at(n['time_s']).sum())/max(src.fissions_at(n['time_s']).sum(),1))
        comparison.append({'time_h':n['time_h'],'same_volume_mean_FIMA':source_mean,
                           'native_minus_uniform_Tmax_K':n['Tmax_K']-u['Tmax_K'],
                           'native_minus_uniform_gap_min_m':n['gap_min_m']-u['gap_min_m'],
                           'native_minus_uniform_fuel_max_ur_m':n['fuel_max_radial_displacement_m']-u['fuel_max_radial_displacement_m'],
                           'native_minus_uniform_fuel_von_mises_max_Pa':n['fuel_von_mises_max_Pa']-u['fuel_von_mises_max_Pa']})
    for mode in ['native','uniform_mean','zero_import','baseline_no_swelling','uniform_analytic','uniform_constant_reference']:
        rows=history(mode)
        scalars=[v for row in rows for k,v in row.items() if isinstance(v,(float,int))]
        assert np.isfinite(scalars).all()
        if mode in ['native','uniform_mean','zero_import']:
            assert rows[0]['FIMA_native_mean']==0
            assert np.all(np.diff([r['FIMA_native_mean'] for r in rows])>=0)
    # Verify the actual exported DG0 cell maps, not just diagnostic arrays.
    for mode in ['native','uniform_mean','zero_import']:
        files=sorted((CASE/'runs'/mode/'output').glob('fields_*.vtu'))
        assert len(files)==21
        field=np.load(CASE/'runs'/mode/'output/final_fields.npz')
        grid=pv.read(files[-1]);fuel_mask=grid.cell_data['MaterialID']==1
        np.testing.assert_allclose(grid.cell_data['FIMA_native_OpenMC'][fuel_mask],field['native_values'],rtol=0,atol=0)
        np.testing.assert_array_equal(grid.cell_data['FIMA_native_OpenMC'][~fuel_mask],0)
        np.testing.assert_array_equal(grid.cell_data['FIMA_native_OpenMC'],grid.cell_data['Swelling_eigenstrain_native'])
    # Neutron outputs and both original regression trees must stay byte-identical.
    before=json.loads((CASE/'protected_manifest_before.json').read_text())
    after={path:digest(ROOT/path) for path in before}
    assert before==after
    (CASE/'protected_manifest_after.json').write_text(json.dumps(after,indent=2)+'\n')
    audit=ROOT/'z3st/cases/verification/neutronics/fe101_5x5_depletion/preparation/high_statistics_40000_40_15/run_3500h_ED/postprocessing/reporting_correction_protected_files.json'
    archived=json.loads(audit.read_text())['sha256']
    unchanged_operators=[]
    for rel in ['z3st/models/mechanical_model.py','z3st/models/thermal_model.py','z3st/core/solver.py','z3st/core/finite_element_setup.py']:
        path=ROOT/rel
        assert digest(path)==archived[str(path)]
        unchanged_operators.append(rel)
    report={'phase':'1 FIMA -> isotropic swelling ONLY','verified':True,
            'FE_binding':json.loads((CASE/'binding_verification.json').read_text()),
            'not_fully_operational_coupled_simulation':True,'power_coupled':False,
            'source_sha256':src.sha256,'source_records':2100,
            'interface':json.loads((CASE/'interface_verification.json').read_text()),
            'C_uniform_regression':state['C'],'D_zero_regression':state['D'],
            'native_final':native[-1],'uniform_same_mean_final':uniform[-1],
            'native_vs_uniform_history':comparison,
            'all_21_times_max_abs_FE_mean_FIMA_error':max(mean_errors),
            'all_21_times_max_relative_fission_sum_error':max(fission_errors),
            'protected_files_unchanged':len(before),'FE_operators_unchanged':unchanged_operators,
            'resources_and_runtimes':{k:{'resources':v['resources_before_solve'],'runtime_s':v['runtime_s']} for k,v in state['solves'].items()},
            'approximations':['DG0 cylindrical overlap averaging on axis-aligned first-order quads; no general unstructured remapper.',
                             'Linear interpolation of cumulative fissions between saved depletion times, no extrapolation.',
                             'Fixed baseline thermal LHR 8780 W/m, insulated ends, outer clad 300 K; no OpenMC power transfer.',
                             'Gap/contact remain mean-interface models; local geometric gaps are diagnostic only.',
                             'U-ZrH swelling correlation unchanged and not independently calibrated by this verification.',
                             'Point stress exports at 31 r=0 axis nodes have u_r/r singularities and are zeroed by the existing writer; reported stress extrema use finite DG0 cell samples.'],
            'technical_retry':'First uniform solve converged but diagnostics used cells instead of cells0 in FEniCSx 0.11; fixed only diagnostics, retained failed attempt, reran with preflight. The initial postprocessing invocation preceded script creation; resumed postprocessing from saved outputs without rerunning solves.'}
    (CASE/'phase1_results.json').write_text(json.dumps(report,indent=2)+'\n')
    fields=np.load(CASE/'runs/native/output/final_fields.npz')
    controls=np.load(CASE/'runs/uniform_mean/output/final_fields.npz')
    bounds=fields['native_bounds_m']
    polys=[[(b[0]*100,b[2]*100),(b[1]*100,b[2]*100),(b[1]*100,b[3]*100),(b[0]*100,b[3]*100)] for b in bounds]
    fig,axes=plt.subplots(1,3,figsize=(13,6),layout='constrained')
    vals=fields['native_values']
    for ax,array,title,label in zip(axes,[vals,vals,vals-controls['native_values']],
                                   ['Imported native FIMA','Isotropic swelling eigenstrain','Native minus equal-mean uniform'],
                                   ['FIMA [fraction]','epsilon_sw,ii [strain]','Delta FIMA [fraction]']):
        artist=PolyCollection(polys,array=array,cmap='viridis' if ax!=axes[2] else 'RdBu_r',edgecolors='none')
        ax.add_collection(artist);ax.autoscale_view();ax.set(xlabel='r [cm]',ylabel='local z [cm]',title=title)
        fig.colorbar(artist,ax=ax,label=label)
    fig.suptitle('PHASE 1 verification, 3500 h — baseline thermal power; no q coupling')
    for extension in ['png','pdf']:fig.savefig(CASE/('phase1_FIMA_swelling_maps.'+extension),dpi=180)
    plt.close(fig)
    n,u=native[-1],uniform[-1]
    lines=['# FE101 phase 1: native FIMA → swelling verification','',
           '**Verification case only; thermal power is the Z3ST baseline, not OpenMC operational power.**','',
           'Coordinate transform: `r_FE_m=r_OpenMC_cm/100`, `z_FE_m=(z_OpenMC_cm−10.20)/100`. The new isolated mesh is 0.3556 m long; the existing 0.356 m regression is unchanged. No stretching, clamping or extrapolation. History: all 21 times through 3500 h, 2100 source records; noncoincident times use linear interpolation of cumulative fissions.','',
           'Imported native FIMA and internally accumulated baseline BU are independent. Explicit `fima_source: native_openmc` prevents a silent legacy fallback. The law remains `epsilon_sw=FIMA*I`, `DeltaV/V=3*FIMA`. No q\'\'\' import, creep or hydrogen redistribution.','',
           'A/B interface tests: PASS (seven tests). Volume relative error {:.3e}; initial atom error {:.3e}; final volume-weighted mean FIMA error {:.3e}; worst fission-sum relative error {:.3e}.'.format(report['interface']['volume_relative_error'],report['interface']['initial_atoms_relative_error'],report['interface']['final_FIMA_mean_relative_error'],report['interface']['max_relative_fission_sum_error']),'',
           'Native radial/axial trends retained. Peak change {:.7f}%; maximum re-binned radial-profile difference {:.7f}%; maximum individual re-binned domain difference {:.7f}%. Integrals are conserved but mixed radial cells average bin jumps.'.format(report['interface']['peak_change_percent'],report['interface']['max_radial_profile_remap_difference_percent'],report['interface']['max_bin_remap_difference_percent']),'',
           'C uniform: PASS. Analytic free-swelling displacement error {:.3e} m; difference from the existing constant volumetric-swelling channel {:.3e} m; tolerance 1e-9 m.'.format(state['C']['analytic_displacement_max_abs_error_m'],state['C']['constant_eigenstrain_control_max_abs_difference_m']),'',
           'D zero: PASS. Temperature difference from no-swelling baseline {:.3e} K; displacement difference {:.3e} m. Internal BU remains nonzero, proving imported zero does not fall back to BU-driven swelling.'.format(state['D']['Tmax_abs_difference_K'],state['D']['umax_abs_difference_m']),'',
           'A–D were PASS before the native solve. Native and uniform-equal-mean solves both completed all 21 times. The uniform control has the same mean at **every saved time**, not merely at the final time.','',
           '| Final metric | Native FIMA | Uniform same-mean FIMA |','|---|---:|---:|']
    for key,scale,label in [('Tmax_K',1,'Fuel Tmax [K]'),('fuel_surface_T_min_K',1,'Fuel surface Tmin [K]'),('fuel_surface_T_mean_K',1,'Fuel surface Tmean [K]'),('fuel_surface_T_max_K',1,'Fuel surface Tmax [K]'),('gap_min_m',1e6,'Local gap minimum [um]'),('gap_mean_m',1e6,'Local gap mean [um]'),('gap_max_m',1e6,'Local gap maximum [um]'),('model_contact_pressure_Pa',1,'Model contact pressure [Pa]'),('fuel_max_radial_displacement_m',1e6,'Fuel maximum radial displacement [um]'),('clad_inner_ur_min_m',1e6,'Clad inner radial displacement min [um]'),('clad_inner_ur_mean_m',1e6,'Clad inner radial displacement mean [um]'),('clad_inner_ur_max_m',1e6,'Clad inner radial displacement max [um]'),('FIMA_native_mean',1,'Imported FIMA volume mean'),('swelling_eigenstrain_max',1,'Maximum isotropic swelling eigenstrain'),('internal_BU_mean_MWd_kgU',1,'Internal baseline BU [MWd/kgU]'),('fuel_von_mises_max_Pa',1e-6,'Fuel cell von Mises max [MPa]'),('clad_von_mises_max_Pa',1e-6,'Clad cell von Mises max [MPa]')]:
        lines.append(f'| {label} | {n[key]*scale:.12g} | {u[key]*scale:.12g} |')
    lines += ['',f"Native maximum FIMA: {n['FIMA_native_max']:.15g}, cell bounds [ri,ro,zlo,zhi] in local metres: {n['FIMA_max_cell_bounds_m']}; corresponding native source domain r10/z6 (27.98–31.536 cm absolute z).",'',
              'Stress component extrema are DG0 cell-centre samples of total elastic stress (thermal and swelling eigenstress included), not continuum extrema or nodal axis values.','',
              '| Region/component | Native min / max [MPa] | Uniform min / max [MPa] |','|---|---:|---:|']
    for material in ['fuel','clad']:
        for component in ['rr','hoop','zz','rz']:
            low=material+'_sigma_'+component+'_min_Pa';high=material+'_sigma_'+component+'_max_Pa'
            lines.append(f'| {material} {component} | {n[low]/1e6:.9g} / {n[high]/1e6:.9g} | {u[low]/1e6:.9g} / {u[high]/1e6:.9g} |')
    lines += ['',f"Contact: native={n['model_contact_pressure_Pa']>0}, uniform={u['model_contact_pressure_Pa']>0}; local geometric closure native={n['local_gap_closed']}, uniform={u['local_gap_closed']}.",'',
              'Resources before each successful solve (OMP=2, MPI=1, BLAS/NumExpr=1; UCX_TLS=self):','',
              '| Solve | RAM available [GiB] | Swap used [GiB] | Disk available [GiB] | CPUs | Runtime [s] |','|---|---:|---:|---:|---:|---:|']
    for mode,v in report['resources_and_runtimes'].items():
        r=v['resources'];lines.append(f"| {mode} | {r['RAM_available_GiB']:.4f} | {r['swap_used_bytes']/1024**3:.4f} | {r['disk_available_GiB']:.4f} | {r['logical_CPUs']} | {v['runtime_s']:.3f} |")
    lines += ['',f"Protected manifest: PASS, {len(before)} unchanged files, including all neutron outputs and both original regression trees. Four FE operator files verified identical to the earlier audit hashes.",'',
              'FE binding tests: PASS. Imported FIMA is independent of BU and q; diagonal swelling equals FIMA, off-diagonal swelling is zero, clad FIMA is zero, adaptive rollback restores coefficient/time state, out-of-range requests fail without state changes, and a missing native field raises instead of falling back to BU. Legacy BU mode remains available.','',
              'Approximations/limits:','']+['- '+text for text in report['approximations']]
    lines += ['',report['technical_retry'],'',
              '**Phase 1 verified within the documented data-transfer, baseline thermal and averaged-gap assumptions.** Native heterogeneity produces the recorded mechanical differences from the equal-mean control; these are deterministic FE responses to a Monte Carlo input field, not a new statistical-significance certification. Phase 2 power coupling is not implemented or started.','',
              'Spatial maps: phase1_FIMA_swelling_maps.png/pdf. Full cell maps and 21-time fields: runs/native/output/fields_*.vtu. Full diagnostics and comparison: phase1_results.json.','',
              'No OpenMC/depletion execution, neutron/regression edits, commit or push.']
    (CASE/'PHASE1_REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'phase1_verified':True,'native_Tmax_K':n['Tmax_K'],'native_gap_min_um':n['gap_min_m']*1e6,'native_FIMA_max':n['FIMA_native_max'],'protected_files':len(before)},indent=2))


if __name__=='__main__':main()
