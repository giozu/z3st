"""Power integrals, comparison controls, native two-field maps and protection."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from prepare import CASE,ROOT,PHASE1,SOURCE
from z3st.coupling.openmc.depletion_fields import HeatingHistory

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(8388608),b''):h.update(x)
    return h.hexdigest()

def records(mode):return [json.loads(x) for x in (CASE/'runs'/mode/'output/diagnostic_history.jsonl').read_text().splitlines()]

def main():
    state=json.loads((CASE/'verification_status.json').read_text());assert state['A_D_before_fully_coupled']=='PASS'
    h=HeatingHistory(SOURCE);r=records('full_native');assert len(r)==21
    power=[]
    for row in r:
        assert row['power_time_s']==row['FIMA_time_s']==row['time_s']
        power.append({k:row[k] for k in ['time_h','fuel_power_source_W','fuel_power_FE_W','power_error_W','power_relative_error']})
    with (CASE/'power_conservation_FE.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=power[0].keys());w.writeheader();w.writerows(power)
    before=json.loads((CASE/'protected_manifest_before.json').read_text());after={p:digest(ROOT/p) for p in before};assert before==after
    (CASE/'protected_manifest_after.json').write_text(json.dumps(after,indent=2)+'\n')
    cases={'phase1_saved':json.loads((PHASE1/'phase1_results.json').read_text())['native_final']}
    for mode in ['zero_native','uniform_native','uniform_baseline','native_thermal','native_uniform_fima','uniform_power_native_fima','full_native','phase1_replay']:
        cases[mode]=records(mode)[-1]
    # Verify both new coefficient exports, including fuel-only support.
    g=pv.read(CASE/'runs/full_native/output/fields_0020.vtu');mask=g.cell_data['MaterialID']==1
    q=g.cell_data['qdot_native_OpenMC_W_m3'];fi=g.cell_data['FIMA_native_OpenMC']
    assert np.isfinite(q).all() and (q>=0).all() and (q[~mask]==0).all()
    np.testing.assert_array_equal(fi,g.cell_data['Swelling_eigenstrain_native'])
    fields=np.load(CASE/'runs/full_native/output/final_fields.npz');therm=np.load(CASE/'runs/full_native/output/thermal_fields.npz')
    np.testing.assert_array_equal(q[mask],therm['qdot_values'])
    np.testing.assert_array_equal(fi[mask],fields['native_values'])
    # Conservative nodal BU source is internally accumulated, not imported BU.
    result={'phase2_verified':True,'one_way_only':True,'source':json.loads((CASE/'source_metadata.json').read_text()),
            'binding_tests':json.loads((CASE/'binding_verification.json').read_text()),
            'source_HDF5_check':json.loads((CASE/'source_heating_verification.json').read_text()),
            'tests':{k:state[k] for k in ['A','B','C','D','A_D_before_fully_coupled','phase1_compatibility']},
            'cases_final':cases,'power_conservation_21_times':power,
            'max_abs_power_error_W':max(abs(x['power_error_W']) for x in power),
            'max_abs_relative_power_error':max(abs(x['power_relative_error']) for x in power),
            'protected_files_unchanged':len(before),'resources_and_runtimes':state['solves'],
            'limitations':['One-way transfer from fixed OpenMC states; no feedback to neutronics.',
                           'Global mean gap/contact conductance and pressure unchanged; local gaps are diagnostics.',
                           'Axis-aligned rectangular DG0 overlap averaging; bin jumps can be smoothed.',
                           'Instantaneous power interpolated linearly; internal BU uses existing right-endpoint accumulation of a conservative positive CG1 projection, not saved BOS BU.',
                           'Unchanged material correlations and outer 300 K/insulated-end boundary conditions; no coolant-energy model.',
                           'Stresses sampled at DG0 cell centres; existing singular axis nodal exports are not used for extrema.']}
    (CASE/'phase2_results.json').write_text(json.dumps(result,indent=2)+'\n')
    polygons=g.cells.reshape(-1,5)[:,1:][mask];xyz=g.points[:,:2]*100
    polys=[xyz[cell] for cell in polygons]
    fig,axes=plt.subplots(2,2,figsize=(12,10),layout='constrained')
    for ax,data,title in [(axes[0,0],q[mask],"q''' native [W/m3]"),(axes[0,1],g.point_data['Temperature'][polygons].mean(axis=1),'Temperature [K], cell-average visualization'),(axes[1,0],fi[mask],'FIMA native [fraction]'),(axes[1,1],fi[mask],'Isotropic swelling eigenstrain [strain]')]:
        a=PolyCollection(polys,array=data,cmap='viridis',edgecolors='none');ax.add_collection(a);ax.autoscale_view();ax.set(xlabel='r [cm]',ylabel='local z [cm]',title=title);fig.colorbar(a,ax=ax)
    fig.suptitle('FE101 — first one-way native FIMA + heating verification, 3500 h')
    for ext in ['png','pdf']:fig.savefig(CASE/('phase2_native_maps.'+ext),dpi=170)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for mode in ['phase1_replay','native_uniform_fima','uniform_power_native_fima','full_native']:
        v=np.load(CASE/'runs'/mode/'output/thermal_fields.npz');z=v['surface_z_m']*100
        axes[0].plot(z,v['gap_m']*1e6,label=mode);axes[1].plot(z,v['fuel_ur_m']*1e6,label=mode)
    for ax,title in zip(axes,['Local geometric gap [um]','Fuel surface radial displacement [um]']):ax.set(xlabel='local z [cm]',ylabel=title);ax.legend(fontsize=8);ax.grid(alpha=.3)
    for ext in ['png','pdf']:fig.savefig(CASE/('phase2_gap_deformation.'+ext),dpi=170)
    plt.close(fig)
    modes=['phase1_saved','uniform_power_native_fima','native_uniform_fima','full_native'];n=cases['full_native']
    lines=['# FE101 FASE 2 — FIMA nativa + heating nativo','',
           '**Coupling one-way dei due campi verificato; nessun nuovo calcolo OpenMC/depletion.**','',
           'Sorgente: `zest_FIMA_2D_10x10.csv`, colonna `deposited_power_W`. Verificati in sola lettura tutti i 21 HDF5: tally `FE101 B1 domain energy deposition`, score `heating-local` in eV/source; normalizzazione con `factor-for-normalization` (`heating-local` globale) a 250000 W di reattore. P_i=250000 H_i/H_core [W]; q_i=P_i/(volume_cm3*1e-6) [W/m3]. Non si usano FIMA o BU per ricostruire potenza.','',
           'Coordinate: r_FE=r_OpenMC/100; z_FE=(z_OpenMC−10.20)/100, lunghezza 0.3556 m. Overlap cilindrico: potenza estensiva trasferita e divisa per volume FE. Nessuno stretching, clamp o rinormalizzazione alla baseline.','',
           'Tempi: 21 snapshot 0–3500 h. Interpolazione lineare della potenza istantanea, non negativa; FIMA cumulativa interpolata alla stessa coordinata temporale. Nessuna extrapolazione. Questa scelta termica non è l’integrazione BOS dell’energia della depletion. BU interno resta distinto e non è usato dallo swelling nativo.','',
           f"Conservazione FE: massimo errore assoluto {result['max_abs_power_error_W']:.3e} W; massimo relativo {result['max_abs_relative_power_error']:.3e}; tolleranza gate 1e-11 relativa (1e-10 W per campo nullo). Tutti i tempi e gli errori sono in power_conservation_FE.csv.",'',
           f"A PASS: dati/overlap/unità/tempi/off-grid/positività. B PASS: sorgente importata nulla, max|T−300 K|={state['B']['max_abs_T_minus_BC_K']:.3e} K anche con LHR baseline non nullo. C PASS: 5000 W uniformi via DG0 contro sorgente Z3ST originale, max ΔT={state['C']['max_abs_T_difference_K']:.3e} K, tolleranza 1e-5 K. D PASS: thermal-only nativo completo, senza swelling. A–D PASS prima del caso completamente accoppiato.",'',
           f"Regressione FASE 1 in copia isolata: PASS; ΔT={state['phase1_compatibility']['T_error_K']:.3e} K, Δu={state['phase1_compatibility']['u_error_m']:.3e} m. Originali FASE 1 invariati.",'',
           '| Metrica finale | Baseline + FI nativa (FASE1) | q uniforme stessa P + FI nativa | q nativa + FI uniforme | q nativa + FI nativa |','|---|---:|---:|---:|---:|']
    for key,scale,label in [('Tmax_K',1,'Tmax fuel [K]'),('fuel_surface_T_min_K',1,'T superficie min [K]'),('fuel_surface_T_mean_K',1,'T superficie media [K]'),('fuel_surface_T_max_K',1,'T superficie max [K]'),('gap_min_m',1e6,'gap min [um]'),('gap_mean_m',1e6,'gap medio [um]'),('gap_max_m',1e6,'gap max [um]'),('fuel_max_radial_displacement_m',1e6,'ur fuel max [um]'),('clad_inner_ur_mean_m',1e6,'ur clad interno medio [um]'),('fuel_von_mises_max_Pa',1e-6,'von Mises fuel max [MPa]'),('clad_von_mises_max_Pa',1e-6,'von Mises clad max [MPa]'),('FIMA_native_mean',1,'FIMA media'),('swelling_eigenstrain_max',1,'eigenstrain swelling max')]:
        lines.append('| '+label+' | '+' | '.join(f'{cases[m][key]*scale:.12g}' for m in modes)+' |')
    lines+=['',f"Caso completo: q min/media/max={n['qdot_min_W_m3']:.12g}/{n['qdot_mean_W_m3']:.12g}/{n['qdot_max_W_m3']:.12g} W/m3; P_fuel={n['fuel_power_FE_W']:.12g} W; FIMA min/media/max={n['FIMA_native_min']:.12g}/{n['FIMA_native_mean']:.12g}/{n['FIMA_native_max']:.12g}.",'',
            f"Tmax in [r,z] locali m: {n['Tmax_position_local_m']}; q max nel box [ri,ro,zlo,zhi] m: {n['qdot_max_cell_bounds_m']}; FIMA max nel box {n['FIMA_max_cell_bounds_m']}; gap minimo a z={n['gap_min_position_local_z_m']} m.",'',
            f"Contatto medio: p={n['model_contact_pressure_Pa']} Pa; chiusura geometrica locale={n['local_gap_closed']}. ur clad interno min/medio/max={n['clad_inner_ur_min_m']*1e6:.9g}/{n['clad_inner_ur_mean_m']*1e6:.9g}/{n['clad_inner_ur_max_m']*1e6:.9g} um.",'',
            'Separazione degli effetti:','',
            '- Potenza totale: uniforme alla potenza OpenMC + FI nativa contro FASE 1, mantenendo FI e forma termica uniforme.','- Eterogeneità termica: q nativa + FI nativa contro q uniforme stessa potenza + FI nativa.','- Eterogeneità swelling: q nativa + FI nativa contro q nativa + FI uniforme alla stessa media a ogni tempo.','- Questi sono effetti deterministici della scelta dei campi nel modello FE; non dimostrano feedback neutronico o significatività statistica delle differenze.','',
            'Thermal-only nativo (D):','',json.dumps(cases['native_thermal'],indent=2),'',
            'Risorse prima di ogni solve (output free/df/nproc/vmstat nei resource_preflight.json):','',
            '| Caso | RAM available GiB | Swap used GiB | Disco GiB | CPU | OMP | Runtime s |','|---|---:|---:|---:|---:|---:|---:|']
    for mode,solve in state['solves'].items():
        rr=solve['resources'];lines.append(f"| {mode} | {rr['RAM_available_GiB']:.4f} | {rr['swap_used_bytes']/1024**3:.4f} | {rr['disk_free_GiB']:.4f} | {rr['logical_CPUs']} | 2 | {solve['runtime_s']:.3f} |")
    lines+=['',f"Protezione PASS: {len(before)} file invariati: neutronica, intera FASE 1, regression 101/PWR e sorgenti non coinvolti.",'',
            'Weak form termo-meccaniche e misure invariate. La sola modifica a thermal_model.py seleziona lo spazio del coefficiente per la stampa diagnostica q_third; il weak form è identico. Swelling e gap/contact non modificati.','',
            'Limiti:','']+['- '+x for x in result['limitations']]
    lines+=['','Mappe: phase2_native_maps.png/pdf; gap/deformazione: phase2_gap_deformation.png/pdf; tutti i campi ai 21 tempi nelle directory runs/*/output.','',
            '**FASE 2 e primo trasferimento FIMA + q verificati entro questi limiti. Nessuna nuova fisica o feedback bidirezionale. Stop per revisione; nessun commit/push.**']
    (CASE/'PHASE2_REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'verified':True,'Tmax_K':n['Tmax_K'],'fuel_power_W':n['fuel_power_FE_W'],'gap_min_um':n['gap_min_m']*1e6,'power_max_relative_error':result['max_abs_relative_power_error'],'protected_files':len(before)},indent=2))

if __name__=='__main__':main()
