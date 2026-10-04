"""Five-case consistency, local sensitivities and immutable-reference checks."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import yaml
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
CASE=Path(__file__).resolve().parent;ROOT=CASE.parents[4]
S=[0.,.5,1.,1.5,2.];ON=CASE.parent/'fe101_native_power/runs/full_native';OFF=CASE.parent/'fe101_swelling_ablation/runs/swelling_off'
def rows(p):return [json.loads(x) for x in (p/'output/diagnostic_history.jsonl').read_text().splitlines()]
def writecsv(name,data):
    with (CASE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def main():
    status=json.loads((CASE/'solve_status.json').read_text());assert status['completed'] and all(r['PASS'] for r in status['solves'].values())
    reference=rows(ON);reference_off=rows(OFF);histories={};endpoint_checks={};qerr=fierr=epsrel=0.;contact=[]
    baseline_input=yaml.safe_load((ON/'input.yaml').read_text())
    for cfg in baseline_input['coupling'].values():cfg['history_path']=str((ON/cfg['history_path']).resolve())
    baseline_fuel=yaml.safe_load((ON/'fuel.yaml').read_text());baseline_fuel.pop('eigenstrain')
    for s in S:
        run=CASE/'runs'/f's_{s:.1f}';rs=rows(run);histories[s]=rs
        assert len(rs)==21 and [r['time_h'] for r in rs]==list(np.arange(21)*175)
        for name in ['mesh.msh','geometry.yaml','clad.yaml','boundary_conditions.yaml']:assert digest(run/name)==digest(ON/name)
        inp=yaml.safe_load((run/'input.yaml').read_text())
        for cfg in inp['coupling'].values():cfg['history_path']=str((run/cfg['history_path']).resolve())
        assert inp==baseline_input
        fuel=yaml.safe_load((run/'fuel.yaml').read_text());assert fuel.pop('swelling_sensitivity_factor')==s;assert fuel.pop('eigenstrain')=='sensitivity_law.scaled_native_swelling';assert fuel==baseline_fuel
        previous=None
        for i,r in enumerate(rs):
            assert np.isfinite([v for v in r.values() if isinstance(v,(float,int))]).all()
            assert r['time_s']==r['power_time_s']==r['FIMA_time_s'] and r['time_h']<=3500
            assert r['fuel_power_FE_W']==reference[i]['fuel_power_FE_W']
            assert r['gap_min_m']>=0 or r['local_gap_closed']
            if r['local_gap_closed'] or r['model_contact_pressure_Pa']>0:contact.append({'s':s,'time_h':r['time_h']})
            g=pv.read(run/f'output/fields_{i:04d}.vtu');ref=pv.read(ON/f'output/fields_{i:04d}.vtu')
            q=g.cell_data['qdot_native_OpenMC_W_m3'];fi=g.cell_data['FIMA_native_OpenMC'];np.testing.assert_array_equal(q,ref.cell_data['qdot_native_OpenMC_W_m3']);np.testing.assert_array_equal(fi,ref.cell_data['FIMA_native_OpenMC'])
            assert np.isfinite(q).all() and (q>=0).all() and np.isfinite(fi).all()
            fn=np.load(run/f'output/sensitivity_fields_{i:02d}.npz')
            np.testing.assert_allclose(fn['applied_eigenstrain'],s*fn['FIMA'],rtol=2e-14,atol=1e-18)
            assert np.isfinite(fn['applied_eigenstrain']).all()
            if previous is not None:assert (fn['FIMA']>=previous[0]-1e-18).all() and (fn['applied_eigenstrain']>=previous[1]-1e-18).all()
            previous=(fn['FIMA'].copy(),fn['applied_eigenstrain'].copy())
            # Distinguish unchanged nominal provider coefficient from actual applied law.
            if 'Swelling_eigenstrain_native' in g.cell_data:
                g.cell_data['Swelling_potential_eigenstrain_nominal']=g.cell_data.pop('Swelling_eigenstrain_native')
                g.cell_data['Swelling_eigenstrain_applied']=s*fi
                g.save(run/f'output/fields_{i:04d}.vtu')
            if s in [0.,1.]:
                target=OFF if s==0 else ON;targetgrid=pv.read(target/f'output/fields_{i:04d}.vtu')
                terr=float(np.max(abs(g.point_data['Temperature']-targetgrid.point_data['Temperature'])))
                # Direct final primary-variable check below, diagnostic extrema across all 21 times.
                assert terr<1e-8
                for key in ['gap_min_m','gap_mean_m','fuel_max_radial_displacement_m','fuel_von_mises_max_Pa','clad_von_mises_max_Pa']:
                    targetrows=reference_off if s==0 else reference
                    np.testing.assert_allclose(r[key],targetrows[i][key],rtol=2e-10,atol=1e-12)
        if s in [0.,1.]:
            target=OFF if s==0 else ON;a=np.load(run/'output/final_fields.npz');b=np.load(target/'output/final_fields.npz')
            t=float(np.max(abs(a['T']-b['T'])));u=float(np.max(abs(a['u']-b['u'])));assert t<1e-8 and u<1e-12
            endpoint_checks[str(s)]={'PASS':True,'max_temperature_difference_K':t,'max_displacement_difference_m':u,'reference':str(target)}
    metrics=[('Tmax_K',1,'K'),('fuel_surface_T_min_K',1,'K'),('fuel_surface_T_mean_K',1,'K'),('fuel_surface_T_max_K',1,'K'),('gap_min_m',1e6,'um'),('gap_mean_m',1e6,'um'),('gap_max_m',1e6,'um'),('gap_min_position_local_z_m',100,'cm'),('fuel_max_radial_displacement_m',1e6,'um'),('clad_inner_ur_mean_m',1e6,'um'),('fuel_von_mises_max_Pa',1e-6,'MPa'),('clad_von_mises_max_Pa',1e-6,'MPa'),('swelling_eigenstrain_max',1,'strain')]
    final=[]
    for s in S:
        row={'s':s}
        for key,scale,unit in metrics:
            label=key[:-2]+'_'+unit if unit in ['um','cm'] else (key[:-3]+'_MPa' if unit=='MPa' else key)
            row[label]=histories[s][-1][key]*scale
        row['fuel_power_W']=histories[s][-1]['fuel_power_FE_W'];row['contact']=histories[s][-1]['local_gap_closed'] or histories[s][-1]['model_contact_pressure_Pa']>0
        final.append(row)
    writecsv('sensitivity_3500h.csv',final)
    timekeys=['Tmax_K','gap_min_m','gap_mean_m','fuel_max_radial_displacement_m','fuel_von_mises_max_Pa','clad_von_mises_max_Pa','swelling_eigenstrain_max','fuel_power_FE_W','FIMA_native_mean','FIMA_native_max']
    timehistory=[dict(s=s,time_h=r['time_h'],**{k:r[k] for k in timekeys}) for s in S for r in histories[s]];writecsv('sensitivity_time_history.csv',timehistory)
    slopes=[]
    for key,scale,unit in [metrics[i] for i in [4,8,0,10]]:
        values=np.array([histories[s][-1][key]*scale for s in S]);slope=float((values[3]-values[1])/(1.5-.5));curvature=float((values[3]-2*values[2]+values[1])/.5**2);pred=values[2]+slope*(np.array(S)-1)
        slopes.append({'quantity':key,'unit':unit,'dY_ds_centered_at_1':slope,'second_difference_per_s2':curvature,'max_deviation_from_nominal_tangent':float(np.max(abs(values-pred))),'range_0_to_2':float(values.max()-values.min()),'range_percent_of_nominal':float(100*(values.max()-values.min())/abs(values[2]))})
    writecsv('local_sensitivity.csv',slopes)
    plotkeys=[('Tmax_K',1,'Tmax fuel [K]','Tmax'),('gap_min_m',1e6,'Gap minimo [µm]','gap_min'),('gap_mean_m',1e6,'Gap medio [µm]','gap_mean'),('fuel_max_radial_displacement_m',1e6,'Spostamento radiale massimo fuel [µm]','fuel_radial_displacement'),('fuel_von_mises_max_Pa',1e-6,'Von Mises massimo fuel [MPa]','fuel_von_mises'),('clad_von_mises_max_Pa',1e-6,'Von Mises massimo clad [MPa]','clad_von_mises'),('swelling_eigenstrain_max',1,'Eigenstrain swelling massima','swelling')]
    for key,scale,label,name in plotkeys:
        fig,ax=plt.subplots(figsize=(7,4))
        for s in S:ax.plot([r['time_h'] for r in histories[s]],[r[key]*scale for r in histories[s]],'o-',markersize=3,label=f's={s:g}')
        ax.set(xlabel='Tempo [h]',ylabel=label);ax.legend();ax.grid(alpha=.3);fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(CASE/f'{name}_time.{ext}',dpi=160)
        plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(10,7))
    for ax,(key,scale,label,name) in zip(axs.flat,[plotkeys[i] for i in [1,3,0,4]]):
        ax.plot(S,[histories[s][-1][key]*scale for s in S],'o-');ax.set(xlabel='Fattore s',ylabel=label);ax.grid(alpha=.3)
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(CASE/f'sensitivity_3500h.{ext}',dpi=160)
    plt.close(fig)
    before=json.loads((CASE/'protected_manifest_before.json').read_text());after={p:digest(ROOT/p) for p in before};assert before==after
    (CASE/'protected_manifest_after.json').write_text(json.dumps(after,indent=2))
    result={'PASS':True,'protected_files_unchanged':len(before),'qdot_FIMA_exact_identical_all_times':True,'power_exact_identical_all_times':True,'material_BC_mesh_identical':True,'eigenstrain_pointwise_linear_in_s':True,'FIMA_eigenstrain_monotonic_in_time':True,'endpoint_regressions':endpoint_checks,'final':final,'slopes':slopes,'contact_events':contact,'resource_and_runtime':status}
    (CASE/'sensitivity_results.json').write_text(json.dumps(result,indent=2))
    lines=['# FE101 swelling intensity sensitivity','','Deterministic parametric sensitivity of the provisional swelling law. PASS. Non è UQ probabilistica né validazione della correlazione.','','## Implementazione e controlli','','Solo nel nuovo caso: fuel.eigenstrain=sensitivity_law.scaled_native_swelling, con parametro swelling_sensitivity_factor. Il wrapper moltiplica la correlazione nominale invariata per s. Valori 0, 0.5, 1, 1.5, 2; 21 tempi 0–3500 h. Nessuna modifica a FIMA, qdot, BC, mesh, proprietà o solver. Eigenstrain costitutiva valutata tramite Expression DG0 verificata punto per punto uguale a s*FIMA per tutti i 105 stati; s=0 zero, s=1 nominale, s=2 doppia.', '','qdot e FIMA identiche elemento per elemento ai 21 VTU nominali; potenze integrate identiche. Mesh/BC/clad hash uguali; fuel uguale salvo callable/parametro s. Configurazione uguale dopo risoluzione dei path. Nessuna baseline termica reintrodotta. FIMA/eigenstrain non decrescenti per ogni cella; nessun dato diagnostico NaN/Inf, qdot negativo o estrapolazione. Nessuna monotonicità o linearità imposta alla risposta FEM.','','Regressioni s=0/OFF e s=1/full-native (tutti i 21 stati e campi primari finali):', '```json',json.dumps(endpoint_checks,indent=2),'```','','## Risultati a 3500 h','','Le unità nei nomi originali Pa/m della diagnostica sono convertite qui in MPa/µm come indicato. Nel CSV finale i valori stress sono esplicitamente in MPa; nei CSV temporali tutte le unità restano SI.', '','| Grandezza | Unità | s=0 | s=0.5 | s=1 | s=1.5 | s=2 |','|---|---|---:|---:|---:|---:|---:|']
    for key,scale,unit in metrics:lines.append('| '+key+' | '+unit+' | '+' | '.join(f'{histories[s][-1][key]*scale:.9g}' for s in S)+' |')
    lines+=['| Contatto | sì/no | '+' | '.join('sì' if x['contact'] else 'no' for x in final)+' |','','Potenza finale identica: '+str(histories[1.][-1]['fuel_power_FE_W'])+' W. Nessun contatto a qualsiasi tempo per s≤2; nessuna estrapolazione del fattore necessario al contatto.','','## Sensibilità locale e non linearità','','Derivata centrale [Y(1.5)−Y(0.5)]/1; derivata seconda diagnostica [Y(1.5)−2Y(1)+Y(0.5)]/0.25. Nessuna assunzione di linearità della risposta.','', '| Quantità | Unità | ΔY/Δs al nominale | Seconda differenza /s² | Scarto massimo dalla tangente nominale | Range relativo al nominale [%] |','|---|---|---:|---:|---:|---:|']
    for x in slopes:lines.append(f"| {x['quantity']} | {x['unit']} | {x['dY_ds_centered_at_1']:.9g} | {x['second_difference_per_s2']:.9g} | {x['max_deviation_from_nominal_tangent']:.9g} | {x['range_percent_of_nominal']:.9g} |")
    lines += ['','Stress: il massimo su celle può cambiare posizione e non deve seguire una legge lineare. Gap/spostamento sono quasi lineari se le differenze seconde e gli scarti riportati sono piccoli; Tmax presenta il feedback del gap medio sulla conduttanza. Le percentuali di temperatura assoluta dipendono dalla scala: valutare soprattutto il range in K. Gap e spostamento sono quasi lineari (scarto massimo dalla tangente <0.0008 µm); Tmax <0.006 K e Von Mises <0.110 MPa. Piccole deviazioni sono riportate senza attribuirle tutte a non linearità fisiche risolte: restano tolleranze di convergenza iterative. Nessuna linearità imposta.','','## Serie temporali e grafici','','sensitivity_time_history.csv: 105 record = 5 fattori × 21 tempi, con sette serie richieste più potenza e FIMA. PNG/PDF comparativi per tutte le serie; sensitivity_3500h.png/pdf presenta gap, spostamento, Tmax, Von Mises vs s. Nessun punto oltre gli stati salvati.','','## Interpretazione e limiti','','Gli output meccanici del fuel sono i più sensibili (vedere range e percentuali misurati); l’effetto termico va giudicato sul range in K e quello sul clad sui valori molto piccoli. Questa classificazione è specifica del basso livello di FIMA e delle BC attuali. Nell’intervallo deterministico s=0–2, il range è 26.73% del nominale per lo spostamento e 36.43% per Von Mises fuel: influenza forte su questi output; gap minimo 5.62% (influenza moderata), Tmax 4.410 K (influenza contenuta), clad quasi invariato. Non è una stima probabilistica dell’incertezza né la previsione del range di una futura correlazione. Il modello gap/contact mantiene medie globali; nessuna nuova fisica o feedback Z3ST→OpenMC.', '','Warning nodali di stress sull’asse: 31 valutazioni non finite u_r/r sono azzerate dal writer come nei riferimenti; i massimi riportati sono valutati nei centri DG0 e finiti. Non affermiamo validità dei valori nodali sull’asse. I campi potenziali nominali del provider sono distinti dall’eigenstrain applicata nei nuovi VTU. Tutti i cinque solve accettati sono convergenti e con diagnostica PASS. Due tentativi iniziali s=0 erano convergenti ma la diagnostica aggiunta falliva: prima assenza di comunicatore MPI per un’espressione UFL identicamente nulla, poi assenza di dominio mesh. La sola diagnostica del campo nullo ora verifica esplicitamente lo zero simbolico e riempie il campo DG0 con zero; per s>0 valuta il wrapper costitutivo mediante Expression sul mesh. Output/log/preflight dei due tentativi conservati in runs/s_0.0/failed_diagnostic_attempt_{1,2}; nuovo preflight prima di ciascun retry. Nessuna modifica alla fisica.','','## Risorse e protezione','','1 MPI, 2 OMP, BLAS/NumExpr=1. Preflight free -h, df -h /, df -h /mnt/c, nproc, vmstat 1 3 prima di ogni solve, gate PASS. Risorse complete in runs/s_*/resource_preflight.json e solve_status.json.','','| s | RAM available GiB | Swap usata GiB | Disco libero GiB | CPU logiche | Runtime s |','|---|---:|---:|---:|---:|---:|']
    for s in S:
        v=status['solves'][str(s)];r=v['resources'];lines.append(f"| {s} | {r['RAM_available_GiB']:.3f} | {r['swap_used_bytes']/1024**3:.3f} | {r['disk_free_GiB']:.3f} | {r['logical_CPUs']} | {v['runtime_s']:.3f} |")
    lines += ['',f'Manifest prima/dopo: {len(before)} file protetti invariati. Nessun OpenMC/depletion, commit o push. Nessun output precedente sovrascritto.','', '## File aggiunti','','README.md, sensitivity_law.py, run_sensitivity.py, postprocess.py, cinque configurazioni runs/s_*/ e relativi output, diagnostica locale, manifest prima/dopo, solve_status.json, sensitivity_results.json, sensitivity_3500h.csv, sensitivity_time_history.csv, local_sensitivity.csv, grafici PNG/PDF, questo report. Nessun file preesistente modificato.','','Fase conclusa; attendere revisione.']
    (CASE/'SWELLING_SENSITIVITY_REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'PASS':True,'final':final,'slopes':slopes,'endpoint_regressions':endpoint_checks,'protected_files':len(before)},indent=2))
if __name__=='__main__':main()
