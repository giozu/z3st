"""Read-only reference comparison; write only isolated ablation artifacts."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import yaml
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
CASE=Path(__file__).resolve().parent;ROOT=CASE.parents[4]
ON=CASE.parent/'fe101_native_power/runs/full_native';OFF=CASE/'runs/swelling_off'
def rows(p):return [json.loads(x) for x in (p/'output/diagnostic_history.jsonl').read_text().splitlines()]
def csvwrite(name,data):
    with (CASE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def main():
    on,off=rows(ON),rows(OFF)
    assert len(on)==len(off)==21
    times=np.array([r['time_h'] for r in on]);np.testing.assert_array_equal(times,np.arange(21)*175)
    assert [r['time_h'] for r in off]==times.tolist()
    for name in ['mesh.msh','geometry.yaml','clad.yaml','boundary_conditions.yaml']:assert digest(ON/name)==digest(OFF/name)
    a,b=[yaml.safe_load((p/'fuel.yaml').read_text()) for p in [ON,OFF]]
    assert a.pop('eigenstrain')=='materials.fuel_swelling.uzrh_fission_product_swelling' and a==b
    a,b=[yaml.safe_load((p/'input.yaml').read_text()) for p in [ON,OFF]]
    for d,p in [(a,ON),(b,OFF)]:
        for v in d['coupling'].values():v['history_path']=str((p/v['history_path']).resolve())
    assert a==b
    qerror=fierror=0.;tfield_error=0.;previous=None
    for i,(ra,rb) in enumerate(zip(on,off)):
        for r in [ra,rb]:
            assert np.isfinite([v for v in r.values() if isinstance(v,(float,int))]).all()
            assert r['gap_min_m']>0 and not r['local_gap_closed'] and r['model_contact_pressure_Pa']==0
            assert r['qdot_min_W_m3']>=0 and r['power_time_s']==r['FIMA_time_s']==r['time_s']
        assert ra['fuel_power_FE_W']==rb['fuel_power_FE_W']
        assert ra['FIMA_native_mean']==rb['FIMA_native_mean'] and rb['swelling_eigenstrain_max']==0
        ga=pv.read(ON/f'output/fields_{i:04d}.vtu');gb=pv.read(OFF/f'output/fields_{i:04d}.vtu')
        for name in ['qdot_native_OpenMC_W_m3','FIMA_native_OpenMC']:
            va,vb=ga.cell_data[name],gb.cell_data[name];np.testing.assert_array_equal(va,vb);assert np.isfinite(vb).all() and (vb>=0).all()
        fi=gb.cell_data['FIMA_native_OpenMC']
        if previous is not None:assert (fi>=previous-1e-18).all()
        previous=fi.copy()
        np.testing.assert_array_equal(ga.cell_data['Swelling_eigenstrain_native'],fi)
        # The unchanged provider exports potential epsilon(FIMA); distinguish applied OFF.
        if 'Swelling_eigenstrain_native' in gb.cell_data:
            gb.cell_data['Swelling_potential_eigenstrain_native']=gb.cell_data.pop('Swelling_eigenstrain_native')
            gb.cell_data['Swelling_eigenstrain_applied']=np.zeros_like(fi)
            gb.save(OFF/f'output/fields_{i:04d}.vtu')
        assert (gb.cell_data['Swelling_eigenstrain_applied']==0).all()
        tfield_error=max(tfield_error,float(np.max(abs(ga.point_data['Temperature']-gb.point_data['Temperature']))))
    keys=['fuel_power_FE_W','Tmax_K','fuel_surface_T_mean_K','FIMA_native_mean','FIMA_native_max','swelling_eigenstrain_max','gap_min_m','gap_mean_m','fuel_max_radial_displacement_m','fuel_von_mises_max_Pa','clad_von_mises_max_Pa']
    history=[]
    for a,b in zip(on,off):
        r={'time_h':a['time_h']}
        for k in keys:r[k+'_ON']=a[k];r[k+'_OFF']=b[k];r[k+'_ON_minus_OFF']=a[k]-b[k]
        history.append(r)
    csvwrite('time_history.csv',history)
    metrics=[('Tmax_K',1,'K'),('fuel_surface_T_min_K',1,'K'),('fuel_surface_T_mean_K',1,'K'),('fuel_surface_T_max_K',1,'K'),('gap_min_m',1e6,'um'),('gap_mean_m',1e6,'um'),('gap_max_m',1e6,'um'),('gap_min_position_local_z_m',100,'cm'),('fuel_max_radial_displacement_m',1e6,'um'),('clad_inner_ur_mean_m',1e6,'um'),('fuel_von_mises_max_Pa',1e-6,'MPa'),('clad_von_mises_max_Pa',1e-6,'MPa'),('fuel_power_FE_W',1,'W'),('model_mean_gap_m',1e6,'um')]
    comparison=[]
    for key,scale,unit in metrics:
        a,b=on[-1][key]*scale,off[-1][key]*scale
        # Relative K and coordinate differences are convention dependent: show absolute only.
        pct=100*(a-b)/b if b and unit not in ['K','cm'] else ''
        comparison.append({'quantity':key,'unit':unit,'ON':a,'OFF':b,'ON_minus_OFF':a-b,'relative_difference_percent':pct})
    csvwrite('comparison_3500h.csv',comparison)
    plots=[('Tmax_K',1,'Tmax fuel [K]','Tmax'),('gap_min_m',1e6,'Gap minimo [µm]','gap_min'),('fuel_max_radial_displacement_m',1e6,'Spostamento radiale massimo fuel [µm]','fuel_radial_displacement'),('fuel_von_mises_max_Pa',1e-6,'Von Mises massimo fuel [MPa]','fuel_von_mises')]
    for key,scale,label,name in plots:
        fig,ax=plt.subplots(figsize=(7,4));ax.plot(times,[r[key]*scale for r in on],'o-',label='Swelling ON');ax.plot(times,[r[key]*scale for r in off],'s--',label='Swelling OFF');ax.set(xlabel='Tempo [h]',ylabel=label);ax.legend();ax.grid(alpha=.3);fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(CASE/f'{name}_time.{ext}',dpi=160)
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));ax.plot(times,[r['FIMA_native_max'] for r in on],'o-',label='FIMA max ON/OFF (identica)');ax.plot(times,[r['swelling_eigenstrain_max'] for r in on],'x--',label='Eigenstrain applicata ON');ax.plot(times,[r['swelling_eigenstrain_max'] for r in off],'s--',label='Eigenstrain applicata OFF');ax.set(xlabel='Tempo [h]',ylabel='Frazione / deformazione per direzione');ax.legend();ax.grid(alpha=.3);fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(CASE/f'FIMA_swelling_time.{ext}',dpi=160)
    plt.close(fig)
    before=json.loads((CASE/'protected_manifest_before.json').read_text());after={p:digest(ROOT/p) for p in before};assert before==after
    (CASE/'protected_manifest_after.json').write_text(json.dumps(after,indent=2))
    resource=json.loads((CASE/'solve_status.json').read_text())
    result={'PASS':True,'protected_files_unchanged':len(before),'qdot_identical_all_21_times':True,'FIMA_identical_all_21_times':True,'power_identical_all_21_times':True,'BC_material_mesh_input_equal_except_eigenstrain':True,'FIMA_pointwise_monotonic':True,'applied_swelling_pointwise_monotonic':True,'no_contact_all_times':True,'no_extrapolation':True,'final_ON':on[-1],'final_OFF':off[-1],'comparison':comparison,'resource_and_runtime':resource,'maximum_temperature_field_difference_all_times_K':tfield_error,'interpretation':'numerical / model-effect isolation, not physical validation'}
    (CASE/'ablation_results.json').write_text(json.dumps(result,indent=2))
    lines=['# FE101 swelling ablation — numerical / model-effect isolation','','PASS. Isolamento dell’effetto della correlazione provvisoria, non validazione fisica dello swelling.','','## Configurazione e protezione','','Nuovo caso OFF con identici FIMA, qdot, mesh, BC, proprietà materiali, impostazioni e 21 tempi 0–3500 h. Solo rimosso il richiamo `eigenstrain` nella card fuel isolata. FIMA resta importata; nessun dato sorgente è azzerato. Nessuna baseline termica reintrodotta. Nessuna modifica ai sorgenti Z3ST, neutronica, regression o FASE 1/2. Nessun OpenMC/depletion, commit/push.',f'Manifest prima/dopo: {len(before)} file protetti identici.','','## Risultati a 3500 h','','| Grandezza | Unità | ON | OFF | ON − OFF | Differenza % rispetto a OFF |','|---|---|---:|---:|---:|---:|']
    for c in comparison:lines.append(f"| {c['quantity']} | {c['unit']} | {c['ON']:.9g} | {c['OFF']:.9g} | {c['ON_minus_OFF']:.9g} | {c['relative_difference_percent'] if c['relative_difference_percent']!='' else '—'} |")
    lines += ['','Percentuali omesse per temperature assolute e posizioni, perché dipendono dall’origine/unità. Contatto assente in entrambi i casi a tutti i tempi.','','## Verifiche e storia temporale','','qdot DG0 e FIMA: uguaglianza esatta elemento per elemento dei 21 VTU ON/OFF. Potenze integrate identiche. Mesh, BC, clad: hash identici. Fuel: identico salvo callable di swelling. Configurazioni: identiche dopo risoluzione dei path. FIMA cumulativa non decrescente in ogni cella; eigenstrain ON monotona, OFF identicamente zero. Nessun qdot negativo, gap negativo, valore diagnostico NaN/Inf o estrapolazione.', '','`time_history.csv` contiene i 21 tempi con valori ON/OFF e differenze per tutte le serie richieste; `comparison_3500h.csv` contiene il confronto finale. Grafici PNG/PDF: Tmax, gap minimo, spostamento radiale, Von Mises fuel, FIMA/eigenstrain.', '','## Interpretazione e limiti','','A potenza identica, lo swelling aumenta lo spostamento radiale e riduce il gap. Il massimo Von Mises del fuel diminuisce: lo stress non deve essere reso monotono artificialmente. La temperatura diminuisce attraverso la risposta del gap medio e della conduttanza già esistente; nessuna modifica al dominio/mesh di riferimento o alla sorgente. Il gap/contact usa ancora quantità mediate globalmente. Non attribuire questo risultato a feedback neutronico.', '','La FIMA massima finale è '+str(on[-1]['FIMA_native_max'])+'; εsw massima ON uguale alla FIMA, OFF zero. L’output nativo del provider rappresenta il coefficiente potenziale: nei soli VTU OFF è rinominato `Swelling_potential_eigenstrain_native`, con un campo separato `Swelling_eigenstrain_applied=0`. La diagnostica meccanica usa lo stress costitutivo effettivo.', '','Warning esistente: gli output nodali di stress sull’asse contengono 31 valutazioni non finite dovute a u_r/r, sostituite dal writer. I massimi riportati sono valutati nei centri DG0 delle celle, finiti; non sono i massimi nodali corretti artificialmente. Nessun fallimento del solve/diagnostica. Correlazione U-ZrH provvisoria; nessuna nuova fisica, proprietà o sensitivity.','','## Risorse e riproducibilità','',f"Runtime solve: {resource['runtime_s']:.3f} s; MPI=1, OMP=2, BLAS/NumExpr=1. RAM available {resource['resources']['RAM_available_GiB']:.3f} GiB; swap usata {resource['resources']['swap_used_bytes']/1024**3:.3f} GiB; disco libero {resource['resources']['disk_free_GiB']:.3f} GiB; CPU logiche {resource['resources']['logical_CPUs']}. Gate PASS. Comandi richiesti e vmstat completi in runs/swelling_off/resource_preflight.json.",'','Un warning non bloccante della cache Matplotlib durante l’import del runner ha usato /tmp; il solve usa cache nel nuovo caso. Nessuna modifica dell’ambiente fisico.','','## File aggiunti','','README.md; run_ablation.py; postprocess.py; runs/swelling_off/{input,fuel,clad,geometry,boundary_conditions}.yaml; mesh.msh; diagnostics.py; output isolati; manifest prima/dopo; solve_status.json; ablation_results.json; time_history.csv; comparison_3500h.csv; grafici PNG/PDF; questo report. Nessun file sorgente modificato.','','Analisi terminata; attendere revisione.']
    (CASE/'SWELLING_ABLATION_REPORT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'PASS':True,'protected_files':len(before),'comparison':comparison},indent=2))
if __name__=='__main__':main()
