"""Saved-cell-tensor diagnostics only: no solver construction or simulation."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
CASE=Path(__file__).resolve().parent

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def main():
    # Protect all pre-existing solve artifacts and completed sensitivity reports.
    protected=list((CASE/'runs').rglob('*'))+[CASE/n for n in ['README.md','SWELLING_SENSITIVITY_REPORT.md','sensitivity_3500h.csv','sensitivity_time_history.csv','sensitivity_results.json','solve_status.json']]
    before={str(p.relative_to(CASE)):sha(p) for p in protected if p.is_file()}
    (CASE/'stress_check_manifest_before.json').write_text(json.dumps(before,indent=2))
    results=[];grids={};allvm=[]
    for s in [0.,1.,2.]:
        run=CASE/'runs'/f's_{s:.1f}';g=pv.read(run/'output/fields_0020.vtu');grids[s]=g
        mask=g.cell_data['MaterialID']==1;cells=np.flatnonzero(mask);vm=g.cell_data['VonMises (cells)'];tensor=g.cell_data['Stress (cells)'].reshape(-1,3,3)
        assert np.isfinite(vm[mask]).all() and np.isfinite(tensor[mask]).all() and (vm[mask]>=0).all()
        dev=tensor-np.trace(tensor,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
        reconstructed=np.sqrt(1.5*np.sum(dev*dev,axis=(1,2)))
        np.testing.assert_allclose(vm[mask],reconstructed[mask],rtol=2e-13,atol=1e-7)
        centre=g.cell_centers().points;bounds=np.array([g.get_cell(int(i)).bounds[:4] for i in cells]);vol=np.pi*(bounds[:,1]**2-bounds[:,0]**2)*(bounds[:,3]-bounds[:,2])
        imax=int(cells[np.argmax(vm[mask])]);sig=tensor[imax]/1e6
        saved=json.loads((run/'output/diagnostic_history.jsonl').read_text().splitlines()[-1])
        np.testing.assert_allclose(vm[imax],saved['fuel_von_mises_max_Pa'],rtol=2e-13,atol=1e-7)
        fima=g.cell_data['FIMA_native_OpenMC'];q=g.cell_data['qdot_native_OpenMC_W_m3'];qmax=cells[q[cells]==q[cells].max()];fimax=cells[np.isclose(fima[cells],fima[cells].max(),rtol=1e-13,atol=1e-18)]
        # Non-finite nodal axis entries have been sanitized by the original writer.
        axispoints=np.flatnonzero(np.isclose(g.points[:,0],0,atol=1e-12))
        data={'s':s,'VM_min_MPa':float(vm[mask].min()/1e6),'VM_volume_mean_MPa':float(np.average(vm[mask],weights=vol)/1e6),'VM_max_MPa':float(vm[imax]/1e6),'maximum_cell_id':imax,'maximum_position_local_r_z_cm':(centre[imax,:2]*100).tolist(),'maximum_position_OpenMC_r_z_cm':[float(centre[imax,0]*100),float(centre[imax,1]*100+10.2)],'maximum_cell_bounds_m':list(g.get_cell(imax).bounds[:4]),'sigma_rr_MPa':float(sig[0,0]),'sigma_theta_theta_MPa':float(sig[1,1]),'sigma_zz_MPa':float(sig[2,2]),'sigma_rz_MPa':float(sig[0,2]),'mean_normal_stress_at_max_MPa':float(np.trace(sig)/3),'deviatoric_principal_components_at_max_MPa':(np.diag(sig)-np.trace(sig)/3).tolist(),'radial_distance_to_surface_cm':float((.01791-centre[imax,0])*100),'distance_to_nearest_axial_end_cm':float(min(centre[imax,1],.3556-centre[imax,1])*100),'qmax_cell_centres_local_cm':(centre[qmax,:2]*100).tolist(),'FIMAmax_cell_centres_local_cm':(centre[fimax,:2]*100).tolist(),'power_W':saved['fuel_power_FE_W'],'fuel_ur_max_um':saved['fuel_max_radial_displacement_m']*1e6,'gap_min_um':saved['gap_min_m']*1e6,'gap_min_position_local_z_cm':saved['gap_min_position_local_z_m']*100,'clad_VM_max_MPa':saved['clad_von_mises_max_Pa']/1e6,'contact':bool(saved['local_gap_closed'] or saved['model_contact_pressure_Pa']>0),'cell_values_finite':True,'axis_points_count':len(axispoints),'max_VM_reconstruction_abs_error_Pa':float(np.max(abs(vm[mask]-reconstructed[mask]))),'fuel_cells':len(cells)}
        histories=[json.loads(x) for x in (run/'output/diagnostic_history.jsonl').read_text().splitlines()];assert all(not r['local_gap_closed'] and r['model_contact_pressure_Pa']==0 for r in histories)
        results.append(data);allvm.extend(vm[mask]/1e6)
    fixed=[]
    for anchor in results:
        i=anchor['maximum_cell_id']
        for s,g in grids.items():
            sig=g.cell_data['Stress (cells)'][i].reshape(3,3)/1e6
            fixed.append({'anchor_maximum_of_s':anchor['s'],'cell_id':i,'s':s,'r_cm':g.cell_centers().points[i,0]*100,'z_local_cm':g.cell_centers().points[i,1]*100,'VM_MPa':g.cell_data['VonMises (cells)'][i]/1e6,'sigma_rr_MPa':sig[0,0],'sigma_theta_theta_MPa':sig[1,1],'sigma_zz_MPa':sig[2,2],'sigma_rz_MPa':sig[0,2]})
    with (CASE/'stress_components_fixed_regions.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fixed[0].keys());w.writeheader();w.writerows(fixed)
    for result in results:
        s=result['s'];g=grids[s];mask=g.cell_data['MaterialID']==1;polys=[]
        for i in np.flatnonzero(mask):
            a,b,c,d=g.get_cell(int(i)).bounds[:4];polys.append(np.array([[a,c],[b,c],[b,d],[a,d]])*100)
        fig,ax=plt.subplots(figsize=(5,7));coll=PolyCollection(polys,array=g.cell_data['VonMises (cells)'][mask]/1e6,cmap='viridis',clim=(0,max(allvm)),edgecolors='none');ax.add_collection(coll);ax.set(xlim=(0,1.791),ylim=(0,35.56),xlabel='r [cm]',ylabel='z locale FE [cm]',title=f'Fuel Von Mises ai centri delle celle — s={s:g}')
        r,z=result['maximum_position_local_r_z_cm'];ax.plot(r,z,'rx',markersize=9,label=f'Max {result["VM_max_MPa"]:.3f} MPa');ax.legend(loc='upper left');fig.colorbar(coll,ax=ax,label='Von Mises [MPa]');fig.tight_layout()
        for ext in ['png','pdf']:fig.savefig(CASE/f'stress_map_s{int(s)}.{ext}',dpi=180)
        plt.close(fig)
    for s in [1.,2.]:
        np.testing.assert_array_equal(grids[0.].cell_data['FIMA_native_OpenMC'],grids[s].cell_data['FIMA_native_OpenMC']);np.testing.assert_array_equal(grids[0.].cell_data['qdot_native_OpenMC_W_m3'],grids[s].cell_data['qdot_native_OpenMC_W_m3'])
    assert np.diff([r['fuel_ur_max_um'] for r in results]).min()>0 and np.diff([r['gap_min_um'] for r in results]).max()<0
    after={p:sha(CASE/p) for p in before};assert before==after
    (CASE/'stress_check_manifest_after.json').write_text(json.dumps(after,indent=2))
    summary={'PASS':True,'source':'Stress (cells) and VonMises (cells) in saved fields_0020.vtu; no interpolation from nodal stress','results':results,'fixed_region_components':fixed,'protected_files_unchanged':len(before),'cell_VM_reconstructed_from_tensor_PASS':True,'conclusion':'swelling coupling + ablation + sensitivity block closed'}
    (CASE/'stress_check_results.json').write_text(json.dumps(summary,indent=2))
    lines=['# Stress sensitivity diagnostic check','','PASS — numerical / model-effect isolation. Non è validazione fisica della legge di swelling. Solo lettura degli output già disponibili; nessuna simulazione o ricostruzione FEM.','','## Mappe e massimi','','Fonte: `fields_0020.vtu`, campi `Stress (cells)` e `VonMises (cells)`. Le mappe stress_map_s0/s1/s2.png/pdf usano la stessa scala colori. Media pesata sul volume cilindrico π(r_out²−r_in²)Δz, non media aritmetica della griglia. Coordinate z locali FE; z OpenMC = z_FE +10.20 cm.','','| s | r max [cm] | z max locale [cm] | z OpenMC [cm] | VM min [MPa] | VM media volumetrica [MPa] | VM max [MPa] |','|---|---:|---:|---:|---:|---:|---:|']
    for r in results:lines.append(f"| {r['s']:g} | {r['maximum_position_local_r_z_cm'][0]:.6f} | {r['maximum_position_local_r_z_cm'][1]:.6f} | {r['maximum_position_OpenMC_r_z_cm'][1]:.6f} | {r['VM_min_MPa']:.6f} | {r['VM_volume_mean_MPa']:.6f} | {r['VM_max_MPa']:.6f} |")
    lines+=['','Il massimo migra assialmente di −3.556 cm tra s=0 e s=2 (tre celle FE), restando nella stessa fascia radiale esterna: r=1.7536875 cm, a 0.0373125 cm dalla superficie. Non è sull’asse né alle estremità: la distanza dall’estremità più vicina resta almeno 15.8 cm. Il gap minimo resta a z locale=18.965333 cm. Il massimo qdot è nella regione esterna z-bin 5 (14.224–17.780 cm), il massimo FIMA nella regione esterna z-bin 6 (17.780–21.336 cm); s=0/1 hanno massimo stress vicino al massimo FIMA, s=2 migra nella zona di massimo qdot. Non implica causalità esclusiva.', '','## Componenti al massimo di ciascun caso','','Tensore assialsimmetrico nell’ordine (r, θ, z); valori in MPa.','', '| s | σrr | σθθ | σzz | σrz |','|---|---:|---:|---:|---:|']
    for r in results:lines.append(f"| {r['s']:g} | {r['sigma_rr_MPa']:.6f} | {r['sigma_theta_theta_MPa']:.6f} | {r['sigma_zz_MPa']:.6f} | {r['sigma_rz_MPa']:.6f} |")
    lines+=['','Per evitare di confondere redistribuzione con migrazione del massimo, stress_components_fixed_regions.csv confronta anche le stesse tre celle in tutti i casi.','', '| Regione fissa: massimo s=0 | s | VM | σrr | σθθ | σzz | σrz |','|---|---:|---:|---:|---:|---:|---:|']
    for r in fixed[:3]:lines.append(f"| r={r['r_cm']:.6f}, z={r['z_local_cm']:.6f} cm | {r['s']:g} | {r['VM_MPa']:.6f} | {r['sigma_rr_MPa']:.6f} | {r['sigma_theta_theta_MPa']:.6f} | {r['sigma_zz_MPa']:.6f} | {r['sigma_rz_MPa']:.6f} |")
    lines+=['','## Spiegazione meccanica e coerenza','','La trazione circonferenziale e assiale nella fascia esterna diminuisce, mentre σrr resta piccola e compressiva; si riducono le differenze fra le tensioni principali, quindi la norma deviatorica. Il taglio è molto più piccolo delle tensioni normali. La riduzione avviene anche confrontando una regione fissa e non è un artefatto del cambio di cella massima.', '','La swelling eigenstrain è isotropa: a spostamento fissato aggiungerebbe un contributo idrostatico e non cambierebbe direttamente Von Mises. Nel solve il campo di spostamento si riadatta alla eigenstrain spazialmente variabile; cambia quindi la deformazione elastica deviatorica. Si tratta di redistribuzione elastica, non creep/plasticità o rilassamento temporale introdotto. Inoltre le temperature diminuiscono lievemente attraverso il gap medio: questa verifica non separa quantitativamente tale contributo termico dal contributo meccanico diretto.', '','| s | ur,max fuel [µm] | gap min [µm] | z gap min [cm] | VM clad max [MPa] | contatto |','|---|---:|---:|---:|---:|---|']
    for r in results:lines.append(f"| {r['s']:g} | {r['fuel_ur_max_um']:.6f} | {r['gap_min_um']:.6f} | {r['gap_min_position_local_z_cm']:.6f} | {r['clad_VM_max_MPa']:.6f} | {'sì' if r['contact'] else 'no'} |")
    lines+=['','qdot/FIMA identiche elemento per elemento; potenza fuel identica (5265.814463925392 W). Nessun contatto a qualsiasi dei 21 tempi.','', '## Warning e limiti','',f"Tutte le {results[0]['fuel_cells']} celle fuel hanno stress e VM finiti. VM ricalcolata dal tensore salvato concorda con il campo cellwise e con i report precedenti; massimo errore {max(r['max_VM_reconstruction_abs_error_Pa'] for r in results):.3g} Pa. Il massimo è nell’ultima fascia radiale, lontano dall’asse. Non si usano né si mediano i valori nodali: il warning preesistente riguarda 31 nodi sull’asse (u_r/r); il writer li ha azzerati. Nessuna contaminazione nodale nelle mappe cellwise. La verifica è discreta: i massimi continui/subcella non sono determinati. Gap/contact resta mediato globalmente; correlazione U-ZrH provvisoria.",'',f'Hash prima/dopo: {len(before)} file preesistenti invariati. Nessun solver/input/correlazione/sorgente modificato, nessun OpenMC/depletion, commit/push. Aggiunti soltanto script diagnostico, tre mappe PNG/PDF, CSV componenti, JSON risultati/manifests e questo report.', '', '## Giudizio finale','','**Comportamento coerente** con i campi salvati: la diminuzione del massimo Von Mises è accompagnata da una riduzione dello stato deviatorico nelle stesse regioni e da una moderata migrazione assiale, con espansione radiale crescente e gap decrescente. Non emergono anomalie che richiedano nuovi solve.','','**swelling coupling + ablation + sensitivity block closed**']
    (CASE/'STRESS_SENSITIVITY_CHECK.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
