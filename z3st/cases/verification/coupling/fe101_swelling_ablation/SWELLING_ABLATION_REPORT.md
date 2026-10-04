# FE101 swelling ablation — numerical / model-effect isolation

PASS. Isolamento dell’effetto della correlazione provvisoria, non validazione fisica dello swelling.

## Configurazione e protezione

Nuovo caso OFF con identici FIMA, qdot, mesh, BC, proprietà materiali, impostazioni e 21 tempi 0–3500 h. Solo rimosso il richiamo `eigenstrain` nella card fuel isolata. FIMA resta importata; nessun dato sorgente è azzerato. Nessuna baseline termica reintrodotta. Nessuna modifica ai sorgenti Z3ST, neutronica, regression o FASE 1/2. Nessun OpenMC/depletion, commit/push.
Manifest prima/dopo: 1500 file protetti identici.

## Risultati a 3500 h

| Grandezza | Unità | ON | OFF | ON − OFF | Differenza % rispetto a OFF |
|---|---|---:|---:|---:|---:|
| Tmax_K | K | 483.862787 | 486.063609 | -2.20082138 | — |
| fuel_surface_T_min_K | K | 367.269984 | 368.624802 | -1.35481819 | — |
| fuel_surface_T_mean_K | K | 393.950893 | 395.809526 | -1.85863301 | — |
| fuel_surface_T_max_K | K | 415.62951 | 417.90577 | -2.27625912 | — |
| gap_min_m | um | 109.912063 | 113.00209 | -3.09002637 | -2.734486042212268 |
| gap_mean_m | um | 114.003927 | 116.515279 | -2.51135173 | -2.155384050846547 |
| gap_max_m | um | 119.171057 | 120.969623 | -1.79856631 | -1.4867916964533725 |
| gap_min_position_local_z_m | cm | 18.9653333 | 18.9653333 | 0 | — |
| fuel_max_radial_displacement_m | um | 23.1186055 | 20.0285643 | 3.09004122 | 15.428171351999884 |
| clad_inner_ur_mean_m | um | 3.00924322 | 3.00924068 | 2.53951344e-06 | 8.439050618318116e-05 |
| fuel_von_mises_max_Pa | MPa | 26.0543799 | 30.8852236 | -4.83084369 | -15.641278031515107 |
| clad_von_mises_max_Pa | MPa | 0.554436831 | 0.554370584 | 6.62473697e-05 | 0.011950015317010676 |
| fuel_power_FE_W | W | 5265.81446 | 5265.81446 | 0 | 0.0 |
| model_mean_gap_m | um | 114.003893 | 116.515219 | -2.51132659 | -2.1553635671708946 |

Percentuali omesse per temperature assolute e posizioni, perché dipendono dall’origine/unità. Contatto assente in entrambi i casi a tutti i tempi.

## Verifiche e storia temporale

qdot DG0 e FIMA: uguaglianza esatta elemento per elemento dei 21 VTU ON/OFF. Potenze integrate identiche. Mesh, BC, clad: hash identici. Fuel: identico salvo callable di swelling. Configurazioni: identiche dopo risoluzione dei path. FIMA cumulativa non decrescente in ogni cella; eigenstrain ON monotona, OFF identicamente zero. Nessun qdot negativo, gap negativo, valore diagnostico NaN/Inf o estrapolazione.

`time_history.csv` contiene i 21 tempi con valori ON/OFF e differenze per tutte le serie richieste; `comparison_3500h.csv` contiene il confronto finale. Grafici PNG/PDF: Tmax, gap minimo, spostamento radiale, Von Mises fuel, FIMA/eigenstrain.

## Interpretazione e limiti

A potenza identica, lo swelling aumenta lo spostamento radiale e riduce il gap. Il massimo Von Mises del fuel diminuisce: lo stress non deve essere reso monotono artificialmente. La temperatura diminuisce attraverso la risposta del gap medio e della conduttanza già esistente; nessuna modifica al dominio/mesh di riferimento o alla sorgente. Il gap/contact usa ancora quantità mediate globalmente. Non attribuire questo risultato a feedback neutronico.

La FIMA massima finale è 0.00022399710727541027; εsw massima ON uguale alla FIMA, OFF zero. L’output nativo del provider rappresenta il coefficiente potenziale: nei soli VTU OFF è rinominato `Swelling_potential_eigenstrain_native`, con un campo separato `Swelling_eigenstrain_applied=0`. La diagnostica meccanica usa lo stress costitutivo effettivo.

Warning esistente: gli output nodali di stress sull’asse contengono 31 valutazioni non finite dovute a u_r/r, sostituite dal writer. I massimi riportati sono valutati nei centri DG0 delle celle, finiti; non sono i massimi nodali corretti artificialmente. Nessun fallimento del solve/diagnostica. Correlazione U-ZrH provvisoria; nessuna nuova fisica, proprietà o sensitivity.

## Risorse e riproducibilità

Runtime solve: 14.774 s; MPI=1, OMP=2, BLAS/NumExpr=1. RAM available 8.279 GiB; swap usata 0.692 GiB; disco libero 923.020 GiB; CPU logiche 8. Gate PASS. Comandi richiesti e vmstat completi in runs/swelling_off/resource_preflight.json.

Un warning non bloccante della cache Matplotlib durante l’import del runner ha usato /tmp; il solve usa cache nel nuovo caso. Nessuna modifica dell’ambiente fisico.

## File aggiunti

README.md; run_ablation.py; postprocess.py; runs/swelling_off/{input,fuel,clad,geometry,boundary_conditions}.yaml; mesh.msh; diagnostics.py; output isolati; manifest prima/dopo; solve_status.json; ablation_results.json; time_history.csv; comparison_3500h.csv; grafici PNG/PDF; questo report. Nessun file sorgente modificato.

Analisi terminata; attendere revisione.
