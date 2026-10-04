# FE101 swelling intensity sensitivity

Deterministic parametric sensitivity of the provisional swelling law. PASS. Non è UQ probabilistica né validazione della correlazione.

## Implementazione e controlli

Solo nel nuovo caso: fuel.eigenstrain=sensitivity_law.scaled_native_swelling, con parametro swelling_sensitivity_factor. Il wrapper moltiplica la correlazione nominale invariata per s. Valori 0, 0.5, 1, 1.5, 2; 21 tempi 0–3500 h. Nessuna modifica a FIMA, qdot, BC, mesh, proprietà o solver. Eigenstrain costitutiva valutata tramite Expression DG0 verificata punto per punto uguale a s*FIMA per tutti i 105 stati; s=0 zero, s=1 nominale, s=2 doppia.

qdot e FIMA identiche elemento per elemento ai 21 VTU nominali; potenze integrate identiche. Mesh/BC/clad hash uguali; fuel uguale salvo callable/parametro s. Configurazione uguale dopo risoluzione dei path. Nessuna baseline termica reintrodotta. FIMA/eigenstrain non decrescenti per ogni cella; nessun dato diagnostico NaN/Inf, qdot negativo o estrapolazione. Nessuna monotonicità o linearità imposta alla risposta FEM.

Regressioni s=0/OFF e s=1/full-native (tutti i 21 stati e campi primari finali):
```json
{
  "0.0": {
    "PASS": true,
    "max_temperature_difference_K": 1.3073986337985843e-11,
    "max_displacement_difference_m": 2.710505431213761e-18,
    "reference": "/home/simone/z3st/z3st/cases/verification/coupling/fe101_swelling_ablation/runs/swelling_off"
  },
  "1.0": {
    "PASS": true,
    "max_temperature_difference_K": 1.1596057447604835e-11,
    "max_displacement_difference_m": 5.5294310796760726e-18,
    "reference": "/home/simone/z3st/z3st/cases/verification/coupling/fe101_native_power/runs/full_native"
  }
}
```

## Risultati a 3500 h

Le unità nei nomi originali Pa/m della diagnostica sono convertite qui in MPa/µm come indicato. Nel CSV finale i valori stress sono esplicitamente in MPa; nei CSV temporali tutte le unità restano SI.

| Grandezza | Unità | s=0 | s=0.5 | s=1 | s=1.5 | s=2 |
|---|---|---:|---:|---:|---:|---:|
| Tmax_K | K | 486.063609 | 484.964184 | 483.862787 | 482.757815 | 481.653364 |
| fuel_surface_T_min_K | K | 368.624802 | 367.947844 | 367.269984 | 366.590286 | 365.911214 |
| fuel_surface_T_mean_K | K | 395.809526 | 394.881002 | 393.950893 | 393.017814 | 392.085338 |
| fuel_surface_T_max_K | K | 417.90577 | 416.768675 | 415.62951 | 414.486546 | 413.344227 |
| gap_min_m | um | 113.00209 | 111.45695 | 109.912063 | 108.367691 | 106.8232 |
| gap_mean_m | um | 116.515279 | 115.259505 | 114.003927 | 112.748738 | 111.493447 |
| gap_max_m | um | 120.969623 | 120.070283 | 119.171057 | 118.272059 | 117.372961 |
| gap_min_position_local_z_m | cm | 18.9653333 | 18.9653333 | 18.9653333 | 18.9653333 | 18.9653333 |
| fuel_max_radial_displacement_m | um | 20.0285643 | 21.5737115 | 23.1186055 | 24.6629826 | 26.2074833 |
| clad_inner_ur_mean_m | um | 3.00924068 | 3.00924193 | 3.00924322 | 3.00924296 | 3.00924649 |
| fuel_von_mises_max_Pa | MPa | 30.8852236 | 28.464777 | 26.0543799 | 23.6946297 | 21.3936119 |
| clad_von_mises_max_Pa | MPa | 0.554370584 | 0.554403629 | 0.554436831 | 0.554454889 | 0.554500929 |
| swelling_eigenstrain_max | strain | 0 | 0.000111998554 | 0.000223997107 | 0.000335995661 | 0.000447994215 |
| Contatto | sì/no | no | no | no | no | no |

Potenza finale identica: 5265.814463925392 W. Nessun contatto a qualsiasi tempo per s≤2; nessuna estrapolazione del fattore necessario al contatto.

## Sensibilità locale e non linearità

Derivata centrale [Y(1.5)−Y(0.5)]/1; derivata seconda diagnostica [Y(1.5)−2Y(1)+Y(0.5)]/0.25. Nessuna assunzione di linearità della risposta.

| Quantità | Unità | ΔY/Δs al nominale | Seconda differenza /s² | Scarto massimo dalla tangente nominale | Range relativo al nominale [%] |
|---|---|---:|---:|---:|---:|
| gap_min_m | um | -3.08925948 | 0.00205470513 | 0.000766892174 | 5.62166673 |
| fuel_max_radial_displacement_m | um | 3.08927111 | -0.00206771459 | 0.000770102294 | 26.7270404 |
| Tmax_K | K | -2.20636937 | -0.0143036654 | 0.00554799584 | 0.911465947 |
| fuel_von_mises_max_Pa | MPa | -4.77014723 | 0.202587685 | 0.109379299 | 36.430004 |

Stress: il massimo su celle può cambiare posizione e non deve seguire una legge lineare. Gap/spostamento sono quasi lineari se le differenze seconde e gli scarti riportati sono piccoli; Tmax presenta il feedback del gap medio sulla conduttanza. Le percentuali di temperatura assoluta dipendono dalla scala: valutare soprattutto il range in K. Gap e spostamento sono quasi lineari (scarto massimo dalla tangente <0.0008 µm); Tmax <0.006 K e Von Mises <0.110 MPa. Piccole deviazioni sono riportate senza attribuirle tutte a non linearità fisiche risolte: restano tolleranze di convergenza iterative. Nessuna linearità imposta.

## Serie temporali e grafici

sensitivity_time_history.csv: 105 record = 5 fattori × 21 tempi, con sette serie richieste più potenza e FIMA. PNG/PDF comparativi per tutte le serie; sensitivity_3500h.png/pdf presenta gap, spostamento, Tmax, Von Mises vs s. Nessun punto oltre gli stati salvati.

## Interpretazione e limiti

Gli output meccanici del fuel sono i più sensibili (vedere range e percentuali misurati); l’effetto termico va giudicato sul range in K e quello sul clad sui valori molto piccoli. Questa classificazione è specifica del basso livello di FIMA e delle BC attuali. Nell’intervallo deterministico s=0–2, il range è 26.73% del nominale per lo spostamento e 36.43% per Von Mises fuel: influenza forte su questi output; gap minimo 5.62% (influenza moderata), Tmax 4.410 K (influenza contenuta), clad quasi invariato. Non è una stima probabilistica dell’incertezza né la previsione del range di una futura correlazione. Il modello gap/contact mantiene medie globali; nessuna nuova fisica o feedback Z3ST→OpenMC.

Warning nodali di stress sull’asse: 31 valutazioni non finite u_r/r sono azzerate dal writer come nei riferimenti; i massimi riportati sono valutati nei centri DG0 e finiti. Non affermiamo validità dei valori nodali sull’asse. I campi potenziali nominali del provider sono distinti dall’eigenstrain applicata nei nuovi VTU. Tutti i cinque solve accettati sono convergenti e con diagnostica PASS. Due tentativi iniziali s=0 erano convergenti ma la diagnostica aggiunta falliva: prima assenza di comunicatore MPI per un’espressione UFL identicamente nulla, poi assenza di dominio mesh. La sola diagnostica del campo nullo ora verifica esplicitamente lo zero simbolico e riempie il campo DG0 con zero; per s>0 valuta il wrapper costitutivo mediante Expression sul mesh. Output/log/preflight dei due tentativi conservati in runs/s_0.0/failed_diagnostic_attempt_{1,2}; nuovo preflight prima di ciascun retry. Nessuna modifica alla fisica.

## Risorse e protezione

1 MPI, 2 OMP, BLAS/NumExpr=1. Preflight free -h, df -h /, df -h /mnt/c, nproc, vmstat 1 3 prima di ogni solve, gate PASS. Risorse complete in runs/s_*/resource_preflight.json e solve_status.json.

| s | RAM available GiB | Swap usata GiB | Disco libero GiB | CPU logiche | Runtime s |
|---|---:|---:|---:|---:|---:|
| 0.0 | 8.289 | 0.692 | 922.996 | 8 | 5.088 |
| 0.5 | 8.250 | 0.692 | 922.989 | 8 | 6.777 |
| 1.0 | 8.256 | 0.691 | 922.982 | 8 | 7.436 |
| 1.5 | 8.247 | 0.691 | 922.974 | 8 | 7.443 |
| 2.0 | 8.245 | 0.691 | 922.967 | 8 | 7.892 |

Manifest prima/dopo: 1722 file protetti invariati. Nessun OpenMC/depletion, commit o push. Nessun output precedente sovrascritto.

## File aggiunti

README.md, sensitivity_law.py, run_sensitivity.py, postprocess.py, cinque configurazioni runs/s_*/ e relativi output, diagnostica locale, manifest prima/dopo, solve_status.json, sensitivity_results.json, sensitivity_3500h.csv, sensitivity_time_history.csv, local_sensitivity.csv, grafici PNG/PDF, questo report. Nessun file preesistente modificato.

Fase conclusa; attendere revisione.

Verifica aggiuntiva: eigenstrain applicata esattamente uguale bit per bit a s*FIMA in tutti i 105 stati (errore assoluto 0 per ogni s). Dettagli: eigenstrain_scaling_checks.json.
