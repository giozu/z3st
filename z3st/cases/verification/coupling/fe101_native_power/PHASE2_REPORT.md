# FE101 FASE 2 — FIMA nativa + heating nativo

**Coupling one-way dei due campi verificato; nessun nuovo calcolo OpenMC/depletion.**

Sorgente: `zest_FIMA_2D_10x10.csv`, colonna `deposited_power_W`. Verificati in sola lettura tutti i 21 HDF5: tally `FE101 B1 domain energy deposition`, score `heating-local` in eV/source; normalizzazione con `factor-for-normalization` (`heating-local` globale) a 250000 W di reattore. P_i=250000 H_i/H_core [W]; q_i=P_i/(volume_cm3*1e-6) [W/m3]. Non si usano FIMA o BU per ricostruire potenza.

Coordinate: r_FE=r_OpenMC/100; z_FE=(z_OpenMC−10.20)/100, lunghezza 0.3556 m. Overlap cilindrico: potenza estensiva trasferita e divisa per volume FE. Nessuno stretching, clamp o rinormalizzazione alla baseline.

Tempi: 21 snapshot 0–3500 h. Interpolazione lineare della potenza istantanea, non negativa; FIMA cumulativa interpolata alla stessa coordinata temporale. Nessuna extrapolazione. Questa scelta termica non è l’integrazione BOS dell’energia della depletion. BU interno resta distinto e non è usato dallo swelling nativo.

Conservazione FE: massimo errore assoluto 9.368e-11 W; massimo relativo 1.739e-14; tolleranza gate 1e-11 relativa (1e-10 W per campo nullo). Tutti i tempi e gli errori sono in power_conservation_FE.csv.

A PASS: dati/overlap/unità/tempi/off-grid/positività. B PASS: sorgente importata nulla, max|T−300 K|=3.132e-11 K anche con LHR baseline non nullo. C PASS: 5000 W uniformi via DG0 contro sorgente Z3ST originale, max ΔT=8.981e-12 K, tolleranza 1e-5 K. D PASS: thermal-only nativo completo, senza swelling. A–D PASS prima del caso completamente accoppiato.

Regressione FASE 1 in copia isolata: PASS; ΔT=5.627e-12 K, Δu=1.952e-18 m. Originali FASE 1 invariati.

| Metrica finale | Baseline + FI nativa (FASE1) | q uniforme stessa P + FI nativa | q nativa + FI uniforme | q nativa + FI nativa |
|---|---:|---:|---:|---:|
| Tmax fuel [K] | 397.633085051 | 454.922710309 | 483.862787179 | 483.862787179 |
| T superficie min [K] | 360.611780725 | 393.726553309 | 367.269984186 | 367.269984186 |
| T superficie media [K] | 360.611780725 | 393.726553309 | 393.950892851 | 393.950892851 |
| T superficie max [K] | 360.611780725 | 393.726553309 | 415.62951044 | 415.62951044 |
| gap min [um] | 118.815596932 | 113.178424869 | 110.541220367 | 109.912063453 |
| gap medio [um] | 119.449407185 | 113.812199701 | 114.003926921 | 114.003926921 |
| gap max [um] | 120.385038922 | 114.854889928 | 118.39932809 | 119.171056594 |
| ur fuel max [um] | 14.1558374457 | 19.8308201886 | 22.4894485748 | 23.1186054889 |
| ur clad interno medio [um] | 2.97143437141 | 3.00924504726 | 3.0092432175 | 3.0092432175 |
| von Mises fuel max [MPa] | 13.5661082629 | 24.3089493615 | 30.9158746104 | 26.0543798778 |
| von Mises clad max [MPa] | 0.26927836309 | 0.45410268343 | 0.554436831144 | 0.554436831144 |
| FIMA media | 0.000153060852823 | 0.000153060852823 | 0.000153060852823 | 0.000153060852823 |
| eigenstrain swelling max | 0.000223997107275 | 0.000223997107275 | 0.000153060852823 | 0.000223997107275 |

Caso completo: q min/media/max=9219738.87059/14694768.7625/21928562.2327 W/m3; P_fuel=5265.81446393 W; FIMA min/media/max=9.31671708113e-05/0.000153060852823/0.000223997107275.

Tmax in [r,z] locali m: [4.1425196606326175e-20, 0.18965333333333348]; q max nel box [ri,ro,zlo,zhi] m: [0.01716375, 0.01791, 0.14224, 0.1540933333333333]; FIMA max nel box [0.01716375, 0.01791, 0.1896533333333332, 0.2015066666666665]; gap minimo a z=0.18965333333333323 m.

Contatto medio: p=0.0 Pa; chiusura geometrica locale=False. ur clad interno min/medio/max=2.98140498/3.00924322/3.03076255 um.

Separazione degli effetti:

- Potenza totale: uniforme alla potenza OpenMC + FI nativa contro FASE 1, mantenendo FI e forma termica uniforme.
- Eterogeneità termica: q nativa + FI nativa contro q uniforme stessa potenza + FI nativa.
- Eterogeneità swelling: q nativa + FI nativa contro q nativa + FI uniforme alla stessa media a ogni tempo.
- Questi sono effetti deterministici della scelta dei campi nel modello FE; non dimostrano feedback neutronico o significatività statistica delle differenze.

Thermal-only nativo (D):

{
  "step": 20,
  "time_s": 12600000.0,
  "time_h": 3500.0,
  "Tmax_K": 486.0636085576469,
  "fuel_surface_T_min_K": 368.62480237127323,
  "fuel_surface_T_mean_K": 395.8095258564917,
  "fuel_surface_T_max_K": 417.90576956116865,
  "gap_min_m": 0.00011300208982627149,
  "gap_mean_m": 0.00011651527865437582,
  "gap_max_m": 0.00012096962290268287,
  "model_mean_gap_m": 0.0001165152192953732,
  "model_contact_pressure_Pa": 0.0,
  "local_gap_closed": false,
  "fuel_max_radial_displacement_m": 2.0028564273437284e-05,
  "clad_inner_ur_min_m": 2.9814322955994388e-06,
  "clad_inner_ur_mean_m": 3.0092406779887576e-06,
  "clad_inner_ur_max_m": 3.0307477336788337e-06,
  "internal_BU_mean_MWd_kgU": 4.272731049389796,
  "Tmax_position_local_m": [
    4.1425196606326175e-20,
    0.18965333333333348
  ],
  "gap_min_position_local_z_m": 0.18965333333333323,
  "qdot_min_W_m3": 9219738.87059204,
  "qdot_mean_W_m3": 14694768.762533339,
  "qdot_max_W_m3": 21928562.232676655,
  "qdot_max_cell_bounds_m": [
    0.01716375,
    0.01791,
    0.14224,
    0.1540933333333333
  ],
  "fuel_power_FE_W": 5265.814463925392,
  "fuel_power_source_W": 5265.8144639253405,
  "power_time_s": 12600000.0,
  "FIMA_time_s": null,
  "power_error_W": 5.184119800105691e-11,
  "power_relative_error": 9.844858446154308e-15,
  "fuel_von_mises_max_Pa": 30885223.56650345,
  "fuel_sigma_rr_min_Pa": -15396094.546902835,
  "fuel_sigma_rr_max_Pa": -194608.7453956604,
  "fuel_sigma_hoop_min_Pa": -15396094.546902835,
  "fuel_sigma_hoop_max_Pa": 30421856.84332156,
  "fuel_sigma_zz_min_Pa": -30963001.051112473,
  "fuel_sigma_zz_max_Pa": 30025570.0863809,
  "fuel_sigma_rz_min_Pa": -3410871.928636302,
  "fuel_sigma_rz_max_Pa": 556991.8952176322,
  "clad_von_mises_max_Pa": 554370.5837742044,
  "clad_sigma_rr_min_Pa": -6915.290781870484,
  "clad_sigma_rr_max_Pa": -1064.452659778297,
  "clad_sigma_hoop_min_Pa": -555514.7734697834,
  "clad_sigma_hoop_max_Pa": 545910.8652911112,
  "clad_sigma_zz_min_Pa": -558031.9816102609,
  "clad_sigma_zz_max_Pa": 543619.5718991533,
  "clad_sigma_rz_min_Pa": -5144.34308999991,
  "clad_sigma_rz_max_Pa": 777.1638152497428
}

Risorse prima di ogni solve (output free/df/nproc/vmstat nei resource_preflight.json):

| Caso | RAM available GiB | Swap used GiB | Disco GiB | CPU | OMP | Runtime s |
|---|---:|---:|---:|---:|---:|---:|
| zero_native | 8.2393 | 0.6933 | 923.0745 | 8 | 2 | 13.610 |
| uniform_native | 8.2431 | 0.6933 | 923.0681 | 8 | 2 | 4.303 |
| uniform_baseline | 8.2389 | 0.6932 | 923.0634 | 8 | 2 | 4.454 |
| native_thermal | 8.2886 | 0.6932 | 923.0586 | 8 | 2 | 7.298 |
| phase1_replay | 8.3013 | 0.6932 | 923.0524 | 8 | 2 | 22.478 |
| uniform_power_native_fima | 8.2986 | 0.6932 | 923.0458 | 8 | 2 | 7.081 |
| native_uniform_fima | 8.2682 | 0.6929 | 923.0401 | 8 | 2 | 6.777 |
| full_native | 8.2617 | 0.6929 | 923.0327 | 8 | 2 | 7.146 |

Protezione PASS: 966 file invariati: neutronica, intera FASE 1, regression 101/PWR e sorgenti non coinvolti.

Weak form termo-meccaniche e misure invariate. La sola modifica a thermal_model.py seleziona lo spazio del coefficiente per la stampa diagnostica q_third; il weak form è identico. Swelling e gap/contact non modificati.

Limiti:

- One-way transfer from fixed OpenMC states; no feedback to neutronics.
- Global mean gap/contact conductance and pressure unchanged; local gaps are diagnostics.
- Axis-aligned rectangular DG0 overlap averaging; bin jumps can be smoothed.
- Instantaneous power interpolated linearly; internal BU uses existing right-endpoint accumulation of a conservative positive CG1 projection, not saved BOS BU.
- Unchanged material correlations and outer 300 K/insulated-end boundary conditions; no coolant-energy model.
- Stresses sampled at DG0 cell centres; existing singular axis nodal exports are not used for extrema.

Mappe: phase2_native_maps.png/pdf; gap/deformazione: phase2_gap_deformation.png/pdf; tutti i campi ai 21 tempi nelle directory runs/*/output.

**FASE 2 e primo trasferimento FIMA + q verificati entro questi limiti. Nessuna nuova fisica o feedback bidirezionale. Stop per revisione; nessun commit/push.**
