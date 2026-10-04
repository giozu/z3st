# FE101 — FASE 2: heating nativo + FIMA nativa

Report completo: `PHASE2_REPORT.md`; dati, confronti e tutti i preflight:
`phase2_results.json`; conservazione FE ai 21 tempi:
`power_conservation_FE.csv`; mappe: `phase2_native_maps.png/pdf` e
`phase2_gap_deformation.png/pdf`. I campi temporali sono in `runs/*/output`.

La directory FASE 1 è protetta e non modificata. Questo è il primo coupling
**one-way** dei due campi disponibili, non un feedback bidirezionale verso
OpenMC e non una nuova validazione delle leggi materiali.

## Sorgente e unità

Si legge esclusivamente `deposited_power_W` dalla tabella validata nativa 10×10
high-stat. Il reader `HeatingHistory` estende `DepletionHistory` senza cambiare
il comportamento FIMA-only. I 21 statepoint sono stati confrontati in lettura
diretta HDF5 (`h5py`, nessun import o esecuzione OpenMC):

- tally materiale: `FE101 B1 domain energy deposition`;
- score: `heating-local`, eV per particella sorgente;
- normalizzazione: `factor-for-normalization`, heating-local integrato nel
  reattore, potenza 250000 W;
- potenza già salvata: `P_i=250000*H_i/H_core` W. La conversione eV→J si
  cancella nel rapporto; equivalentemente S=250000/(H_core*J_per_eV),
  P_i=H_i*J_per_eV*S;
- densità volumetrica: `q_i=P_i/(volume_cm3*1e-6)` W/m³.

Non si ricostruisce potenza da fissioni, FIMA o BU. Nessuna rinormalizzazione
alla baseline 8780 W/m. Gli stati salvati forniscono 5387.615085 W iniziali e
5265.814464 W finali al fuel, non 250 kW al singolo elemento.

## Trasferimento e ownership

Coordinate invariate dalla FASE 1: r_FE=r_cm/100,
z_FE=(z_cm−10.20)/100. Mesh copiata in questo caso isolato, lunghezza 0.3556 m,
nessuna estrapolazione o stretching. Con overlap cilindrico V_kj:

```
P_FE,k = sum_j [V_kj/V_source,j * P_source,j]
q_FE,k = P_FE,k / V_FE,k
```

`qdot_native` è DG0, nullo nel cladding; `q_third` lo referenzia per la weak
form esistente. `q_third_baseline` conserva separatamente il vecchio coefficiente,
che viene bypassato quando il canale è attivo. `fima_native` mantiene la sua
ownership indipendente. Non sono cambiate weak form, strain, stiffness,
swelling o gap/contact. L'unica modifica a `thermal_model.py` riguarda lo
spazio usato per la **stampa diagnostica** dei valori q_third.

Il BU interno resta CG1 e viene accumulato dal codice esistente con regola
right-endpoint. Per evitare indexing DG0/CG1 incompatibile, l'adapter fornisce
solo a questo diagnostico una proiezione CG1 mass-lumped positiva e conservativa:
q_node,a = integral(2πr N_a q_DG0)/integral(2πr N_a). La weak form termica
continua a usare DG0. Lo swelling nativo non usa il BU interno. Non si importa
BU OpenMC e non si pretende che l'integrazione right-endpoint coincida con
l'energia BOS della depletion.

Configurazione opt-in:

```yaml
coupling:
  native_power:
    enabled: true
    material: fuel
    axial_origin_cm: 10.20
    history_path: <CSV salvato con deposited_power_W>
```

È indipendente da `native_fima`, quindi si può verificare il solo termico.
Il supporto è limitato a quad lineari rettangolari r-z, con un solo materiale
sorgente selezionato. SCIANTIX/porosity non sono abilitati o introdotti.

## Tempo, test e risorse

Potenza **istantanea** interpolata linearmente fra i 21 snapshot; FIMA
cumulativa valutata allo stesso tempo fisico. Preservata positività e potenza
totale interpolata, nessuna extrapolazione oltre 3500 h. Snapshot/rollback
sono verificati per entrambi i campi. Questa interpolazione termica è una
scelta esplicita distinta dall'integrazione BOS dell'energia di depletion.

Cinque test dati PASS (provenienza HDF5, conservazione, tempi, nullo/uniforme,
input non fisici rifiutati); binding FE PASS (ownership, baseline ignorata,
proiezione conservativa, tempi simultanei, rollback, export DG0 XDMF).
Controlli prima del caso completo:

- A: conservazione potenza;
- B: import nullo, nessun uso accidentale della baseline;
- C: 5000 W uniformi contro sorgente uniforme Z3ST originale;
- D: storia q nativa senza swelling, tutti i 21 tempi;
- replay FASE 1 in copia isolata, confrontato con gli output originali.

Poi: q uniforme stessa potenza + FI nativa; q nativa + FI uniforme stessa
media; q nativa + FI nativa. I controlli hanno uguale potenza/FIMA media
a **ogni snapshot**, non soltanto a fine storia.

Prima di ogni solve sono eseguiti `free -h`, `df -h /`, `df -h /mnt/c`,
`nproc`, `vmstat 1 3`; output completi nei `resource_preflight.json`.
Gate: available RAM e disco ≥2 GiB, swap used ≤max(1 GiB,50% totale), nessuno
scambio attivo vmstat si/so >1024 KiB/s negli intervalli campionati. OMP=2,
MPI=1, BLAS/NumExpr=1, UCX_TLS=self; cache nel nuovo caso.

Riproducibilità, da ambiente z3st con PYTHONPATH alla root e
PYTHONDONTWRITEBYTECODE=1:

1. Su una copia nuova: `python prepare.py` (rifiuta output esistenti).
2. `python run_verification.py` esegue gate/test/controlli/completo/report.
3. `python run_verification.py --resume` verifica dati e binding, riutilizza
   soltanto solve già PASS e completa l'analisi senza sovrascriverli.

Non cancellare o sovrascrivere i risultati esistenti per ripetere la pipeline.
Manifest prima/dopo proteggono neutronica, FASE 1, regression e sorgenti non
coinvolti. `IMPLEMENTATION_FILES.json` elenca le modifiche esatte.

## Interpretazione e limiti

Il confronto uniforme nuova potenza vs FASE 1 isola il cambio di rating;
nativo vs uniforme alla stessa potenza isola la distribuzione termica;
nativo FI vs uniforme FI, con q identica, isola lo swelling eterogeneo.
Nessuna differenza viene attribuita automaticamente a feedback neutronico.

Gap/contact conservano **medie globali**. Il gap locale e la sua posizione
sono diagnostici dagli spostamenti, non un nuovo modello di contatto locale.
Restano outer clad 300 K, estremi isolati e correlazioni materiali esistenti.
Le tensioni riportate sono campioni al centro delle celle; i nodi r=0 del
writer legacy non sono usati per gli estremi. Le mappe T visualizzano medie
nodali di cella; Tmax è estratto direttamente dai nodi FE.

**FASE 2 verificata entro questi limiti. Nessuna nuova fisica, commit, push
o fase successiva automatica. Stop per revisione.**
