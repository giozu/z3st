# FE101 — FASE 1: FIMA nativa OpenMC → swelling

Caso isolato di **verifica dell'interfaccia FIMA/swelling**. Non è una simulazione
operativa completamente accoppiata. Nessun import di potenza, q''', BU OpenMC,
creep o nuova legge costitutiva. La potenza Z3ST resta 8780 W/m; il burnup
interno resta una variabile distinta in MWd/kgU.

Il report completo è `PHASE1_REPORT.md`; dati numerici e confronto uniforme sono
in `phase1_results.json`, stato della pipeline in `verification_status.json`,
mappe in `phase1_FIMA_swelling_maps.png/pdf` e campi temporali VTU in
`runs/native/output/`. `IMPLEMENTATION_FILES.json` elenca i file esatti.

## Coordinate, trasferimento e tempi

Sorgente: tabella nativa 10×10 high-stat, 100 domini ×21 tempi, conservata
immutabile nel caso neutronico. Un hash SHA256 identifica il file letto.

Il nuovo caso usa un fuel lungo **0.3556 m**, con gli stessi raggi e la stessa
topologia 24×30 quad del regression. Il vecchio caso lungo 0.356 m resta
immutato. La differenza di 0.04 cm è quindi risolta nella geometria di questo
nuovo caso, senza stretching del campo o assegnazione per indice.

```
r_FE [m] = r_OpenMC [cm] / 100
z_FE [m] = (z_OpenMC [cm] - 10.20) / 100
t_FE [s] = t_OpenMC [h] * 3600
```

Le sovrapposizioni sono volumi cilindrici: π(ro²−ri²)Δz. Per cella FE k e dominio
sorgente j si calcola a_kj=V_overlap/V_source. Si trasferiscono separatamente
gli atomi iniziali N0_k=Σa_kj N0_j e le fissioni C_k(t)=Σa_kj C_j(t), poi si
ricostruisce FIMA_k=C_k/N0_k. L'inventario iniziale uniforme rende tale media
equivalente alla media volumetrica. Il campo FE è DG0, nullo nel cladding.
Il reader non richiede OpenMC, pandas o l'ambiente neutronico.

Interpolazione lineare delle **fissioni cumulative** fra tempi salvati, con
denominatori iniziali invariati. Tutti i 21 tempi sono letti; l'orizzonte è
3500 h. Tempi fuori dall'intervallo o celle fuori dal volume attivo generano
errore: nessuna estrapolazione, clamp o estensione a 4290 h.

Le celle attraversate da confini radiali mediano i salti. Su questa mesh il
massimo finale è preservato, ma il rebinning del campo rimappato introduce
fino a 0.684434% nel profilo radiale e 0.794318% nei singoli domini. Questi
scarti di trasferimento sono distinti dalle incertezze Monte Carlo.

## Ownership e configurazione opzionale

```yaml
coupling:
  native_fima:
    enabled: true
    material: fuel
    axial_origin_cm: 10.20
    history_path: <percorso alla tabella salvata>
```

La card fuel deve selezionare esplicitamente `fima_source: native_openmc` e la
legge `materials.fuel_swelling.uzrh_fission_product_swelling`. In quel modo
usa solo `model.fima_native`; se il campo manca, fallisce. La modalità legacy
di default resta `legacy_burnup`. La relazione è condivisa e indipendente
dall'origine: ε_sw=FIMA I, ΔV/V=3 FIMA. Il campo esportato
`Swelling_eigenstrain_native` è il valore di ciascuna delle tre componenti
diagonali, con componenti fuori diagonale nulle; non è lo strain termico totale.

L'interfaccia aggiorna i coefficienti persistenti in-place, senza modificare
q_third o burnup. I passi adattivi ricevono il tempo fisico del sottopasso;
snapshot/restore preservano coefficienti e tempo del provider. Sono supportate
attualmente celle quad lineari rettangolari assialsimmetriche, non triangoli,
quad distorti o trasferimenti 3D generici.

## Riproducibilità e gate

Usare l'ambiente con FEniCSx (`conda activate z3st`) e `PYTHONPATH` alla radice
del repository. La pipeline si rifiuta di sovrascrivere output esistenti;
una ripetizione completa deve usare un checkout/copia separata e proteggere
i risultati precedenti. Non eseguire `prepare.py` sopra solve già prodotti.

Per una preparazione nuova:

1. `python prepare.py`: costruisce solamente questo caso e i controlli.
2. `python run_verification.py`: test dati, binding FE, controlli C/D; solo
   dopo PASS esegue uniforme alla stessa media e campo nativo; genera report.
3. `python run_verification.py --postprocess-only`: verifica i solve salvati,
   ripete solo test di binding senza solve e rigenera l'analisi.

Prima di **ogni** solve la pipeline registra RAM, swap, disco, CPU e gate nel
sottocaso. Richiede RAM ≥2 GiB, disco ≥2 GiB e swap usata ≤max(1 GiB,50% totale).
Un gate fallito impedisce l'avvio. OMP=2, MPI=1, BLAS/NumExpr=1,
UCX_TLS=self; cache matplotlib/FFCX in directory scrivibili del nuovo caso.

Controlli eseguiti:

- A/B: identità, conservazione, tutti i tempi, interpolazione fuori griglia,
  monotonicità, indici riordinati, input corrotti/range/volumi rifiutati,
  profili e picco; sette test PASS.
- Binding FE: ownership, isotropia, rollback, mancato fallback, clad nullo,
  potenza immutata e comportamento legacy; PASS senza solve.
- C: campo uniforme contro espansione libera analitica e contro il canale
  esistente di swelling volumetrico costante; PASS entro 1 nm.
- D: storia nulla contro stesso baseline con swelling disabilitato; PASS.
- E: storia vera 10×10 fino a 3500 h e controllo uniforme con uguale media
  **a ogni tempo**; entrambi completati.

L'analisi verifica gli hash prima/dopo di tutti i risultati neutronici e dei
due regression originali. Le weak form, misure, rigidezza e strain restano
identici all'audit. Le tensioni riportate sono campioni DG0 al centro delle
celle; il writer legacy azzera valori singolari u_r/r ai nodi dell'asse.
Le energie elastiche esportate dal baseline non sono state validate come
energie comprensive dello swelling in questa fase.

La pressione di contatto e la conductance restano medie globali; il gap locale
è un diagnostico ricostruito dagli spostamenti. Le temperature uguali fra
nativo/uniforme non dimostrano un coupling termico locale risolto.

I tentativi tecnici sono documentati: correzione del keyword `cells0` nel
diagnostico e ripresa del postprocessing senza ripetere i solve. Nessun
cambiamento della fisica per aggirare problemi di esecuzione.

**FASE 1 verificata entro questi limiti. Nessuna FASE 2 automatica.**
