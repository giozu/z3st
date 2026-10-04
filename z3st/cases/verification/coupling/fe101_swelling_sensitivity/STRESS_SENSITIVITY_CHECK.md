# Stress sensitivity diagnostic check

PASS — numerical / model-effect isolation. Non è validazione fisica della legge di swelling. Solo lettura degli output già disponibili; nessuna simulazione o ricostruzione FEM.

## Mappe e massimi

Fonte: `fields_0020.vtu`, campi `Stress (cells)` e `VonMises (cells)`. Le mappe stress_map_s0/s1/s2.png/pdf usano la stessa scala colori. Media pesata sul volume cilindrico π(r_out²−r_in²)Δz, non media aritmetica della griglia. Coordinate z locali FE; z OpenMC = z_FE +10.20 cm.

| s | r max [cm] | z max locale [cm] | z OpenMC [cm] | VM min [MPa] | VM media volumetrica [MPa] | VM max [MPa] |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 1.753688 | 19.558000 | 29.758000 | 1.943970 | 13.939570 | 30.885224 |
| 1 | 1.753688 | 18.372667 | 28.572667 | 1.674887 | 12.077573 | 26.054380 |
| 2 | 1.753688 | 16.002000 | 26.202000 | 1.406296 | 10.244363 | 21.393612 |

Il massimo migra assialmente di −3.556 cm tra s=0 e s=2 (tre celle FE), restando nella stessa fascia radiale esterna: r=1.7536875 cm, a 0.0373125 cm dalla superficie. Non è sull’asse né alle estremità: la distanza dall’estremità più vicina resta almeno 15.8 cm. Il gap minimo resta a z locale=18.965333 cm. Il massimo qdot è nella regione esterna z-bin 5 (14.224–17.780 cm), il massimo FIMA nella regione esterna z-bin 6 (17.780–21.336 cm); s=0/1 hanno massimo stress vicino al massimo FIMA, s=2 migra nella zona di massimo qdot. Non implica causalità esclusiva.

## Componenti al massimo di ciascun caso

Tensore assialsimmetrico nell’ordine (r, θ, z); valori in MPa.

| s | σrr | σθθ | σzz | σrz |
|---|---:|---:|---:|---:|
| 0 | -0.662461 | 30.416136 | 30.025570 | 0.034268 |
| 1 | -0.699730 | 25.533273 | 25.172274 | -0.003898 |
| 2 | -0.374205 | 20.999224 | 21.034366 | -0.192219 |

Per evitare di confondere redistribuzione con migrazione del massimo, stress_components_fixed_regions.csv confronta anche le stesse tre celle in tutti i casi.

| Regione fissa: massimo s=0 | s | VM | σrr | σθθ | σzz | σrz |
|---|---:|---:|---:|---:|---:|---:|
| r=1.753688, z=19.558000 cm | 0 | 30.885224 | -0.662461 | 30.416136 | 30.025570 | 0.034268 |
| r=1.753688, z=19.558000 cm | 1 | 26.001730 | -0.471573 | 25.656635 | 25.401025 | 0.081562 |
| r=1.753688, z=19.558000 cm | 2 | 21.119173 | -0.281604 | 20.896458 | 20.775806 | 0.128855 |

## Spiegazione meccanica e coerenza

La trazione circonferenziale e assiale nella fascia esterna diminuisce, mentre σrr resta piccola e compressiva; si riducono le differenze fra le tensioni principali, quindi la norma deviatorica. Il taglio è molto più piccolo delle tensioni normali. La riduzione avviene anche confrontando una regione fissa e non è un artefatto del cambio di cella massima.

La swelling eigenstrain è isotropa: a spostamento fissato aggiungerebbe un contributo idrostatico e non cambierebbe direttamente Von Mises. Nel solve il campo di spostamento si riadatta alla eigenstrain spazialmente variabile; cambia quindi la deformazione elastica deviatorica. Si tratta di redistribuzione elastica, non creep/plasticità o rilassamento temporale introdotto. Inoltre le temperature diminuiscono lievemente attraverso il gap medio: questa verifica non separa quantitativamente tale contributo termico dal contributo meccanico diretto.

| s | ur,max fuel [µm] | gap min [µm] | z gap min [cm] | VM clad max [MPa] | contatto |
|---|---:|---:|---:|---:|---|
| 0 | 20.028564 | 113.002090 | 18.965333 | 0.554371 | no |
| 1 | 23.118605 | 109.912063 | 18.965333 | 0.554437 | no |
| 2 | 26.207483 | 106.823200 | 18.965333 | 0.554501 | no |

qdot/FIMA identiche elemento per elemento; potenza fuel identica (5265.814463925392 W). Nessun contatto a qualsiasi dei 21 tempi.

## Warning e limiti

Tutte le 720 celle fuel hanno stress e VM finiti. VM ricalcolata dal tensore salvato concorda con il campo cellwise e con i report precedenti; massimo errore 7.45e-09 Pa. Il massimo è nell’ultima fascia radiale, lontano dall’asse. Non si usano né si mediano i valori nodali: il warning preesistente riguarda 31 nodi sull’asse (u_r/r); il writer li ha azzerati. Nessuna contaminazione nodale nelle mappe cellwise. La verifica è discreta: i massimi continui/subcella non sono determinati. Gap/contact resta mediato globalmente; correlazione U-ZrH provvisoria.

Hash prima/dopo: 330 file preesistenti invariati. Nessun solver/input/correlazione/sorgente modificato, nessun OpenMC/depletion, commit/push. Aggiunti soltanto script diagnostico, tre mappe PNG/PDF, CSV componenti, JSON risultati/manifests e questo report.

## Giudizio finale

**Comportamento coerente** con i campi salvati: la diminuzione del massimo Von Mises è accompagnata da una riduzione dello stato deviatorico nelle stesse regioni e da una moderata migrazione assiale, con espansione radiale crescente e gap decrescente. Non emergono anomalie che richiedano nuovi solve.

**swelling coupling + ablation + sensitivity block closed**
