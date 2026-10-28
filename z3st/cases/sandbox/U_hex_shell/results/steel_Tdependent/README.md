# steel_Tdependent

IAEA VP1 with axial gradient (400 °C below the core, ramp over 1.5-2.5 m, cold/hot corner
500/550 °C above; clamped bottom) with the T-dependent 15-15Ti card of `Materials.txt`
(Pizzocri 2013): E, eta, k, cp from the Sobolev or the Di Marcello set, T_ref = 400 °C.
Same mesh, BCs and T in every run; each run is compared with its own beam (hex_post VP1Tdep).

| prefix | runs | script (from `bowing/corner-corner/IAEA/VP1/`) |
|---|---|---|
| `IAEA_VP1_DiMarcello_vs_Sobolev_` | z3st Di Marcello vs z3st Sobolev | `python 1515Ti_Tdependent/compare_laws.py` |
| `IAEA_VP1_Sobolev_vs_Ansys_` | z3st vs Ansys, Sobolev | `python compare_ansys.py --case 1515Ti_Tdependent --ansys <A>/Sobolev` |
| `IAEA_VP1_DiMarcello_vs_Ansys_` | z3st vs Ansys, Di Marcello | `python compare_ansys.py --case 1515Ti_Tdependent/DiMarcello --ansys <A>/DiMarcello` |

`<A>` = `/mnt/c/Users/Mat/Desktop/Newcleo/FEM/ansys_1515Ti` (Mechanical, material from `1515Ti_Tdep.xml`,
SOLID185 enhanced strain, `ds.dat` and the export of `export_ansys_vp1.mac` per set).

Each prefix: `vtb_bowing`, `curvature`, `displacement_z`, `stress_z`, `along_s_{stress,strain_sz,normal}`,
`max_stress`, `differences` (table) and `summary.json`.

|  | z3st Sobolev | Ansys Sobolev | z3st Di Marcello | Ansys Di Marcello |
|---|---|---|---|---|
| bow at the top (mm) | 13.207 | 13.212 | 12.306 | 12.312 |
| beam (mm) | 13.227 | | 12.325 | |
| max von Mises (MPa) | 2.815 | 2.814 | 2.664 | 2.664 |

- Di Marcello / Sobolev: bow and curvature 0.932 (alpha from 400 °C, ~17.9 vs 19.2e-6 1/K);
  stresses 0.947 (E Di Marcello +1.6 %).
- Ansys - z3st at the top: +5.7 um (Sobolev), +5.3 um (Di Marcello), as with the constant card
  (+5.2 um): element formulation (Q1 vs ES), not the material.
- T rebuilt from the Ansys thermal strain = z3st T within 0.9 mK (Sobolev) and 0.3 mK (Di Marcello).
