<!-- # --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- -->

# 1515Ti_Tdep

IAEA VP1 with the axial gradient (z-gradient option): 400 °C below the core, a ramp over 1.5-2.5 m, and the cold/hot corner at 500/550 °C above it. The bottom is clamped. The material is the T-dependent 15-15Ti card of `Materials.txt` (Pizzocri 2013): E, eta, k and cp from the Sobolev or the Di Marcello set, T_ref = 400 °C. Every run uses the same mesh, BCs and T.

| folder | runs | script (from `bowing/corner-corner/IAEA/VP1/`) |
|---|---|---|
| `vs_const/` | z3st Sobolev, z3st Di Marcello, z3st constant card (../) | `python 1515Ti_Tdependent/compare_const.py`, `python 1515Ti_Tdependent/compare_const_corners.py` |
| `Sobolev_vs_Ansys/` | z3st vs Ansys, Sobolev | `python compare_ansys.py --case 1515Ti_Tdependent --ansys <A>/Sobolev` |
| `DiMarcello_vs_Ansys/` | z3st vs Ansys, Di Marcello | `python compare_ansys.py --case 1515Ti_Tdependent/DiMarcello --ansys <A>/DiMarcello` |

`<A>` = `/mnt/c/Users/Mat/Desktop/Newcleo/FEM/ansys_1515Ti`. This holds the Mechanical runs: material from `1515Ti_Tdep.xml`, SOLID185 with enhanced strain, plus `ds.dat` and the export of `export_ansys_vp1.mac` for each set.

## vs_const: error of the constant card

Each figure shows the three runs together. The two T-dependent sets are never compared alone. The figures carry no closed form, benchmark or tolerance: eta(T) and E(T) make the problem non-linear in T, and the differences between the runs are the result.

Δ = constant − T-dependent (unit of the quantity). On each profile along z, and on each curve along s as drawn, max |Δ| is given for each set, and as a % of max |T-dependent| on that profile.

Figures: `vtb_bowing`, `curvature`, `displacement_z`, `stress_z` (runs and Δ along z), `along_s_{stress,strain_sz,normal}` (the band under each figure gives max |Δ| for each height and set), `corners_z`, `max_stress`, `errors` (table) and `summary.json`.

|  | Sobolev | Di Marcello | constant | const − Sobolev | const − Di Marcello |
|---|---|---|---|---|---|
| bow at the top (mm) | 13.21 | 12.31 | 12.23 | −0.98 (−7.4 %) | −0.075 (−0.6 %) |
| κ at z = 3 m (1/m) | 6.52e-3 | 6.07e-3 | 6.00e-3 | −8.0 % | −1.2 % |
| max σ_vm (MPa) | 2.815 | 2.664 | 2.731 | −3.0 % | +2.5 % |
| max over z of the section max σ_vm (% of max) | | | | 5.0 % | 9.8 % |

## <law>_vs_Ansys

The figures and table are those of `compare_ansys.py`. Each run is compared with its own beam (hex_post VP1Tdep).

|  | z3st Sobolev | Ansys Sobolev | z3st Di Marcello | Ansys Di Marcello |
|---|---|---|---|---|
| bow at the top (mm) | 13.207 | 13.212 | 12.306 | 12.312 |
| beam (mm) | 13.227 | | 12.325 | |
| max von Mises (MPa) | 2.815 | 2.814 | 2.664 | 2.664 |

- Ansys − z3st at the top: +5.7 µm (Sobolev) and +5.3 µm (Di Marcello). The constant card gives +5.2 µm. The gap comes from the element formulation (Q1 vs ES), not from the material.
