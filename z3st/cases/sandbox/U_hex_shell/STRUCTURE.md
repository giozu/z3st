<!-- # --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Bianca Funaro
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- -->

# U_hex_shell — folder map

Hexagonal 15-15Ti duct (wrapper). Every runnable subfolder is a case: `cd <subfolder> && ./Allrun`.
Skip `output/`, `__pycache__/`, `*.msh`, `*.vtu`, `*.h5` and `*.xdmf` when you search.

## Top level

| File                    | Role                                                                    |
|-------------------------|-------------------------------------------------------------------------|
| `hex_post.py`           | common post-processing: `python hex_post.py <subfolder>`, reads `<subfolder>/post.yaml` |
| `Allpost`               | runs `hex_post.py` on every `post.yaml` at any depth, then `VP1/compare_ansys.py` (default case `z_uniform/`) and `results_organize.py` |
| `results_organize.py`   | groups `results/` figures by problem (run after `./Allpost`)            |
| `k_1515.py`, `15_15Ti_local.yaml` | case-local 15-15Ti k(T) and material card                     |
| `results/`              | one folder per problem plus a `README.md` index (see below)             |

## Cases

| Subfolder                                         | Mesh                  | Allrun after the solve                                      | Gold |
|---------------------------------------------------|-----------------------|-------------------------------------------------------------|------|
| `2d/`                                             | `2d/mesh.geo`         | no Allrun (`input.yaml` paths are relative to `U_hex_shell/`) | yes |
| `inner_q_outer_T/`                                | gmsh -3 (via `utils/allrun.sh`, `DIM=3`) | `plots.py`; also `non-regression.py`, `case_params.py` | yes |
| `bowing/flat-flat/T_asymmetries/`                 | gmsh -3               | `plots.py`                                                  | no   |
| `bowing/corner-corner/IAEA/VP1/`                  | gmsh -3               | `non-regression.py`                                         | yes  |
| `bowing/corner-corner/IAEA/VP1/z_uniform/`        | `../mesh.msh`         | `hex_post.py .`, `../compare_ansys.py --case .`             | no   |
| `bowing/corner-corner/IAEA/VP1/1515Ti_Tdependent/`| `../mesh.msh`         | `non-regression.py`; also `compare_const.py`                | yes  |
| `bowing/corner-corner/IAEA/VP1/cylinder/`         | gmsh -3               | `vp1_cylinder_compare.py`                                   | no   |
| `bowing/corner-corner/IAEA/VP1/cylinder/fine/`    | own `mesh.geo`        | no Allrun; `fine_check.py`                                  | no   |
| `bowing/corner-corner/IAEA/VP1/cylinder/fine_hex/`| own `mesh.geo`        | no Allrun                                                   | no   |
| `bowing/corner-corner/IAEA/VP2/`                  | gmsh -3               | `vp2_post.py`, `vp2_benchmark.py`, `vp2_corners.py`, `../along_y.py .`, `../iwgfr_table.py VP2` | no |

Each case has its own `input.yaml`, `geometry.yaml`, `boundary_conditions.yaml` and `post.yaml` where present. `z_uniform/` and `1515Ti_Tdependent/` reuse the mesh, geometry and BCs of `VP1/`.

## IAEA scripts

- `VP1/`: `vp1_temperature.py` (imposed T; `Z_UNIFORM` switches the axial ramp off), `vp1_benchmark.py`, `vp1_corners.py`, `compare_ansys.py`, `gui_fields.py`, `VP1_notes.txt`, `15_15Ti_const.yaml` / `15_15Ti_local.yaml`
- `VP2/`: `vp2_temperature.py`, `vp2_post.py`, `vp2_benchmark.py`, `vp2_corners.py`, `VP2_notes.txt`
- `IAEA/`: shared `along_y.py`, `iwgfr.py`, `iwgfr_table.py`

## results/

`inner_q_flux/`, `bowing_flat-flat_T_asymmetries/`, `bowing_corner-corner_IAEA_VP1/{z-gradient,z-uniform,1515Ti_Tdep_vs_Ansys,1515Ti_Tdep_vs_const}/`, `bowing_corner-corner_IAEA_VP2/z-gradient/`. The z-gradient and z-uniform options are kept apart.
