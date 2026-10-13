# Single-edge notched shear test (SENS, plane strain)

Test case based on the SENS-shear benchmark of Ambati et al. 2015
(Comput Mech 55:383-405, §4.2), using the AT2 phase-field formulation with
the hybrid constraint (their Eq. 27).

Excluded from the local suite (`suite_exclude.txt`: 796 steps, about 1.3 h). Run it by hand.

## Geometry and loading

- 1 mm × 1 mm square plate, 2D plane strain.
- Horizontal notch of length 0.5 mm from the left edge to the centre, at
  `y = Ly/2`. Modelled as a zero-width slit, with a `D = 1` Dirichlet condition
  on the notch (`region: crack`).
- Bottom edge clamped: `u = (0, 0)`.
- Top edge: `u = (u_x, 0)`. Left and right edges are free.
- `u_x` ramps from 0 to 30 µm in 796 steps (`n_steps: 796`), in four segments
  (`bc_generator.py`):

| Segment | Δu |
|---|---|
| 0 to 5 µm | 1 µm |
| 5 to 8 µm | 0.2 µm |
| 8 to 15 µm | 0.01 µm |
| 15 to 30 µm | 0.2 µm |

## Material (`z3st/materials/high_carbon_steel.yaml`)

| Parameter | Value | Ambati §4.2 |
|---|---|---|
| `E`  | 210 GPa | same (λ = 121.15 GPa, μ = 80.77 GPa) |
| `nu` | 0.3 | same |
| `Gc` | 2700 J/m² | same (2.7×10⁻³ kN/mm) |

AT2 analytical threshold (derived in `spine.py` at load time):
`σ_c = sqrt(27·Gc·E / (256·ℓ_c)) ≈ 3.87 GPa` at `ℓ_c = 4 µm`.

## Mesh (`mesh.geo`)

- 2D triangulation, graded:
  - `h_fine = ℓ_c / 5 = 0.8 µm` at the notch tip and mouth (Points 5 and 6).
  - `h_coarse = Lx / 100 = 10 µm` at the corners.
  - Gmsh interpolates linearly between the per-point sizes.
- No `D = 0` condition is applied at the corners.

## Phase-field formulation (`input.yaml`)

- AT2 crack-density functional, `lc = 4 µm`.
- Spectral energy split (Miehe), as in Ambati: `split: miehe`.
- Hybrid constraint (Ambati Eq. 27): `hybrid_constraint: true`.

## Solver (`input.yaml`)

- Staggered scheme with the history field frozen within each step at the
  displacement of the previous step (Miehe 2010, Box 1): `staggered: lagged`.
  `max_iters: 200`, `relax_u: 1.0`, `relax_D: 0.8`,
  adaptive relaxation off (`relax_adaptive: false`).
- Anderson acceleration of the damage update, `anderson_depth: 5`
  (`anderson_switch` at its default 5).
- Mechanical: linear, `direct_mumps`.
- Damage: linear, `direct_mumps`.

## Expected results

Ambati Fig. 12 (hybrid row, p. 398): a curved Mode-II crack starts at the notch
tip `(0.5 mm, 0.5 mm)` and runs down-right to the bottom edge near the
lower-right corner. Fig. 13, hybrid curve: first peak about 0.49 kN at about
12.5 µm, minimum about 0.35 kN near 16 µm, about 0.47 kN at 30 µm.

Stored run (`energies.txt`, `force_displacement.txt`, gold), 796 steps,
1 h 16 min, every step converged in 2 or 3 iterations:
- First peak 0.531 kN at `u_x = 10.6 µm`.
- Minimum 0.442 kN at 15 µm, then the force rises to 0.538 kN at 20 µm and
  0.588 kN at 30 µm.
- The crack starts at the notch tip, kinks down-right and at 30 µm reaches
  about `(0.74, 0.10) mm`. The damage stays below 0.14 on the top and bottom
  edges.
- Last step, per unit depth: `E_el` = 8.74 J/m, `E_frac` = 4.30 J/m.
  `Gc` times the crack length (about 0.97 mm) is 2.6 J/m.
- The largest von Mises stress, 14.5 GPa, is at the clamped lower-right
  corner `(1, 0) mm`, a corner singularity, where the damage stays 0.

`output/damage_steps.png` shows the damage at steps 300, 550 and 795, and
`output/damage_evolution.gif` the crack growth.

## Differences from Ambati

Measured in `INVESTIGATION.md` (section 20):
- The curve has Ambati's shape (peak, minimum, rise) but lies above it: first
  peak 8 % higher and 1.9 µm earlier, minimum 26 % higher and 1 µm earlier,
  force at 30 µm 25 % higher.
- `staggered: lagged` delays the history by one increment. On Ambati's mesh
  this raised the peak by 2.4 % against the coupled scheme (section 18).
- The notch is a line with `D = 1` here and an initial history field in
  Ambati; on Ambati's mesh the notch representation sets the initial slope
  (sections 13 and 14).
- At 30 µm Ambati's crack reaches the bottom edge; here its tip is at
  `y = 0.10 mm`.
- With the coupled scheme (`staggered: coupled`) the spectral split almost
  never converges (section 10). With the Amor split it converges in most
  steps, but the bottom ligament breaks at 25 µm (section 19).

## Files

- `mesh.geo`, `geometry.yaml`: geometry and label map.
- `input.yaml`: physics, regime, solver options.
- `boundary_conditions.yaml`: clamped bottom, sliding top, `D = 1` on the notch.
- `bc_generator.py`: writes the displacement list in `boundary_conditions.yaml`.
- `diagnostics.py`: per-step hook, writes `force_displacement.txt` at the case root.
- `non-regression.py`: field plots of the last step, energy balance, gold checks.
- `plot_energy_balance.py`, `plot_force_displacement.py`: standalone plots.
- `plot_damage_evolution.py`: damage at three steps and the crack-growth animation;
  needs all the VTU files and `output/force_displacement.csv`.
- `INVESTIGATION.md`: comparison with Ambati, one factor at a time.
- `sweep_gamma.sh`: star-convex study, runs the case for `gamma_star` = 0, 1, 5 into
  `output_starconvex_<tag>/`. It needs `split: star_convex` and a `gamma_star`
  line in `input.yaml`, which the stored case does not have.
- `Allrun`, `Allclean`.

## Outputs

Written by the solver: `output/fields_NNNN.vtu` (one per step), `energies.txt`
at the case root. Written by `diagnostics.py` at the case root:
`force_displacement.txt`. `Allrun` writes the log `log_z3st.md`.

Written by `non-regression.py` in `output/`:
- `damage_field.png`: D at the last step.
- `stress_vm_field.png`: von Mises stress at the last step.
- `stress_xy_field.png`: σ_xy at the last step.
- `stress_xx_field.png`: σ_xx at the last step.
- `crack_driving_force_field.png`: H at the last step, averaged from cells to points.
- `energy_balance.png`: E_el, E_frac, E_tot against step.
- `non-regression.json`: tracked values `D_max`, `sigma_vm_max_MPa`, `E_el_final`,
  `E_frac_final`, `E_tot_final`.

Written by `plot_force_displacement.py`: `output/force_displacement.csv` and
`output/force_displacement.png`.

Written by `plot_damage_evolution.py`: `output/damage_steps.png` and
`output/damage_evolution.gif`.

## Running

```bash
./Allclean      # remove output and logs; keeps the mesh, the last VTU,
                # force_displacement.csv and the gold
./Allrun        # mesh -> z3st -> non-regression
```
