# Single-edge notched shear test (SENS, plane strain)

Test case based on the SENS-shear benchmark of Ambati et al. 2015
(Comput Mech 55:383-405, §4.2), using the AT2 phase-field formulation with
the hybrid constraint (their Eq. 27).

Excluded from the local suite (`suite_exclude.txt`: 796 steps, about 2 h). Run it by hand.

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
- Star-convex energy split: `split: star_convex`, `gamma_star: 5.0`.
- Hybrid constraint (Ambati Eq. 27): `hybrid_constraint: true`.

## Solver (`input.yaml`)

- Staggered scheme, `max_iters: 200`, `relax_u: 1.0`, `relax_D: 0.8`,
  adaptive relaxation off (`relax_adaptive: false`).
- Mechanical: linear, `direct_mumps`.
- Damage: linear, `direct_mumps`.

## Expected results

Ambati Fig. 12 (hybrid row, p.398): a curved Mode-II crack starts at the notch tip
`(0.5 mm, 0.5 mm)`, arcs down-right and arrests in the lower-right region without
reaching the corner. The energy plots draw the Ambati Fig. 12d arrest target
`Gc · (Dn + 0.55 mm) ≈ 2.84 J`.

Stored run (`energies.txt`, `force_displacement.txt`, gold):
- `E_frac` starts at 1.358 J (regularised notch, `Gc · Dn ≈ 1.35 J`) and reaches
  4.669 J at the last step.
- `E_el` at the last step: 8.079 J.
- Peak reaction force 0.556 kN at `u_x = 27.4 µm` (step 782).
- `D_max = 1.0`.

<!-- [TBC] the stored final E_frac (4.67 J) exceeds the Ambati arrest target (2.84 J): crack arrest is not reproduced at 30 µm, or the target is not comparable -->

## Files

- `mesh.geo`, `geometry.yaml`: geometry and label map.
- `input.yaml`: physics, regime, solver options.
- `boundary_conditions.yaml`: clamped bottom, sliding top, `D = 1` on the notch.
- `bc_generator.py`: writes the displacement list in `boundary_conditions.yaml`.
- `diagnostics.py`: per-step hook, writes `force_displacement.txt` at the case root.
- `non-regression.py`: field plots of the last step, energy balance, gold checks.
- `plot_energy_balance.py`, `plot_force_displacement.py`: standalone plots.
- `sweep_gamma.sh`: runs the case for `gamma_star` = 0, 1, 5 into `output_starconvex_<tag>/`.
- `Allrun`, `Allclean`.

## Outputs

Written by the solver at the case root: `energies.txt`, `force_displacement.txt`.

Written by `non-regression.py` in `output/`:
- `damage_field.png`: D at the last step.
- `stress_vm_field.png`: von Mises stress at the last step.
- `stress_xy_field.png`: σ_xy at the last step.
- `stress_xx_field.png`: σ_xx at the last step.
- `crack_driving_force_field.png`: H at the last step (written only when the VTU carries H).
- `energy_balance.png`: E_el, E_frac, E_tot against step.
- `non-regression.json`: tracked values `D_max`, `sigma_vm_max_MPa`, `E_el_final`,
  `E_frac_final`, `E_tot_final`.

Written by `plot_force_displacement.py`: `output/force_displacement.png`.

## Running

```bash
./Allclean      # remove output and logs (mesh preserved)
./Allrun        # mesh -> z3st -> non-regression
```
