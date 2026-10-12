# Single-edge notched shear test (SENS, plane strain)

Test case based on the SENS-shear benchmark of Ambati et al. 2015
(Comput Mech 55:383-405, §4.2), using the AT2 phase-field formulation with
the hybrid constraint (their Eq. 27).

Excluded from the local suite (`suite_exclude.txt`: 796 steps, about 8 h). Run it by hand.

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
- Amor energy split: `split: amor`.
- Hybrid constraint (Ambati Eq. 27): `hybrid_constraint: true`.

## Solver (`input.yaml`)

- Staggered scheme, `max_iters: 200`, `relax_u: 1.0`, `relax_D: 0.8`,
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
8 h 7 min, 41 steps not converged within 200 iterations:
- First peak 0.460 kN at `u_x = 10.5 µm`.
- After the peak the force falls to 0.42 kN at 12 µm, then stays between 0.36
  and 0.41 kN (0.393 kN at 15 µm, 0.402 kN at 20 µm, 0.363 kN at 24.8 µm).
- At 25.0 µm (step 770) the ligament between the crack and the lower-right
  corner breaks along the bottom edge and the force drops to 0.003 kN. The
  main crack has then reached about `(0.71, 0.10) mm`.
- Last step, per unit depth: `E_el` = 0.053 J/m, `E_frac` = 6.07 J/m.

## Differences from Ambati

Measured in `INVESTIGATION.md`:
- Peak 6 % below Ambati's and 2 µm earlier. With the star-convex split
  (`gamma_star: 5`) the peak is 0.455 kN at 10.2 µm, so the split moves it by
  1 %. The notch is a line with `D = 1` here and an initial history field in
  Ambati; on Ambati's mesh the notch representation sets the initial slope
  (sections 13 and 14).
- No post-peak minimum and rise, and the bottom ligament breaks. Ambati drives
  the damage with the spectral split. Amor puts the whole deviatoric energy
  in `psi+`, also under compression: at the lower-right corner, where the
  strain trace is about -2 %, Amor's `psi+` exceeds `psi-` and the damage
  grows, while the spectral `psi+` stays 7 to 60 times below `psi-`. The
  spectral split is not used because the staggered loop almost never converges
  with it (section 10).
- 41 steps do not converge: Gauss points switch side of the hybrid
  constraint (`psi- > psi+`) between iterations and the damage increment
  stalls at about 1e-6, above `stag_tol: 1e-7` (sections 6 and 19). The last
  iterate is kept.

## Files

- `mesh.geo`, `geometry.yaml`: geometry and label map.
- `input.yaml`: physics, regime, solver options.
- `boundary_conditions.yaml`: clamped bottom, sliding top, `D = 1` on the notch.
- `bc_generator.py`: writes the displacement list in `boundary_conditions.yaml`.
- `diagnostics.py`: per-step hook, writes `force_displacement.txt` at the case root.
- `non-regression.py`: field plots of the last step, energy balance, gold checks.
- `plot_energy_balance.py`, `plot_force_displacement.py`: standalone plots.
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

## Running

```bash
./Allclean      # remove output and logs; keeps the mesh, the last VTU,
                # force_displacement.csv and the gold
./Allrun        # mesh -> z3st -> non-regression
```
