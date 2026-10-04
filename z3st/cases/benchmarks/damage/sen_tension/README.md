# Single-edge notched tension test (SENT, plane strain)

Test case based on the SENT-tension benchmark of Ambati et al. 2015
(Comput Mech 55:383-405, §4.1), using the AT2 phase-field formulation.
Companion to `benchmarks/damage/sen_shear/` (same plate and notch, different load).

Excluded from the local suite (`suite_exclude.txt`: 701 steps). Run it by hand.

## Geometry and loading

- 1 mm × 1 mm square plate, 2D plane strain.
- Horizontal notch of length 0.5 mm from the left edge to the centre, at
  `y = Ly/2`. Modelled as a zero-width slit, with a `D = 1` Dirichlet condition
  on the notch (`region: crack`).
- Bottom edge clamped: `u = (0, 0)`.
- Top edge displaced vertically: `u = (0, u_y)`. Left and right edges are free.
- `u_y` ramps from 0 to 7 µm in steps of 0.01 µm, 701 steps (`n_steps: 701`).
  The list is written by `bc_generator.ipynb`.

## Material (`z3st/materials/high_carbon_steel.yaml`)

| Parameter | Value | Ambati §4.1 |
|---|---|---|
| `E`  | 210 GPa | same (λ = 121.15 GPa, μ = 80.77 GPa) |
| `nu` | 0.3 | same |
| `Gc` | 2700 J/m² | same (2.7×10⁻³ kN/mm) |

AT2 analytical threshold (derived in `spine.py` at load time):
`σ_c = sqrt(27·Gc·E / (256·ℓ_c)) ≈ 3.87 GPa` at `ℓ_c = 4 µm`.

## Mesh (`mesh.geo`)

- 2D triangulation, graded:
  - `h_fine = ℓ_c / 7 ≈ 0.57 µm` at the notch tip and mouth (Points 5 and 6).
  - `h_coarse = Lx / 75 ≈ 13.3 µm` at the corners.
  - Gmsh interpolates linearly between the per-point sizes.

## Phase-field formulation (`input.yaml`)

- AT2 crack-density functional, `lc = 4 µm`.
- No `split` key, so AT2 uses the default Miehe spectral split
  (`damage_model.py`).
- No `hybrid_constraint` key.

<!-- [TBC] the case is described as the Ambati hybrid formulation; input.yaml sets neither hybrid_constraint nor a split, check which is intended -->

## Solver (`input.yaml`)

- Staggered scheme, `max_iters: 500`, `relax_u: 1.0`, `relax_D: 0.8`,
  adaptive relaxation off.
- Mechanical: linear, `direct_mumps`, `stag_tol: 1e-5`.
- Damage: linear, `iterative_hypre`, `stag_tol: 1e-5`.
- No time adaptivity: a step that reaches `max_iters` is accepted with its last iterate.

## Convergence of the stored run (`log_z3st.md`)

- 698 of 701 steps converge. Steps 539, 573 and 650 reach the 500-iteration cap
  and are accepted. They coincide with the largest jumps of fracture energy.
- 20 650 staggered iterations in total: median 7 per step, 90th percentile 80,
  maximum 399. Thirty steps take more than 200 iterations, all in the crack-growth phase.
- Wall time 9 h 09 min.
- In the three capped steps the residuals fall for tens of iterations, then rise
  to a plateau. The step after (651) converges in 167 iterations with
  falling residuals.

Missing capabilities related to this behaviour are listed in `z3st/ai/CONTEXT.md` §10.

## Expected results

Ambati Fig. 8 and Fig. 9 (p.396): a straight horizontal Mode-I crack from the notch
tip `(0.5 mm, 0.5 mm)` to the right edge `(1.0 mm, 0.5 mm)`. Force-displacement
linear up to `u_y ≈ 5.5 µm`, peak at `u_y ≈ 5.6 µm` and about 0.7 kN, then a
near-vertical drop.

The energy plot draws two references: notch baseline `Gc · Dn = 1.35 J` and full
ligament `Gc · Lx = 2.7 J`.

Stored run (`energies.txt`, gold):
- `E_frac`: 1.357 J at step 0, 2.521 J at step 700.
- `E_el` at step 700: 1.943 J.
- `D_max = 1.0`, `sigma_vm_max_MPa = 3212`.

## Files

- `mesh.geo`, `geometry.yaml`: geometry and label map.
- `input.yaml`: physics, regime, solver options.
- `boundary_conditions.yaml`: clamped bottom, vertically displaced top, `D = 1` on the notch.
- `bc_generator.ipynb`: notebook that writes the displacement list. Sets `u0`, `u1`, `delta_u1`.
- `non-regression.py`: field plots of the last VTU, energy balance, gold checks.
- `plot_force_displacement.py`: force-displacement plot.
- `Allrun`, `Allclean`.

## Outputs

Written by `non-regression.py` in `output/`:
- `damage_field.png`: D at the last step.
- `stress_yy_field.png`: σ_yy at the last step.
- `stress_xx_field.png`: σ_xx at the last step.
- `stress_vm_field.png`: von Mises stress at the last step.
- `crack_driving_force_field.png`: H at the last step (written only when the VTU carries H).
- `energy_balance.png`: E_el, E_frac, E_tot against step with the two references above.
- `non-regression.json`: tracked values `D_max`, `sigma_vm_max_MPa`, `E_el_final`,
  `E_frac_final`, `E_tot_final`.

Written by `plot_force_displacement.py`: `output/force_displacement.png`.

## Running

```bash
# to regenerate the BC list: run all cells of bc_generator.ipynb,
# then set n_steps in input.yaml to the printed count
./Allclean
./Allrun
```
