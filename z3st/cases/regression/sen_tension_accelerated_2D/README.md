# Coarse SENT with lagged staggering, Anderson acceleration and notch history

Regression case for three damage options that no other suite case uses:
`damage.staggered: lagged`, `solver_settings.anderson_depth: 5` and
`damage.notch_history`. The geometry and material are those of
`benchmarks/damage/sen_tension/` on a coarse mesh, so the case runs in about
two minutes. It is not a comparison with Ambati et al. (2015).

## Setup

- 1 mm × 1 mm plate, 2D plane strain, `high_carbon_steel.yaml`, AT2, `lc = 4 µm`,
  no energy split.
- Notch from `(0, 0.5)` to `(0.5, 0.5)` mm, a line in the mesh with no `D = 1`
  condition. The crack is set by `notch_history` with the default
  `B = 1000` and `half_width = lc / 2`.
- Bottom edge clamped, top edge `u = (0, u_y)`: `u_y` from 0 to 5 µm in steps
  of 1 µm, then to 10 µm in steps of 0.1 µm, 56 steps.
- Mesh (`mesh.geo`): 5806 nodes, 2 µm at the notch tip, 12.5 µm at the notch
  mouth, 25 µm at the corners.
- Solver: `staggered: lagged`, `anderson_depth: 5`, `direct_mumps` for both
  fields, `stag_tol: 1e-5`.

## Stored run

- Every step converges in 2 or 3 staggered iterations.
- The crack grows from the notch tip along `y = 0.5 mm` by about 0.1 mm.
- Last step, per unit depth: `E_el` = 6.78 J/m, `E_frac` = 3.28 J/m.

## Files

- `mesh.geo`, `geometry.yaml`, `input.yaml`, `boundary_conditions.yaml`.
- `non-regression.py`: field plots of the last step (damage, stresses, H),
  energy balance, gold check of `D_max`, `sigma_vm_max_MPa`, `E_el_final`,
  `E_frac_final`, `E_tot_final`.
- `Allrun`, `Allclean`.
