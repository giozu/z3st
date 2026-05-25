# Case 19 — Single-edge notched tension test (SENT, plane strain)

Reproduces the **Ambati et al. 2015** (Comput Mech 55:383–405) §4.1 SENT-tension
benchmark using the AT2 hybrid phase-field formulation (their Eq. 27). Companion
to `19_single-edge_notched_shear_test/` (same plate and notch, different BC).

## Geometry & loading

- 1 mm × 1 mm square plate, 2D plane strain.
- Horizontal notch of length 0.5 mm from the left edge to the centre, at
  `y = Ly/2`. Modelled as a zero-width slit, with `D = 1` Dirichlet along the
  notch to anchor the crack from step 0.
- Bottom edge fully clamped: `u = (0, 0)`.
- Top edge displaced vertically: `u = (0, u_y)`. Left and right edges are free.
- `u_y` ramps from 0 to 7 µm in steps of 0.01 µm (`bc_generator.ipynb`
  defaults). Ambati's hybrid peak (Fig. 9) lands at `u_y ≈ 5.6 µm`; the
  brittle drop is essentially vertical at that displacement.

## Material (`../../materials/high_carbon_steel.yaml`)

| Parameter | Value | Matches Ambati §4.1 |
|---|---|---|
| `E`  | 210 GPa | ✓ (λ = 121.15 GPa, μ = 80.77 GPa) |
| `nu` | 0.3 | ✓ |
| `Gc` | 2700 J/m² | ✓ (2.7×10⁻³ kN/mm) |

AT2 analytical threshold (derived in `spine.py` at load time):
`σ_c = sqrt(27·Gc·E / (256·ℓ_c)) ≈ 3.87 GPa` at `ℓ_c = 4 µm`.

## Mesh

- 2D triangulation, graded:
  - `h_fine = ℓ_c / 5 = 0.8 µm` at the notch tip and mouth (Points 5 and 6).
  - `h_coarse = Lx / 75 ≈ 13.3 µm` at the corners.
  - Gmsh interpolates linearly between the per-point sizes.
- The fine zone at the notch tip is required to resolve the Mode-I singular
  stress concentration that initiates the straight crack.
- The coarse mesh at the corners (h ≈ 2·ℓc) smears the BC singularity at
  the top/bottom-edge-meets-side-edge corners below the damage-growth rate.

## Phase-field formulation

- **AT2** crack-density functional.
- **Amor split** for `ψ⁺ / ψ⁻` decomposition (under tension, the volumetric
  and deviatoric parts both feed `ψ⁺` — the correct Mode-I driver).
- **Hybrid constraint** (Ambati Eq. 27): in cells where `ψ⁻ > ψ⁺`, the
  contribution of `ψ⁺` to `H` is zeroed.

## Solver

- Staggered scheme with adaptive relaxation off, `relax_u = 1.0`, `relax_D = 0.8`.
- Mechanical: linear, `direct_mumps` (robust against the `g(D) → K = 10⁻⁶`
  heterogeneity once cracks open).
- Damage: linear, `iterative_hypre` (AT2's mild `(H + 1)` mass coefficient).

## Expected results (Ambati Fig. 8 / Fig. 9, p.396)

- **Damage field**: straight horizontal Mode-I crack from the notch tip
  `(0.5 mm, 0.5 mm)` to the right edge `(1.0 mm, 0.5 mm)`. Ambati Fig. 8c at
  u = 6 µm shows full ligament traversal.
- **Force-displacement** (Fig. 9): linear up to `u_y ≈ 5.5 µm`, peak at
  `u_y ≈ 5.6 µm` and ~0.7 kN, near-vertical brittle drop as the crack
  traverses the ligament.
- **Energy balance**: `E_frac` starts at `Gc · Dn = 1.35 J` (regularised
  notch baseline) and ramps to `Gc · Lx = 2.7 J` once the crack reaches the
  right edge (full ligament traversal).

## Files

- `mesh.geo`, `geometry.yaml`         — geometry and label map.
- `input.yaml`                        — physics, regime, solver options.
- `boundary_conditions.yaml`          — clamped bottom, vertically displaced top,
                                         D=1 pre-crack.
- `bc_generator.ipynb`                — notebook to regenerate the displacement
                                         list. Sets `u0`, `u1`, `delta_u1`.
- `non-regression.py`                 — quick diagnostic on the last VTU.
- `Allrun`, `Allclean`                — case-14-style drivers.

## Diagnostic outputs (in `output/`)

- `damage_field.png`                  — 2D map of D at the final step
                                         (Ambati Fig. 8c reproducer); notch
                                         slit overlaid in cyan.
- `stress_yy_field.png`               — σ_yy (Mode-I tensile driver) at final.
- `stress_xx_field.png`               — σ_xx (transverse; Poisson response).
- `stress_vm_field.png`               — Von Mises at final.
- `crack_driving_force_field.png`     — H = (2·ℓc/Gc)·ψ⁺ at final.
- `energy_balance.png`                — E_el, E_frac, E_tot vs step with the
                                         notch-baseline (Gc·Dn = 1.35 J) and
                                         full-ligament (Gc·Lx = 2.7 J) refs.

## Running

```bash
# regenerate BC list if needed (open bc_generator.ipynb and run all cells,
# then update `n_steps` in input.yaml to match the printed count)
./Allclean
./Allrun
```
