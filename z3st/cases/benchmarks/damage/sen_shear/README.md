# Case 19 — Single-edge notched shear test (SENS, plane strain)

Reproduces the **Ambati et al. 2015** (Comput Mech 55:383–405) §4.2 SENS-shear
benchmark using the AT2 hybrid phase-field formulation (their Eq. 27).

## Geometry & loading

- 1 mm × 1 mm square plate, 2D plane strain.
- Horizontal notch of length 0.5 mm from the left edge to the centre, at
  `y = Ly/2`. Modelled as a zero-width slit, with `D = 1` Dirichlet along the
  notch to anchor the crack from step 0.
- Bottom edge fully clamped: `u = (0, 0)`.
- Top edge slides horizontally: `u = (u_x, 0)` (the `u_y = 0` clamp prevents
  vertical lifting). Left and right edges are free.
- `u_x` ramps from 5 µm (step 0) to 12 µm (step 700) in 0.01 µm increments.
  Ambati's hybrid peak (Fig. 13) lands at `u_x ≈ 13–14 µm`; the present
  window catches the initiation and early propagation phase.

## Material (`../../../materials/high_carbon_steel.yaml`)

| Parameter | Value | Matches Ambati §4.2 |
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
- The fine zone at the notch tip is required to resolve the Mode-II singular
  stress concentration that initiates the curved crack.
- The coarse mesh at the corners (h ≈ 2·ℓc) smears the BC singularity at
  the top/bottom-edge-meets-side-edge corners below the damage-growth rate,
  avoiding the need for explicit `D = 0` corner BCs.

## Phase-field formulation

- **AT2** crack-density functional.
- **Amor split** for `ψ⁺ / ψ⁻` decomposition (under shear, the deviatoric
  component goes into `ψ⁺` — the correct Mode-II driver).
- **Hybrid constraint** (Ambati Eq. 27): in cells where `ψ⁻ > ψ⁺`, the
  contribution of `ψ⁺` to `H` is zeroed, suppressing damage in compression.

## Solver

- Staggered scheme with adaptive relaxation (defaults).
- Mechanical: linear, `iterative_hypre` (CG + BoomerAMG).
- Damage: linear, `iterative_hypre` (CG + BoomerAMG on AT2's mild `(H + 1)`
  mass coefficient).

## Expected results (Ambati Fig. 12 hybrid row, p.398)

- **Damage field**: curved Mode-II crack initiating at the notch tip
  `(0.5 mm, 0.5 mm)`, arcing ~45° down-right, **arresting in the lower-right
  region without reaching the corner**. Ambati p.398: *"no further evolution
  of the phase-field in the lower-right corner is possible … the subsequent
  behavior corresponds to the linearly elastic response of the cracked
  specimen clamped at the undamaged lower-right portion of the boundary."*
  Arc length at arrest ≈ 0.55 mm.
- **Shear stress** `τ_xy` at the notch tip (Ambati Fig. 13): linear up to
  `u_x ≈ 13–14 µm` at peak force ~0.5 kN, then a sharp drop.
- **Energy balance**: `E_frac` starts at `Gc · Dn ≈ 1.35 J` (regularised
  notch baseline) and ramps toward `Gc · 0.55 mm ≈ 1.49 J` at arrest.

## Files

- `mesh.geo`, `geometry.yaml`         — geometry and label map.
- `input.yaml`                        — physics, regime, solver options.
- `boundary_conditions.yaml`          — clamped bottom, sliding top, D=1 pre-crack.
- `bc_generator.ipynb`                — notebook used to generate the 700-step
                                         displacement list.
- `non-regression.py`                 — diagnostic plots and pass/fail checks.
- `Allrun`, `Allclean`                — case-14-style drivers.

## Diagnostic outputs (in `output/`)

- `shear_response.png`                — τ_xy(tip) vs γ = u_x / Ly, with
                                         G·γ linear-elastic reference.
- `damage_evolution.png`              — max-D, normalised H, τ_xy at notch tip
                                         vs step, with the AT2 threshold line.
- `crack_profile_diagonal.png`        — D along the notch-tip → bottom-right
                                         diagonal across all loading steps.
- `energy_balance.png`                — E_el, E_frac, E_tot vs step with the
                                         Ambati arc-arrest reference.
- `damage_field.png`                  — 2D map of D at the final step
                                         (Ambati Fig. 12 hybrid reproducer);
                                         notch slit overlaid in cyan.
- `damage_evolution_panels.png`       — 4-panel D snapshot across the loading
                                         window.
- `stress_xy_field.png`               — σ_xy (Mode-II driver) at final.
- `stress_xx_field.png`               — σ_xx (shows the tensile/compressive
                                         diagonals of shear loading).
- `stress_vm_field.png`               — Von Mises at final.
- `crack_driving_force_field.png`     — H = (2·ℓc/Gc)·ψ⁺ at final.

## Running

```bash
./Allclean      # remove output and logs (mesh preserved)
./Allrun        # mesh -> z3st -> non-regression
```
