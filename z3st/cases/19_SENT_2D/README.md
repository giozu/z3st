# Case 19 — SENT 2D (Single-Edge Notched Tension, plane strain)

Reproduces the **Ambati et al. 2015** (Comput Mech 55:383–405) §4.1 SENT-tension
benchmark using the AT2 hybrid phase-field formulation (their Eq. 27).

## Geometry & loading

- 1 mm × 1 mm square plate, 2D plane strain (`regime: 2d`).
- Horizontal notch of length `Dn = 0.5 mm` running from the left edge to the
  centre, at `y = Ly/2`. Modelled as a zero-width slit (curve traversed twice
  in the curve loop), with `D = 1` Dirichlet along the notch to anchor the
  crack from step 0.
- Bottom edge clamped; top edge displaced upward in 100 steps of 0.01 µm
  (from 5.01 µm to 6.0 µm). This window spans the crack-initiation transition
  of the Ambati Fig. 6 force-displacement curve.

## Material (`materials/high_carbon_steel.yaml`)

| Parameter | Value |
|---|---|
| `E`  | 210 GPa |
| `nu` | 0.3 |
| `Gc` | 2700 J/m² |

AT2 critical stress (derived in `spine.py` at load time):
`σ_c = sqrt(27·Gc·E / (256·ℓ_c)) ≈ 3.87 GPa` at `ℓ_c = 4 µm`.

## Mesh

- True 2D triangulation (`gmsh -2`), graded:
  - `h_fine = ℓ_c / 5 = 0.8 µm` at the notch tip (h/ℓ_c = 0.2 — Borden/Miehe ideal).
  - `h_coarse = Lx / 75 ≈ 13.3 µm` in the bulk.
- This is the key difference from the older `19_single-edge_notched_tension_test/`
  case, whose mesh was generated in 3D mode and produced a thin extruded layer.

## Phase-field formulation

- **AT2** crack-density functional, `γ = ½(D²/ℓ_c + ℓ_c |∇D|²)`.
- **Amor split** (volumetric/deviatoric) for the elastic-energy decomposition
  into `ψ⁺` and `ψ⁻`.
- **Hybrid constraint** (Ambati 2015): in cells where `ψ⁻ > ψ⁺`, the
  contribution of `ψ⁺` to the history field `H` is zeroed — suppresses crack
  growth in compression while keeping the mechanical block linear in `u`.

## Solver choices

- `mechanical.linear_solver: direct_mumps` — once cracks open and `g(D) → K ≈
  10⁻⁶` in cracked cells, the mechanical operator picks up 6–7 orders of
  magnitude heterogeneity that BoomerAMG handles poorly. MUMPS is robust and
  fast at this 2D system size (~10 k DOFs).
- `damage.linear_solver: iterative_hypre` — AT2's mass coefficient is `(H + 1)`,
  much milder than AT1's `2H`, so BoomerAMG is fine here.

## Expected results (Ambati Fig. 5–6)

- **Damage field**: a straight horizontal mode-I crack extending from the notch
  tip `(0.5 mm, 0.5 mm)` to the right edge `(1.0 mm, 0.5 mm)`.
- **Force-displacement curve**: linear up to `u_y ≈ 5.5 µm`, sharp peak around
  `5.5–6 µm`, near-vertical brittle drop as the crack traverses the ligament.
- **Crack profile** `D(y)` at mid-ligament: matches the AT2 analytical
  `exp(-|y - y_c| / ℓ_c)` decay around the crack centre.
- **Energy balance**: `E_frac` saturates near `Gc · (Lx - Dn) · 1 m = 1.35 J`
  per unit out-of-plane thickness once the crack fully traverses the ligament.

## Files

- `mesh.geo`, `geometry.yaml`         — geometry and label map.
- `input.yaml`                        — physics, regime, solver options.
- `boundary_conditions.yaml`          — clamped bottom, ramped top, D=1 pre-crack.
- `non-regression.py`                 — diagnostic plots and pass/fail checks.
- `Allrun`, `Allclean`                — case-14-style drivers.

## Diagnostic outputs (all in `output/`)

- `stress_strain_curve.png`           — σ_yy vs ε_yy at the notch tip,
                                         overlaid with the plane-strain
                                         analytical line.
- `damage_evolution.png`              — max-D, normalised H, and σ_yy at the
                                         notch tip vs step (the σ_c reference
                                         line is the AT2 analytical threshold).
- `crack_profile_evolution.png`       — D(y) at x = Dn/2 across all loading
                                         steps, compared with the AT2
                                         analytical decay `exp(-|y - y_c|/ℓ_c)`.
- `energy_balance.png`                — E_el, E_frac, E_tot vs step, with the
                                         `Gc · Lx` saturation reference.
- `damage_field.png`                  — ParaView-style 2D map of D at the final
                                         step (Ambati Fig. 5 reproducer); notch
                                         slit marked in cyan.
- `damage_evolution_panels.png`       — 4-panel D snapshot across the loading
                                         window.
- `stress_yy_field.png`               — σ_yy (Mode-I tensile driver) at final.
- `stress_xx_field.png`               — σ_xx (transverse; compressive lobes at
                                         the crack tip indicate shielding).
- `stress_vm_field.png`               — von Mises equivalent stress at final.
- `crack_driving_force_field.png`     — non-dimensional H = (2ℓ_c/Gc)·ψ⁺ at
                                         final, showing where damage is being
                                         driven.
