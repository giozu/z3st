# UO2 pellet thermal-shock fracture - 2D transverse cross-section
## z3st application case (McClenny Fig. 8 reproducer)

### Role

This case is the **2D Cartesian (x, y) reproducer** of the McClenny
et al. (JNM 565, 2022) UO2 pellet thermal-shock experiment, modeled as
a transverse cross-section under plane strain. It is the natural 2D
analogue of the McClenny Fig. 8 (top row) result: a circular pellet
with a 60-deg cold-contact arc on its perimeter, producing discrete
radial cracks.

Unlike the axisymmetric `14_full_cylinder_thermal_2D_rz/` variant
(which can only produce an unphysical annular damage band), this case
can break azimuthal symmetry and capture the localized tensile-hoop
stress pattern that drives the experimentally observed radial cracks.

### Why this 2D idealization is faithful to the experiment

McClenny's experimental setup was deliberately designed so that:
- The cold bath was in contact with **only 1/6 of the lateral surface
  area** of the pellet (Fig. 4), creating a 60-deg azimuthal contact
  wedge.
- An **alumina spacer at the capsule bottom** eliminated axial thermal
  contact (Section 3), so heat transfer was confined to the radial
  direction.

The "no axial gradient" condition is exactly the plane-strain
assumption: every transverse cross-section sees the same physics.
The plane-strain (x, y) cross-section of the pellet, with a 60-deg
contact arc on its perimeter, therefore captures the dominant physics
of the experiment at much lower compute cost than the full 3D
simulation. McClenny themselves use this 2D representation for their
parametric study (Fig. 8 top row, Fig. A.13, Fig. A.14).

### Geometry choice: half-disc with mirror symmetry

The contact wedge in the model is centered around the +x axis and
spans -30 deg to +30 deg (full 60 deg). This makes the y = 0 plane a
mirror symmetry plane: the loading and the geometry are symmetric
across it. The mesh therefore covers only the **upper half** of the
disc (y >= 0); the lower half is implicit by mirror symmetry.

This halves the compute cost and lets `Clamp_y` on the diameter act
as the symmetry boundary condition (removing both rigid-body
y-translation and rotation around z). A small 50-um pin segment near
(-R, 0) carries a `Clamp_x` to remove rigid x-translation.

### Geometry, material, BCs

| Parameter | Value |
|---|---|
| Pellet radius `R` | 10 mm |
| Domain | upper half-disc (y >= 0); contact arc 0 to +30 deg in the upper half (= 60 deg full) |
| Material | UO2 from `materials/uo2.yaml`: E = 358 GPa, nu = 0.23, alpha = 1e-5 /K, sigma_c = 2 GPa |
| `T_initial` | 1023.15 K (= 750 deg C), uniform |
| `T_quench` | 263.15 K (= -10 deg C), Dirichlet on the 60-deg contact arc |
| Insulated rest of perimeter | natural Neumann (zero heat flux) |
| Symmetry plane y = 0 | `Clamp_y` (mirror BC) |
| Pin segment at (-R, 0) | `Clamp_x` (50 um, removes rigid x-translation) |
| Damage | AT1, Ambati hybrid, Amor split, `lc = 50 um`, hybrid_constraint = true |
| Time window | 0.0001 to 0.1 s, n_steps = 100, dt = 1 ms |

`sigma_c = 2 GPa` is calibrated above the bulk-artifact threshold of
plane strain (about 5.6 MJ/m^3 from the blocked z-thermal-expansion,
suppressed by the regime-aware eigenstrain fix in `damage_model.py`)
and below the peak surface tensile-strain energy at the rim, so that
damage initiates only along the cold contact arc. The corresponding
`Gc` is auto-derived in `spine.py` (Gc = (8/3) lc sigma_c^2 / E ~ 1490 J/m^2).

### Phase-field formulation

Z3ST uses the **Ambati hybrid formulation** (Comput. Mech. 55 (2015)
383-405, Eq. 27): linear isotropic stress degradation `sigma = (1-D)^2
* dPsi0/de`, AT1 surface energy with Amor (volumetric/deviatoric)
elastic-energy split, hybrid constraint setting `H = 0` in compression-
dominated cells. The damage driving force is evaluated on the
**elastic strain** `eps_el = eps(u) - alpha (T - T_ref) I` so that
uniform thermal expansion doesn't appear as a damage driver. In 2D
plane strain, the z-component of the eigenstrain is suppressed (it
would otherwise create a uniform compressive bulk artifact via the
geometrically-blocked thermal expansion in z; see `damage_model.py::
_thermal_eigenstrain` docstring).

McClenny instead uses the Miehe anisotropic formulation with viscous
Allen-Cahn evolution; reproducing their crack pattern with the hybrid
formulation is the methodological contribution of this case.

### Expected results

- Discrete radial cracks emerge from the cold contact arc within the
  -30 deg to +30 deg wedge.
- Per McClenny (p. 7): "two major (longer) radial cracks" plus a fan
  of shorter surface cracks.
- The long cracks **do not reach the centre**: they are arrested by
  the central compression zone (the Amor split sends the compressive
  hoop strain into psi_neg, and the hybrid constraint sets H = 0 in
  those cells).
- Crack initiation around `t ~ 1e-2 s` (McClenny p. 8: "the cracks
  immediately appear on the pellet outer surface at about 10^-2 s
  right after the instantaneous drop in temperature").
- `E_frac` should grow rapidly once cracks initiate; `E_el` should
  ramp up before the crack threshold and partially relax after.

Crack bands will appear ~50x wider than McClenny's because of the lc
coarsening (50 um vs 1 um in the reference). The topology and timing
are the diagnostic targets, not the band width.

### Running

```bash
cd 14_full_cylinder_cracking_2D_xy/
./Allrun
```

### Outputs (in `output/`)

- `damage_field.png` -- ParaView-style 2D colormap of `D` on the half-disc,
  with the cold contact arc highlighted; the McClenny Fig. 8 (top, right)
  reproduction.
- `temperature_field.png` -- same paraview-like rendering of `T` at the
  final time.
- `stress_vm_field.png` -- von-Mises equivalent stress at the final time
  (clipped to the 99th percentile to keep the colorbar informative).
- `stress_hoop_field.png` -- hoop stress sigma_theta_theta with a
  symmetric red/blue colormap (red = tensile = crack-driver).
- `damage_angular.png` -- D_max(theta) along the outer ring, peaks count
  individual radial cracks.
- `thermal_shock_results.png` -- T(r) along contact midline, T(t) at
  reference points, D(r) within the wedge.
- `stress_evolution.png` -- sigma_rr(r), sigma_tt(r) along the contact
  midline.
- `energy_balance.png` -- E_el(t), E_frac(t).

### References

- McClenny et al., JNM 565 (2022) 153719.
- Ambati, Gerasimov, De Lorenzis, Comput. Mech. 55 (2015) 383-405.
