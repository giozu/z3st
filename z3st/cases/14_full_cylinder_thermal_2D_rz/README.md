# Axisymmetric thermal-shock verification (UO2 pellet, 2D r-z)

## z3st verification case

### Role

This case is a **thermal + linear-elastic verification** of z3st in
axisymmetric (r, z) mode, using the McClenny et al. (JNM 565, 2022)
UO2 pellet quench parameters as a *thermal* loading scenario. **Damage
is intentionally disabled.**

### Why no damage in the axisymmetric variant

McClenny's experiment was **explicitly designed** to be azimuthally
asymmetric: only 1/6 of the lateral surface area was placed in contact
with the cold bath, and an alumina spacer at the capsule bottom was
used to *eliminate* axial thermal contact, so heat transfer was
deliberately confined to the radial direction (McClenny et al., Fig. 4
and §3). The fracture pattern observed (Fig. 8: two long radial cracks
plus a fan of shorter surface cracks) is a direct consequence of that
60-deg azimuthal asymmetry.

Axisymmetric (r, z) modeling mathematically prohibits any variation in
the azimuthal direction: every field is a function of (r, z) alone.
Consequently:

- The 60-deg contact wedge **cannot be represented**. The only
  physically valid axisymmetric idealization is full-circumference
  cooling (Dirichlet `T_quench` on the entire `r = Ro` lateral
  surface), which is *not* the experimental scenario.
- Any "single-face quench" axisymmetric variant (e.g. cooling only the
  top face) would *contradict* McClenny's design, since the alumina
  spacer was specifically there to avoid axial gradients.
- With damage enabled, full-circumference cooling can only produce an
  axisymmetric annular damage band at `r = Ro` extending the full
  axial height -- not the discrete radial cracks of the experiment.
  The staggered solver also becomes ill-conditioned once the band
  reaches `D ~ 1` (the effective stiffness `g(D) ~ K = 1e-6` makes
  the mechanical block singular).

For these reasons, **damage is disabled in this case** and it serves
only as:

1. a verification of the axisymmetric thermal solver against the
   analytic Bessel-series solution for transient cylindrical cooling
   under uniform Dirichlet on `r = Ro`;
2. a verification of the linear thermo-elastic stress field
   (sigma_rr, sigma_zz, sigma_theta_theta);
3. a sanity check that axisymmetric integration weights `2*pi*r` are
   applied consistently in the energy and stress integrals.

For the actual McClenny Fig. 8 cracking reproduction, use:

- `14_full_cylinder_cracking_2D_xy/` -- 2D Cartesian disc, plane strain,
  60-deg contact arc -- the closest 2D analogue of the experiment;
- `14_full_cylinder_cracking/` -- 3D, gold-standard model.

### Geometry, material, BCs

| Parameter | Value |
|---|---|
| Pellet radius `R` | 10 mm |
| Axial height `H` | 10 mm |
| Material | UO2 (E = 358 GPa, nu = 0.23, k = 5 W/m-K, alpha = 1e-5 /K) |
| `T_initial` | 1023.15 K (= 750 deg C) |
| `T_quench` | 263.15 K (= -10 deg C), Dirichlet on the **entire** lateral surface `r = Ro` |
| Top, bottom, axis | natural Neumann (zero flux) thermal |
| Mechanical | `Clamp_y` on bottom, `Clamp_x` on axis |
| Damage | **disabled** |
| Time window | 0.0001 to 0.05 s, n_steps = 50, dt = 1 ms |

### Verification targets

- The thermal radial profile `T(r, t)` at mid-height should match the
  analytic series sum `T(r, t) = T_quench + (T_initial - T_quench) *
  sum_n C_n J_0(lambda_n r) exp(-alpha_th lambda_n^2 t)` to within
  about 1 % (the `non-regression.py` script reports the L2 error).
- The radial stress profile under thermal load should match Kingery's
  thin-cylinder analytic estimate qualitatively: tensile at `r = Ro`,
  compressive at the centre, zero net radial force on each cross-section.

### Running

```bash
cd 14_full_cylinder_thermal_2D_rz/
./Allrun
```

`Allrun` chains `Allclean -> gmsh -2 mesh.geo -> python -m z3st -> non-regression.py`.

### References

- McClenny, Butt, Abdoelatef et al., *Experimentally validated multiphysics
  modeling of fracture induced by thermal shocks in sintered UO2 pellets*,
  JNM 565 (2022) 153719.
- Ambati, Gerasimov, De Lorenzis, *A review on phase-field models of brittle
  fracture and a new fast hybrid formulation*, Comput. Mech. 55 (2015) 383-405.
