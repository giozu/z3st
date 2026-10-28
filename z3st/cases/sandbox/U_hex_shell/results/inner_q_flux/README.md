# inner_q_flux

Radial through-wall gradient: inner heat flux q, outer T_o fixed. D6 symmetric: no bowing.

Geometry and coordinates: `geometry_2d.png` (section at the along-s heights) and `geometry_3d.png`. s runs at mid-wall from V1 counter-clockwise; ξ through the wall, outward.

## zmaxvm (z = 12.5 mm)

```
zmaxvm, z = 12.5 mm: every flat at the same T (spread 3.2e-02 K, D6 symmetric): no hot or cold flat/vertex
through the wall: inner face 228.72 °C, outer face 226.85 °C (hot face: inner)
```

No along-s panel at the statistical zero (< 1e-4).

## zmid (z = 487.5 mm)

```
zmid, z = 487.5 mm: every flat at the same T (spread 3.2e-02 K, D6 symmetric): no hot or cold flat/vertex
through the wall: inner face 228.72 °C, outer face 226.85 °C (hot face: inner)
```

**FLAG convergence/statistics:** ε_zz: constant along s, variation 1.3e-10 (5.4e-08 relative) — y-axis zoomed on solver noise  

## Bowing figures

**FLAG convergence/statistics** `*vtb_bowing.png`: max |w| = 8.5e-12 mm ≈ 0  
**FLAG convergence/statistics** `*curvature.png`: max |kappa| = 7.4e-14 1/m ≈ 0  

## Files

- `3d_along_s_zmaxvm.png`
- `3d_along_s_zmid.png`
- `3d_curvature.png`
- `3d_displacement_z.png`
- `3d_max_stress.png`
- `3d_stress_z.png`
- `3d_summary.json`
- `3d_vtb_bowing.png`
- `3d_vtb_temperature_bowing.png`
- `geometry_2d.png`
- `geometry_3d.png`
