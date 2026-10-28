# bowing_corner-corner_IAEA_VP2/z-gradient

Corner-to-corner bowing, IAEA VP2 (IWGFR/75 example 2) WITH axial gradient: as VP1, section non-linear (corners 1, 2, 6 at 500 °C, 3, 4, 5 at 550 °C above the core, linear between corners); clamped bottom. Comparison with the published benchmark and the CRP codes.

Geometry and coordinates: `geometry_2d.png` (section at the along-s heights) and `geometry_3d.png`. s runs at mid-wall from V1 counter-clockwise; ξ through the wall, outward.

## inlet core (z = 1500.0 mm)

```
inlet core, z = 1500.0 mm: T uniform (400.0 °C): no hot or cold flat/vertex
```

**FLAG convergence/statistics:** ε_ss: max |v| = 1.5e-05 ≈ 0  
**FLAG convergence/statistics:** ε_zz: max |v| = 1.3e-05 ≈ 0  
**FLAG convergence/statistics:** ε_nn: max |v| = 9.9e-06 ≈ 0  

## mid core (z = 2000.0 mm)

```
mid core, z = 2000.0 mm: hot flat normal 300° (475.0 °C), cold flat normal 120° (450.0 °C); hot vertex V5 (475.0 °C), cold vertex V2 (450.0 °C)
through the wall: inner face 462.50 °C, outer face 462.50 °C (hot face: none, ΔT = 1.1e-13 K negligible against the 25 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## outlet core (z = 2500.0 mm)

```
outlet core, z = 2500.0 mm: hot flat normal 300° (550.0 °C), cold flat normal 120° (500.0 °C); hot vertex V5 (550.0 °C), cold vertex V2 (500.0 °C)
through the wall: inner face 525.00 °C, outer face 525.00 °C (hot face: none, ΔT = 3.4e-13 K negligible against the 50 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## tip (z = 4000.0 mm)

```
tip, z = 4000.0 mm: hot flat normal 300° (550.0 °C), cold flat normal 120° (500.0 °C); hot vertex V5 (550.0 °C), cold vertex V2 (500.0 °C)
through the wall: inner face 525.00 °C, outer face 525.00 °C (hot face: none, ΔT = 3.4e-13 K negligible against the 50 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## Files

- `IAEA_VP2_z-gradient_along_s_normal.png`
- `IAEA_VP2_z-gradient_along_s_strain_sz.png`
- `IAEA_VP2_z-gradient_along_s_stress.png`
- `IAEA_VP2_z-gradient_along_y.png`
- `IAEA_VP2_z-gradient_corners_z.png`
- `IAEA_VP2_z-gradient_curvature.png`
- `IAEA_VP2_z-gradient_displacement_z.png`
- `IAEA_VP2_z-gradient_max_stress.png`
- `IAEA_VP2_z-gradient_stress_z.png`
- `IAEA_VP2_z-gradient_summary.json`
- `IAEA_VP2_z-gradient_vs_benchmark.png`
- `IAEA_VP2_z-gradient_vs_benchmark_summary.json`
- `IAEA_VP2_z-gradient_vs_benchmark_tables.png`
- `IAEA_VP2_z-gradient_vs_codes.csv`
- `IAEA_VP2_z-gradient_vs_codes.png`
- `IAEA_VP2_z-gradient_vtb_bowing.png`
- `IAEA_VP2_z-gradient_vtb_temperature_bowing.png`
- `geometry_2d.png`
- `geometry_3d.png`
