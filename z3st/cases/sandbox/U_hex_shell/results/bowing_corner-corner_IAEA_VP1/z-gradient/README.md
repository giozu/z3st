# bowing_corner-corner_IAEA_VP1/z-gradient

Corner-to-corner bowing, IAEA VP1 WITH axial gradient: T = 400 °C below the core, ramp across the core (1.5-2.5 m), cold/hot corner 500/550 °C above; clamped bottom. Comparison with the published benchmark and with Ansys when exported (same option).

Geometry and coordinates: `geometry_2d.png` (section at the along-s heights) and `geometry_3d.png`. s runs at mid-wall from V1 counter-clockwise; ξ through the wall, outward.

## inlet core (z = 1500.0 mm)

```
inlet core, z = 1500.0 mm: T uniform (400.0 °C): no hot or cold flat/vertex
```

**FLAG convergence/statistics:** ε_ss: max |v| = 1.4e-05 ≈ 0  
**FLAG convergence/statistics:** ε_zz: max |v| = 1.3e-05 ≈ 0  
**FLAG convergence/statistics:** ε_nn: max |v| = 9.4e-06 ≈ 0  

## mid core (z = 2000.0 mm)

```
mid core, z = 2000.0 mm: hot flat normal 300° (472.3 °C), cold flat normal 120° (453.1 °C); hot vertex V5 (475.0 °C), cold vertex V2 (450.0 °C)
through the wall: inner face 462.50 °C, outer face 462.50 °C (hot face: none, ΔT = 0 K negligible against the 25.6 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## outlet core (z = 2500.0 mm)

```
outlet core, z = 2500.0 mm: hot flat normal 300° (544.5 °C), cold flat normal 120° (506.2 °C); hot vertex V5 (550.0 °C), cold vertex V2 (500.0 °C)
through the wall: inner face 525.00 °C, outer face 525.00 °C (hot face: none, ΔT = 4.5e-13 K negligible against the 51.2 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## tip (z = 4000.0 mm)

```
tip, z = 4000.0 mm: hot flat normal 300° (544.5 °C), cold flat normal 120° (506.2 °C); hot vertex V5 (550.0 °C), cold vertex V2 (500.0 °C)
through the wall: inner face 525.00 °C, outer face 525.00 °C (hot face: none, ΔT = 4.5e-13 K negligible against the 51.2 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## Files

- `IAEA_VP1_z-gradient_along_s_normal.png`
- `IAEA_VP1_z-gradient_along_s_strain_sz.png`
- `IAEA_VP1_z-gradient_along_s_stress.png`
- `IAEA_VP1_z-gradient_along_y.png`
- `IAEA_VP1_z-gradient_corners_z.png`
- `IAEA_VP1_z-gradient_curvature.png`
- `IAEA_VP1_z-gradient_displacement_z.png`
- `IAEA_VP1_z-gradient_max_stress.png`
- `IAEA_VP1_z-gradient_stress_z.png`
- `IAEA_VP1_z-gradient_summary.json`
- `IAEA_VP1_z-gradient_vs_Ansys_along_s_normal.png`
- `IAEA_VP1_z-gradient_vs_Ansys_along_s_strain_sz.png`
- `IAEA_VP1_z-gradient_vs_Ansys_along_s_stress.png`
- `IAEA_VP1_z-gradient_vs_Ansys_curvature.png`
- `IAEA_VP1_z-gradient_vs_Ansys_differences.png`
- `IAEA_VP1_z-gradient_vs_Ansys_displacement_z.png`
- `IAEA_VP1_z-gradient_vs_Ansys_max_stress.png`
- `IAEA_VP1_z-gradient_vs_Ansys_stress_z.png`
- `IAEA_VP1_z-gradient_vs_Ansys_summary.json`
- `IAEA_VP1_z-gradient_vs_Ansys_vtb_bowing.png`
- `IAEA_VP1_z-gradient_vs_benchmark.png`
- `IAEA_VP1_z-gradient_vs_benchmark_fig2.png`
- `IAEA_VP1_z-gradient_vs_benchmark_summary.json`
- `IAEA_VP1_z-gradient_vs_benchmark_tables.png`
- `IAEA_VP1_z-gradient_vs_codes.csv`
- `IAEA_VP1_z-gradient_vs_codes.png`
- `IAEA_VP1_z-gradient_vtb_bowing.png`
- `IAEA_VP1_z-gradient_vtb_temperature_bowing.png`
- `geometry_2d.png`
- `geometry_3d.png`
