# bowing_corner-corner_IAEA_VP1/z-uniform

Corner-to-corner bowing, IAEA VP1 step 1 WITHOUT axial gradient (Z_UNIFORM): the core-top section (cold/hot corner 500/550 °C) over the whole height; clamped bottom. Comparison with Ansys (same option).

Geometry and coordinates: `geometry_2d.png` (section at the along-s heights) and `geometry_3d.png`. s runs at mid-wall from V1 counter-clockwise; ξ through the wall, outward.

## inlet core (z = 1500.0 mm)

```
inlet core, z = 1500.0 mm: hot flat normal 300° (544.5 °C), cold flat normal 120° (506.2 °C); hot vertex V5 (550.0 °C), cold vertex V2 (500.0 °C)
through the wall: inner face 525.00 °C, outer face 525.00 °C (hot face: none, ΔT = 4.5e-13 K negligible against the 51.2 K range)
```

No along-s panel at the statistical zero (< 1e-4).

## mid core (z = 2000.0 mm)

```
mid core, z = 2000.0 mm: hot flat normal 300° (544.5 °C), cold flat normal 120° (506.2 °C); hot vertex V5 (550.0 °C), cold vertex V2 (500.0 °C)
through the wall: inner face 525.00 °C, outer face 525.00 °C (hot face: none, ΔT = 4.5e-13 K negligible against the 51.2 K range)
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

- `IAEA_VP1_z-uniform_along_s_normal.png`
- `IAEA_VP1_z-uniform_along_s_strain_sz.png`
- `IAEA_VP1_z-uniform_along_s_stress.png`
- `IAEA_VP1_z-uniform_curvature.png`
- `IAEA_VP1_z-uniform_displacement_z.png`
- `IAEA_VP1_z-uniform_max_stress.png`
- `IAEA_VP1_z-uniform_stress_z.png`
- `IAEA_VP1_z-uniform_summary.json`
- `IAEA_VP1_z-uniform_vs_Ansys_along_s_normal.png`
- `IAEA_VP1_z-uniform_vs_Ansys_along_s_strain_sz.png`
- `IAEA_VP1_z-uniform_vs_Ansys_along_s_stress.png`
- `IAEA_VP1_z-uniform_vs_Ansys_curvature.png`
- `IAEA_VP1_z-uniform_vs_Ansys_differences.png`
- `IAEA_VP1_z-uniform_vs_Ansys_displacement_z.png`
- `IAEA_VP1_z-uniform_vs_Ansys_max_stress.png`
- `IAEA_VP1_z-uniform_vs_Ansys_stress_z.png`
- `IAEA_VP1_z-uniform_vs_Ansys_summary.json`
- `IAEA_VP1_z-uniform_vs_Ansys_vtb_bowing.png`
- `IAEA_VP1_z-uniform_vtb_bowing.png`
- `IAEA_VP1_z-uniform_vtb_temperature_bowing.png`
- `geometry_2d.png`
- `geometry_3d.png`
