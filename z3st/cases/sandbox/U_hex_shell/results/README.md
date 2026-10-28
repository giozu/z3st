# results

One subfolder per problem; each has `geometry_2d.png`, `geometry_3d.png` and a README with hot/cold sides and the statistical-zero flags.

- [inner_q_flux](inner_q_flux/README.md) — Radial through-wall gradient: inner heat flux q, outer T_o fixed. D6 symmetric: no bowing.
- [bowing_flat-flat_T_asymmetries](bowing_flat-flat_T_asymmetries/README.md) — Flat-to-flat bowing: Dirichlet T on the nord (cold) and sud (hot) flats, free duct.
- [bowing_corner-corner_IAEA_VP1/z-gradient](bowing_corner-corner_IAEA_VP1/z-gradient/README.md) — Corner-to-corner bowing, IAEA VP1 WITH axial gradient: T = 400 °C below the core, ramp across the core (1.5-2.5 m), cold/hot corner 500/550 °C above; clamped bottom. Comparison with the published benchmark and with Ansys when exported (same option).
- [bowing_corner-corner_IAEA_VP1/z-uniform](bowing_corner-corner_IAEA_VP1/z-uniform/README.md) — Corner-to-corner bowing, IAEA VP1 step 1 WITHOUT axial gradient (Z_UNIFORM): the core-top section (cold/hot corner 500/550 °C) over the whole height; clamped bottom. Comparison with Ansys (same option).
- [bowing_corner-corner_IAEA_VP2/z-gradient](bowing_corner-corner_IAEA_VP2/z-gradient/README.md) — Corner-to-corner bowing, IAEA VP2 (IWGFR/75 example 2) WITH axial gradient: as VP1, section non-linear (corners 1, 2, 6 at 500 °C, 3, 4, 5 at 550 °C above the core, linear between corners); clamped bottom. Comparison with the published benchmark and the CRP codes.
