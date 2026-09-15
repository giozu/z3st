# box_notch_2D

AT2 phase-field plate (1 m x 0.5 m, high-carbon steel, lc = 2 mm) with a 100 mm deep V-notch at mid-span of the top edge. Left edge Clamp_x, bottom edge Clamp_y.

## Loading (changed 2026-09-15)

Until 2026-09-15 the right edge carried a single 180 MPa traction, which has no equilibrium past the peak once damage is coupled within the step (commit 75a44af). The right edge now carries a prescribed displacement ramp `Dirichlet_x` 0, 50, ..., 300, 400, 600, 800 um (10 steps), same supports and direction. `mesh.geo` adds a band of lc/2.5 elements below the notch tip. Relaxation off, `max_iters` 2000. `notched_plate_2D` runs the same specimen on a finer ramp.

## Checks

Final damage (1), largest sigma_xx over the history against the AT2 strength, crack offset from x = Lx/2, tracked peak reaction (integral of sigma_xx on the right edge).

## Status

Not re-blessed. The run of 2026-09-15 converged at every step (497 iterations, 2051 s): peak 21.15 MN/m at 200 um, then full separation along x = Lx/2 (offset 7.5 mm) and a residual force of 0.04 MN/m. The summary fails on `max_stress_xx`: 585 MPa against sigma_c = 173 MPa. The notch-tip stress exceeds the homogeneous AT2 strength on this mesh, so the check is not a valid analytic reference and needs an owner decision (track it, or compare a different quantity).
