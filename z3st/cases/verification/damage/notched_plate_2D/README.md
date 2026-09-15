# notched_plate_2D

AT2 phase-field plate (1 m x 0.5 m, high-carbon steel, lc = 2 mm) with a 100 mm deep V-notch at mid-span of the top edge. Left edge Clamp_x, bottom edge Clamp_y.

## Loading (changed 2026-09-15)

Until 2026-09-15 the right edge carried a traction ramp to 300 MPa, which has no equilibrium past the peak once damage is coupled within the step (commit 75a44af). The right edge now carries a prescribed displacement ramp `Dirichlet_x` from 0 to 800 um in 21 steps (20 um steps to 300 um, then 100 um steps). `mesh.geo` adds a band of lc/2.5 elements below the notch tip. Relaxation off, `max_iters` 2000.

## Checks

The old metrics (u, sigma_xx and eps_xx at the final traction) were replaced by: crack offset from x = Lx/2 over Lx (1 %), and tracked peak reaction (integral of sigma_xx on the right edge), displacement at peak, final/peak force ratio, final crack depth and final fracture energy.

## Reference run (2026-09-15, gold)

Converged at every step, 867 staggered iterations, 3135 s. Peak 23.17 MN/m at 220 um, the crack runs to the bottom edge (depth 500 mm) along x = Lx/2 (offset 1.5 mm) in steps 12-13, final force 0.6 % of peak.
