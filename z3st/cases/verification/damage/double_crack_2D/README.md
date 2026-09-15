# double_crack_2D

AT2 phase-field plate (1 m x 0.5 m, steel with sigma_c = 600 MPa, Gc from the AT2 identity, lc = 2 mm) with two edge pre-cracks of 0.25 m at y = Ly/2. Bottom edge Clamp_y, right edge Clamp_x.

## Loading (changed 2026-09-15)

Until 2026-09-15 the top edge carried a traction ramp to 200 MPa, which has no equilibrium past the peak once damage is coupled within the step (commit 75a44af). The top edge now carries a prescribed displacement ramp `Dirichlet_y` from 0 to 600 um in 19 steps (20 um steps to 300 um, then 100 um steps). `mesh.geo` adds a band of lc/2.5 elements between the crack tips. Relaxation off, `max_iters` 2000.

## Checks

The old metrics (u, sigma_yy and eps_yy at the final traction) have no meaning under displacement control and were replaced by: AT2 energy of the two prescribed cracks against Gc*2*Dn (analytic, 1 %), offset of the grown cracks from y = Ly/2 over Ly (1 %), and tracked peak reaction (integral of sigma_yy on the top edge), displacement at peak, final/peak force ratio and final fracture energy.

## Reference run (2026-09-15, gold)

Converged at every step, 1219 staggered iterations, 3999 s. Peak 100.5 MN/m at 400 um, final force 0.5 % of peak after the ligament separates (step 18, 761 iterations). E_frac initial within 0.6 % of Gc*2*Dn, crack offset 4.45 mm.
