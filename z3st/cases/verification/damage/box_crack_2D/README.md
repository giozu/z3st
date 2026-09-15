# box_crack_2D

AT2 phase-field plate (1 m x 0.5 m, high-carbon steel, lc = 2 mm) with an edge pre-crack D = 1 on 0 < x < 0.25 m at y = Ly/2. Bottom edge Clamp_y, right edge Clamp_x.

## Loading (changed 2026-09-15)

Until 2026-09-15 the top edge carried a single 200 MPa traction. Since the damage iterate is coupled into the mechanical solve within the step (commit 75a44af), no equilibrium exists under traction control past the peak and the run diverged. The top edge now carries a prescribed displacement ramp `Dirichlet_y` from 0 to 400 um in 17 steps (10 um steps to 100 um, then 50 um steps), with the same supports and loading direction. `mesh.geo` adds a band of elements of size lc/2.5 along the expected crack path, since the crack cannot propagate through the 10 mm background elements. Relaxation is off (`relax_u = relax_D = 1`, plain alternate minimisation) and `max_iters` is 2000.

## Checks

`non-regression.py` reads the per-step VTUs: maximum damage (1 on the pre-crack), the largest sigma_yy over the history against the AT2 strength, the AT2 energy of the prescribed crack against Gc*Dn, the offset of the grown crack from y = Ly/2, and tracked peak reaction force (integral of sigma_yy on the top edge), displacement at peak, final/peak force ratio, crack tip position and localisation length.

## Status

Not re-blessed. On 2026-09-15 the refined run passed the peak (35.7 MN/m at 100 um on the coarse mesh) but the crack-jump step (12) did not converge within the 1.5 h cap (more than 1000 staggered iterations). The gold still holds the traction-controlled values.
