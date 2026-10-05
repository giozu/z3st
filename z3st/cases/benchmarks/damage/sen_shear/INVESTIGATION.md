# Investigation plan: crack arrest in sen_shear

## The open point

Beyond the regularised notch the crack dissipates 4.669 - 1.358 = 3.31 J. The
arrest length of Ambati et al. (2015) Fig. 12d, 0.55 mm, corresponds to
`Gc · 0.55 mm = 1.49 J`, about 2.2 times less. Two explanations are possible:

- (A) the crack runs further than in the reference,
- (B) the crack has the reference length, and AT2 on this mesh overestimates the
  fracture energy.

A full run takes about 2 hours, so the steps below are ordered by cost.

## Step 0: measurements on the stored run (no new run)

1. From the last VTU, extract the crack path as the ridge of `D >= 0.5` starting
   at the notch tip `(0.5 mm, 0.5 mm)`. Measure its length `L` and its end point.
2. Compute `G_eff = 3.31 J / L`. The discretisation overestimate expected for
   AT2 is about `Gc (1 + h / (2 lc))`, with `h` the mesh size along the path.
3. `L` of about 1.2 mm or more points to (A). `L` close to 0.55 mm with
   `G_eff` close to `2 Gc` points to (B).

## Step 1: setup against Ambati §4.2 (reading, no run)

| Item | This case | To check in Ambati |
|---|---|---|
| Energy split | `star_convex`, `gamma_star: 5.0` | the split from which the hybrid formulation (Eq. 27) takes psi+, expected to be the spectral split of Miehe |
| Mesh along the path | 0.8 um at the two notch points only, graded linearly to 10 um at the corners, so several um along the arc, above `lc = 4 um` | mesh size in the crack region |
| Boundary conditions | left and right edges free | constraints on the left and right edges (some versions of the Miehe shear test constrain `u_y` there) |
| Load and length scale | `u_x` up to 30 um, `lc = 4 um` | final displacement and `lc` |

Also check where the arrest length of 0.55 mm comes from, since it is read from
a figure.

## Step 2: cheap tests before any 2-hour run

- Shortened load path: same mesh and parameters, fewer steps between 8 and
  15 um (700 steps of 0.01 um at present). If the crack path does not depend on
  the load step, every later test becomes cheaper.
- Split: `sweep_gamma.sh` runs `gamma_star` = 0, 1 and 5. Add the spectral split
  with the hybrid constraint, the setup of Ambati.

## Step 3: mesh along the path

Refine to `h <= lc / 2` along the expected path with a Gmsh `Box` or `Distance`
field on the arc, not only at the notch tip, and check that `G_eff · L`
converges.

## Step 4: outcome

- If the crack arrests as in Ambati: update the case, the gold (owner approval)
  and the README, and close the open point.
- If it does not: the README states the cause found (split, mesh or boundary
  conditions) with the measured numbers, and the case stays a regression test.
