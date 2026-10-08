# Investigation: sen_shear against Ambati et al. (2015)

## Why

The README of this case recorded an open point: beyond the regularised notch
the crack dissipates 4.669 - 1.358 = 3.31 J/m, against `Gc · 0.55 mm = 1.49 J/m`
for an arrest length of 0.55 mm attributed to Ambati et al. (2015) Fig. 12d.
The aim is to find out whether the case reproduces the hybrid AT2 shear test
of Ambati §4.2 and, where it does not, why. A full run takes about 2 hours, so
the work is ordered by cost: measurements on the stored run first, new runs
after.

## Status

- The 0.55 mm arrest target is not from the paper. In Ambati Fig. 12d the
  hybrid crack reaches the bottom edge near the lower-right corner. The open
  point is restated: this case matches initiation and the first 0.2 mm of the
  crack, then the crack slows down and the response stays stiffer.
- The fracture energy per unit length grows from 1.6 to 3 Gc along the crack
  path and follows the local mesh size (3 to 9 um, with `lc = 4 um`). Ambati
  used a uniform mesh with h about 7 um, coarser than most of this path, and
  his crack reaches the bottom edge. The suspected cause is therefore the
  grading of the mesh, which makes the crack more expensive as it moves into
  coarser cells, more than its size. On a uniform 7 um mesh, close to
  Ambati's, the crack path matches Ambati and Hirshikesh et al. and reaches the
  bottom edge (section 8). The force does not match: the peak is 0.82 kN
  against about 0.49 kN, and the specimen fails completely at 30 um where
  Ambati's force rises again. 136 of the 796 steps did not converge.
- What is shown: the graded mesh slows the crack down; on a uniform mesh the
  path matches. What is not shown: that a different mesh makes the force
  match. The remaining differences from Ambati that do not depend on the mesh
  are the notch (`D = 1` on a line here, an initial history field in Ambati),
  the energy split, the implementation of the hybrid constraint, the staggered
  strategy (tolerance here, a fixed number of iterations in Ambati) and the
  element type (P1 triangles here, bilinear quadrilaterals in Ambati).
  Ambati's own mesh (h about 7 um, lc = 4 um) is not shown to be converged, so
  the aim is to reproduce his test qualitatively and to measure where the
  differences come from, one factor at a time, on his mesh.
- On Ambati's mesh (uniform bilinear quadrilaterals, 144 x 144) the crack does
  not start from the notch at all within 30 um; damage forms instead in bands
  along the top and bottom edges (section 9). The mesh alone does not explain
  the differences; the notch representation and the energy split are the next
  factors.
- On every refined mesh tried, the staggered loop stops converging within the
  first 8 um of loading. The cause is an exact period-2 oscillation of the
  damage, of amplitude about 1e-5, at nodes along the notch line behind the
  tip. It keeps `||dD||/||D||` near 1e-6, above the damage `stag_tol` of 1e-7.
- The switch is identified for the first non-converged steps: the hybrid
  test `psi- > psi+` of a single compressed cell next to the notch face flips
  at every iteration, and the history and the damage follow it.
- The loss of convergence is not limited to refined meshes: the uniform 7 um
  mesh has 136 non-converged steps out of 796.
- A looser damage tolerance does not help. With `stag_tol: 1e-5` for the damage
  the stored mesh itself loses convergence (29 of 303 steps), while with 1e-7 it
  converges in all 796 steps. The fix belongs in the code.
- No reference iterates the hybrid staggered loop to a tight tolerance on this
  test: Ambati uses a fixed 1, 2, 4 or 8 staggered iterations per step.
  Vicentini et al. (2024), who introduced the star-convex split, report loss of
  iterative convergence under mesh refinement with that split and leave it
  open. The combination used here (hybrid with star-convex, `gamma_star: 5`)
  has no published precedent.

## What has been done

### 1. Crack length and fracture energy on the stored run (no new run)

Crack polyline from the `D > 0.9` nodes outside the notch, binned in y, and
the AT2 energy density `Gc/(2 lc) (D^2 + lc^2 |grad D|^2)` integrated on the
P1 triangles of `mesh.msh`. The integral reproduces the solver's `E_frac` at
the last step, 4.669 J/m.

Crack length beyond the notch `L = 0.491 mm`, end point (0.732, 0.114) mm. The
crack is not longer than the reference.

| Region | E_frac (J/m) | Sharp-crack value |
|---|---|---|
| notch band | 1.492 | `Gc · 0.5 mm` = 1.350 |
| crack band (within 10 lc of the path) | 2.727 | `Gc · L` = 1.326 (ratio 2.06) |
| diffuse halo | 0.450 | 0 |

The crack band is born wide and does not widen after the tip has passed: at
y = 0.423 mm its width (`D >= 0.5`) is 13.5 um from step 500 to step 795, about
twice the ideal AT2 profile. The energy per unit length follows the mesh:

| y (mm) | median h (um) | h / lc | 1 + h / (2 lc) | measured |
|---|---|---|---|---|
| 0.423 | 3.3 | 0.82 | 1.41 | 1.58 Gc |
| 0.346 | 5.1 | 1.27 | 1.64 | 1.80 Gc |
| 0.268 | 6.1 | 1.52 | 1.76 | 1.95 Gc |
| 0.191 | 7.2 | 1.79 | 1.89 | 2.95 Gc |

### 2. Setup and results against Ambati §4.2 (reading the paper)

- `lc = 4.0e-3 mm`, as in this case. Boundary conditions (Fig. 7b) as in this
  case: bottom clamped, top displaced horizontally, sides free.
- Uniform mesh of 20592 quadrilaterals, `h` about 7 um, increments of `1e-5 mm`
  throughout.
- Fig. 12, hybrid row: at `u = 0.030 mm` the crack reaches the bottom edge near
  the lower-right corner. No arrest at 0.55 mm.
- Fig. 13, hybrid curve, read from the plot: first peak about 0.49 kN at about
  12.5 um, minimum about 0.35 kN near 16 um, about 0.47 kN at 30 um.
- The split from which the hybrid model takes psi+ in the examples is not
  stated in the text read so far.

| | Ambati (hybrid) | This case |
|---|---|---|
| first peak | about 0.49 kN at 12.5 um | 0.480 kN at 11.4 um |
| post-peak minimum | about 0.35 kN near 16 um | 0.422 kN at 15 um |
| force at 30 um | about 0.47 kN | 0.550 kN |
| crack tip at 15 um | Fig. 12b | y = 0.31 mm, similar |
| crack tip at 20 um | near the bottom edge (Fig. 12c) | y = 0.24 mm |
| crack tip at 30 um | at the bottom edge (Fig. 12d) | y = 0.11 mm |

### 3. Mesh refined along the path: the staggered loop does not converge

Mesh with a Gmsh `Distance` + `Threshold` band of half-width 60 um along the
expected path, from the notch tip through (0.62, 0.25) and (0.73, 0.11) mm to
the bottom edge at x = 0.84 mm. Each run changes one item against the stored
settings (`relax_D: 0.8`, `relax_adaptive: false`, `hybrid_constraint: true`,
`split: star_convex`, `gamma_star: 5`, damage `stag_tol: 1e-7`). Staggered
iterations per step from step 0, X for 200 iterations without convergence:

| Run | Band h | Change | Iterations per step | Result |
|---|---|---|---|---|
| stored case | none (0.8 um at the notch, graded to 10 um) | none | 12 7 8 8 8 ... 6 | converges in all 796 steps |
| a | 2 um (lc/2) | none | 12 7 19 22 32 42 98 31 19 ... | X from step 23 (u = 6.6 um) |
| b | 2 um | `relax_D: 0.5` | 25 15 X | X from step 2 |
| c | 1 um (lc/4) | none | 12 7 8 10 12 19 19 22 29 32 X | X from step 11 (u = 6.0 um) |
| d | 1 um | `hybrid_constraint: false` | 12 7 8 10 12 19 19 22 27 33 45 63 X | X from step 13 |
| e | 2 um | increments 0.2 um to 5 um, 0.02 um to 8 um | 12 5 6 6 6 7 ... 9 9 9 X X X | X at u = 4.0-4.4 um |
| f | 2 um | `relax_adaptive: true` | 5 19 X X 97 X | worse |
| g | 2 um | `gamma_star: 0` | 12 7 8 10 13 21 22 28 37 46 78 X X X | X from u = 6.6 um |

Up to the loss of convergence the damage field matches the stored case (max D
outside the notch about 0.12 at (0.512, 0.498) mm). Element quality is
equivalent (minimum angle 33 deg against 31 deg). Relaxation, load increment,
hybrid constraint and split do not remove the problem. The point where it
appears changes between otherwise identical runs (step 2 in one run, step 23
in another), so it is sensitive to round-off.

### 4. Where the damage iterate oscillates

A wrapper around `Spine._damage_step`, in a debug copy of the case outside the
repository, saved `D` at every staggered iteration of two non-converged steps
on the 2 um band with the stored settings:

| Step | Nodes above 10 % of the max amplitude | Location | Distance from the notch tip | Amplitude of D |
|---|---|---|---|---|
| 2 | 55 | x 0.188-0.195 mm, y 0.4945-0.4995 mm (just below the notch) | 305-312 um | 1.5e-5 |
| 3 | 39 | x 0.388-0.393 mm, y 0.5004-0.5046 mm (just above the notch) | 107-112 um | 4.4e-5 |

The oscillation is an exact period-2 cycle (`D_k - D_(k-2)` is zero to float32
precision) at nodes of the regularised notch profile (D 0.58 to 0.91) along
the notch line, behind the tip, outside the refined band. About 50 nodes at
an amplitude of 1e-5 keep `||dD||/||D||` near 1e-6. An exact period-2 cycle
points to a switch that flips at every iteration: the sign of `tr(eps)` in the
split, or the hybrid test `psi- > psi+`, on the sheared notch faces where
`tr(eps)` is close to zero. Why the stored mesh does not show it is not known.

### 5. Looser damage tolerance

Damage `stag_tol` from 1e-7 to 1e-5, everything else as the stored settings,
full loading, two runs in parallel with 2 threads each, stopped after 2 hours:

| Run | Steps done | Non-converged steps | First at |
|---|---|---|---|
| stored mesh | 303 of 796 | 29 | u = 5.4 um (steps 7 to 19, then 87, 120, 241 to 296) |
| 2 um band | 279 of 796 | 60 | u = 5.0 um (steps 5 to 23, 81 to 88, 107 to 114, 244 to 278) |

With the looser tolerance the stored mesh, which converges in all 796 steps at
1e-7, loses convergence too. Stopping the iteration earlier leads to slightly
different states, which then fall into the same cycle. A looser tolerance is
not a workaround.

### 6. The switch that cycles

A second wrapper, around `Spine.update_history`, recorded at every staggered
iteration of steps 1 to 4 (2 um band, stored settings, damage `stag_tol`
1e-7), for the DG0 cells within 20 um of the notch: the history at step start,
the new value before the ratchet, the hybrid mask `psi- > psi+` and
`tr(eps)`. Steps 2, 3 and 4 did not converge.

| Step | Cell | Hybrid mask, last 6 iterations | tr(eps) | H of the cell |
|---|---|---|---|---|
| 2 | (0.1916, 0.4989) mm | 1 0 1 0 1 0 | -0.0045, constant sign | alternates, jump 2.4e-2 |
| 3 | (0.3907, 0.5007) mm | 0 1 0 1 0 1 | -0.013, constant sign | alternates, jump 0.157 |

In steps 2 and 3 one cell, and only one in the recorded band, changes its
hybrid mask at every iteration, in phase with the period-2 cycle of the damage
at the nodes around it. The cell is compressed (`tr(eps) < 0`) next to the
notch face, where `psi-` and `psi+` are close. When the mask is set, the new
history is zeroed and the ratchet returns `H` to its step-start value, the
local damage drops, the displacement changes slightly and at the next
iteration `psi- < psi+`. The ratchet `max(H_start, H_new)` and the sign of
`tr(eps)` never change in the recorded band. In step 4 the oscillation lies
outside the recorded band and its switch is not identified. Run d, without
the hybrid constraint, also lost convergence (step 13), so a second switch may
exist later in the loading.

### 7. What the literature says

Read for this investigation: Ambati et al. (2015), Bourdin et al. (2000),
Gerasimov and De Lorenzis (2016, 2019, 2022), Hirshikesh et al. (2019, Front.
Struct. Civ. Eng., FEniCS), Vicentini et al. (2024), Sidharth and Rao (2024).

Mesh and effective toughness.
- Bourdin et al. (2000) give only an order of magnitude: the excess surface
  energy is O(h/c), with c = lc/2 in the AT2 convention and h the radius of a
  circle inscribed in or containing an element, with no constant. The explicit
  form `Gc (1 + h/(4 eps))`, which with `eps = lc/2` gives `Gc (1 + h/(2 lc))`,
  is attributed to Bourdin et al. (2008) and has not been checked in that
  paper. The factor `1 + h/(2 lc)` in section 1 is an estimate; the measured
  growth of the energy with h is the result.
- Meshes used for this test: Ambati, uniform, 20592 quadrilaterals, h about
  7 um, `lc = 4 um` ("seems to be enough to eliminate the mesh-related
  effects", p. 397). Gerasimov and De Lorenzis (2016), `lc = 10 um`,
  `h < lc/2` in a pre-refined band. Hirshikesh et al. (2019), uniform
  triangles, `lc = 2 h` with `lc = 11 um`.

Hybrid constraint.
- Ambati Eq. 27c is a pointwise condition on d: "for all x: psi0+ < psi0- implies
  d := 0". Its numerical implementation is not described. Z3ST instead zeroes
  the new history contribution per DG0 cell. Ambati does not state which split
  psi0+ the hybrid runs use. Hirshikesh et al. give no code for the hybrid
  constraint either.

Staggered convergence.
- Ambati's tables and iteration study (Table 3, Figs. 14 and 15) use a fixed
  1, 2, 4 or 8 staggered iterations per step. His energy-based stopping
  criterion (Sect. 3.4) needs an energy functional, which the hybrid model
  does not have (p. 392).
- Gerasimov and De Lorenzis (2016) iterate the staggered loop to convergence
  (Tol 1e-4 on the phase-field residual), needing up to about 300 iterations
  per step, with no history field and no hybrid constraint. Stopping early
  underestimates the crack (their Fig. 10).
- Vicentini et al. (2024), Sect. 5.4: with the star-convex split "iterative
  convergence is lost when varying the mesh size" (h = lc/3 converges,
  h = lc/5 does not), an issue they leave open.
- None of these papers treats an on/off switch cycling inside the staggered
  loop. Constraint sets built from the previous step and held fixed within the
  step (Bourdin's crack set, used by Gerasimov and De Lorenzis) are the closest
  precedent for evaluating the hybrid mask once per step.
- Gerasimov and De Lorenzis (2019) replace the history field with a penalty on
  `<d - d_n>_-`, with `gamma_opt = (Gc/lc)(1/TOL_ir^2 - 1)` for AT2 (Eq. 50,
  `TOL_ir = 0.01`). They show that the history field thickens the damage band
  and overestimates the fracture energy. A penalty would remove the history
  field but not the hybrid switch.

### 8. Uniform mesh with h about 7 um, as Ambati

Stored settings (damage `stag_tol` 1e-7, star-convex `gamma_star: 5`, hybrid
constraint) on a uniform mesh of 47845 triangles with edges of about 7 um,
including the notch tip. 796 steps in 4 h 10 min, 136 of them not converged.

| | Ambati (hybrid, Fig. 12-13) | Uniform 7 um | Stored graded mesh |
|---|---|---|---|
| first peak | about 0.49 kN at 12.5 um | 0.822 kN at 15.8 um | 0.480 kN at 11.4 um |
| after the peak | about 0.35 kN, then rising | 0.317 kN at 17 um, 0.368 kN at 25 um | 0.422 kN at 15 um |
| force at 30 um | about 0.47 kN | 0.015 kN | 0.550 kN |
| crack tip at 20 um | near the bottom edge | about (0.78, 0.09) mm | y = 0.24 mm |
| crack at 30 um | reaches the bottom edge near the lower-right corner | reaches the bottom edge, then runs along it to the corner | y = 0.11 mm |

Hirshikesh et al. (2019) report the hybrid crack tip at about (0.77, 0.08) mm
at u = 21 um, read from their Fig. 11b. The path on the uniform mesh matches
both references, which supports the grading of the stored mesh as the cause
of the slowdown. Two differences remain. Initiation is late and the peak is
68 % above Ambati's, because the notch, imposed as `D = 1` on a line, is
resolved by 7 um cells at the tip; Ambati models the notch with an initial
history field. At 30 um the specimen fails completely, with the crack running
along the bottom edge, while Ambati's force rises again.

### 9. Ambati's mesh: uniform bilinear quadrilaterals

Stored settings on a structured mesh of 144 x 144 bilinear quadrilaterals
(20736, against Ambati's 20592), the notch being the line shared by the two
left blocks with `D = 1` imposed on it. 796 steps in 3 h 29 min, 117 of them
not converged.

| | Ambati | Uniform triangles, 7 um (section 8) | Uniform quadrilaterals, 7 um |
|---|---|---|---|
| first peak | about 0.49 kN at 12.5 um | 0.822 kN at 15.8 um | 0.928 kN at 17.8 um |
| force at 30 um | about 0.47 kN | 0.015 kN | 0.734 kN |
| crack from the notch tip | yes, to the bottom edge | yes, to the bottom edge | no |
| damage along the edges | not shown | small spots at the corners | bands along the top edge (x from about 0.65 to 1 mm) and the bottom edge |

With the element type and size of Ambati's mesh the case does not reproduce
his result: the notch tip does not initiate a crack and the specimen is damaged
along the loaded and clamped edges instead. Two observations, not yet tested:
with `D = 1` imposed only at the nodes of the notch line and cells of 7 um
(above lc = 4 um) the damage drops from 1 to 0 within one cell, so the notch is
represented weakly; and Ambati shows no damage at the edges, which may come
from the star-convex split with `gamma_star: 5` or from the hybrid constraint
near the constrained edges.

### 10. Energy split on Ambati's mesh

Same quadrilateral mesh and settings as section 9, one change each:

| Run | Split | Result |
|---|---|---|
| S1 | `amor` | 796 steps in 1 h 44 min, 21 not converged. Peak 0.916 kN at 17.8 um. No crack from the notch; damage along the top edge, complete failure along it from about 25 um (force zero) |
| S2 | `miehe` (spectral) | 71 of the first 74 steps not converged, about 1.6 min per step; stopped |

The split does not change the outcome of section 9: with Amor as with the
star-convex split the notch does not initiate a crack on this mesh and the
top edge fails instead. With the spectral split the loop almost never
converges.

A difference of implementation remains untested. Z3ST stores the history `H`
per cell (`self.Q`, DG0), with one value interpolated at the cell centre.
Ambati uses fully integrated bilinear quadrilaterals, so the history lives at
the 2 x 2 Gauss points. At the notch tip the strain is singular, and with
cells of 7 um (above lc = 4 um) one value at the cell centre underestimates
the peak of `psi+` near the tip more than four Gauss points do; on
quadrilaterals, whose cells have twice the area of the triangles, more so.
This is consistent with the three meshes: on the stored mesh (0.8 um at the
tip) the crack starts close to Ambati's load, on uniform triangles it starts
late, on uniform quadrilaterals it does not start and a constrained edge,
where the strain is high but smooth, fails first. It is a hypothesis.

## Next steps

1. Done (section 8): on a uniform 7 um mesh the crack reaches the bottom
   edge; initiation is late on that mesh.
1b. Done (section 9): on Ambati's mesh the crack does not start from the
   notch and the edges are damaged.
1c. History at quadrature points instead of DG0 cells (code change), on
   Ambati's mesh. Isolates the representation of `H`, the suspected reason why
   the notch does not initiate on coarse cells.
1d. Notch as in Ambati: an initial history field along the notch instead of
   `D = 1` (needs a small addition to the code). Isolates the notch
   representation, the suspected cause of the high peak.
2. Fix the cycling switch in the code (issue to open). Options: evaluate the
   hybrid mask once per step from the start-of-step displacement (closest
   precedent), or make it irreversible within the step (a cell switched off
   stays off until the next step, which ends the cycle by construction).
   Report the number of switched cells per step in the log. Check the damage
   cases of the suite against their golds.
3. Decide the energy split of the case: the star-convex split with
   `gamma_star: 5` has no precedent with the hybrid model and can make psi+
   negative in compression. On Ambati's mesh Amor behaves like star-convex
   (section 10) and the spectral split almost never converges.
4. With the loop converging, run a mesh fine at the notch tip (0.8 um) and
   uniform along the whole expected path (2 to 4 um, no coarsening towards the
   bottom) to the end of the loading and compare path, force-displacement and energy per unit length
   with Ambati Fig. 12 and 13.
5. Outcome: update the case (mesh, split), the gold with owner approval, and
   the README, or state in the README the measured differences and their
   cause.
