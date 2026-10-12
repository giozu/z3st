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
- With the notch as a band of fully damaged cells on Ambati's mesh (section
  12) the crack starts from the tip at 13.3 um (peak 0.567 kN, Ambati about
  0.49 kN at 12.5 um), follows Ambati's and Hirshikesh's path to the bottom
  edge, and the edges are no longer damaged. The line notch caused the
  missing initiation. Still different: the initial slope is 4 % above
  Ambati's and the peak 16 % above (section 13), and at 28.2 um, in a non-converged step, the
  bottom ligament breaks and the force drops to zero.
- The initial slope depends on how the notch is represented (section 13).
  An elastic run with the notch as a geometric slit, the exact sharp notch,
  gives 0.046 kN/um on Ambati's mesh and on a mesh twice as fine. Ambati's
  slope, about 0.043 kN/um, is 7 % below that of the sharp slit, so his
  initial history field softens the specimen more than a sharp notch would.
  The line notch is 18 % above the slit and the band notch 4 % below it.
  An initial history field narrower than one cell (Borden's width, lc / 2)
  is as stiff as the line notch on this mesh; from a half-width of 8 um the
  first row of nodes is broken and the slope is within 4 % of Ambati's
  (section 14).
  Plane strain, the boundary conditions and the material constants are the
  same as Ambati's.
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

### 11. History at quadrature points

Option `damage.history: quadrature` (added for this test): `H` at the points
of a degree-2 rule, 2 x 2 Gauss points on quadrilaterals, and the damage form
integrated with the same rule. Default `cell` unchanged (checked on
`regression/two_elliptical_cavities_2D` against its gold). On Ambati's mesh
with the settings of section 9, at step 620 (u about 14 um):

| Run | max D ahead of the tip (within 30 um) | max D at the top edge | max D at the bottom edge |
|---|---|---|---|
| H per cell (section 9) | 0.106 | 0.240 | 0.230 |
| H at quadrature points | 0.108 | 0.300 | 0.287 |

The force at 14 um is the same (0.74 and 0.75 kN). Storing the history at
the Gauss points does not make the notch initiate earlier; it increases the
damage at the edges. The hypothesis of section 10 is not supported.

Remaining candidate for the missing initiation: the notch. With `D = 1`
imposed only at the nodes of the notch line and cells of 7 um (above lc), the
damage drops from 1 to 0 within one cell on each side, so the cells across
the notch are only partly degraded and the notch still carries part of the
shear. On the stored mesh (0.8 um at the notch) the profile is resolved and
the faces are separated. Ambati models the notch with an initial history
field, which produces a damage band a few lc wide also on a coarse mesh.

### 12. Notch as a band of damaged cells

Same quadrilateral mesh and settings as section 9 (star-convex split,
`gamma_star: 5`, hybrid constraint, H per cell), one change: `D = 1` on the
notch line and on the two grid lines one cell above and below it (y = 0.5 mm
and 0.5 +- 1/144 mm, 0 <= x <= 0.5 mm), so the two rows of cells across the
notch are fully damaged and the faces are separated. This is a stand-in for
Ambati's initial history field, with no change to the code. 796 steps in
9 h 23 min, 224 of them not converged (69 up to the peak, the rest during
propagation).

| | Ambati | Line notch (section 9) | Band notch |
|---|---|---|---|
| force at 10 um | about 0.41 kN | 0.54 kN | 0.439 kN |
| first peak | about 0.49 kN at 12.5 um | 0.928 kN at 17.8 um | 0.567 kN at 13.3 um |
| crack from the notch tip | yes | no | yes, from 13.3 um |
| after the peak | about 0.35 kN, then rising | - | 0.433 kN at 20 um, rising to 0.498 kN at 28 um |
| crack tip at 21 um | - | - | about (0.71, 0.10) mm (Hirshikesh: (0.77, 0.08) mm) |
| crack at 25 um | - | - | reaches the bottom edge at x about 0.76 mm |
| force at 30 um | about 0.47 kN | 0.734 kN | 0.001 kN (bottom ligament broken at 28.2 um) |
| max D at the top / bottom edge | not shown | bands along both edges | 0.16 / 0.08 until the crack arrives |

With the faces separated the notch tip initiates a crack on Ambati's mesh, at
about Ambati's load, and the crack follows his path and Hirshikesh's to the
bottom edge, with the force rising again after the drop as in his Fig. 13.
The damage bands along the top and bottom edges of sections 9 to 11 are gone.
The line notch, not the split or the storage of `H`, caused the missing
initiation. Two differences remain. The initial slope is 4 % above Ambati's
(0.045 against about 0.043 kN/um) and the peak is about 16 % higher; section
13 traces the slope to the representation of the notch. At 28.2 um, in one step that did not
converge, the ligament between the crack and the lower-right corner breaks
along the bottom edge and the force drops to zero; the uniform triangles of
section 8 end in the same way. Ambati's force is still rising at 30 um.

### 13. Initial slope and the notch

Ambati states plane strain and fully integrated bilinear quadrilaterals
(p. 395); his Fig. 7b clamps the bottom edge and moves the top edge
horizontally with the vertical displacement fixed, the lateral edges free.
These are the conditions of this case, and his lambda = 121.15 and
mu = 80.77 kN/mm2 are E = 210 GPa and nu = 0.3, the constants of
`high_carbon_steel.yaml`. The force is per 1 mm of thickness in both.
Before the peak his curves overlap. Digitising the Miehe curve, the one
that can be separated from the grid in the image, gives F = 0.22 kN at 5 um
and 0.43 kN at 10 um, a slope of about 0.043 kN/um (a first reading by eye
gave 0.041).

Elastic runs (damage off, two steps of 1 um) with the notch as a geometric
slit, the two faces being separate lines that share only the tip, give the
slope of a sharp notch:

| Notch | Mesh | Slope at 1 um (kN/um) | against the slit |
|---|---|---|---|
| geometric slit, no damage | 144 x 144 quadrilaterals | 0.0470 | - |
| geometric slit, no damage | 288 x 288 quadrilaterals | 0.0465 | - |
| `D = 1` on a line (section 9) | 144 x 144 quadrilaterals | 0.0552 | +18 % |
| `D = 1` on a line (section 8) | uniform 7 um triangles | 0.0544 | +17 % |
| `D = 1` on a band (section 12) | 144 x 144 quadrilaterals | 0.0448 | -4 % |
| Ambati, initial history field | 20592 quadrilaterals | about 0.043 | about -7 % |

The slit converges with the mesh (1 % between the two meshes). The line
notch is stiffer than a sharp notch: the faces stay partly bonded, as found
in section 12. The band notch is slightly softer than the slit, as expected
from two rows of cells with `D = 1` and the damage profile around them.
Ambati's specimen is about 7 % softer than a sharp notch from the first
step, which points to his initial history field, a damage band several lc
wide along the notch, rather than to the elasticity. Z3ST's elastic
response is not in question; the initial slope measures the notch.

### 14. Notch as an initial history field

Option `damage.notch_history` (added for this test): `H0 = B (1 - d / w)`
within a distance `w` of the notch segment, with `B = 1000`, as in Borden et
al. (2012), whom Ambati follows. Z3ST normalises `H` by `Gc / (2 lc)`, so
`H0` needs no material constant. Borden's half-width is his `l0`, which is
`lc / 2` in the convention of Miehe and Ambati, so the default is
`w = lc / 2 = 2 um`. Ambati's mesh (section 9) with the notch line in the
mesh but no `D = 1` on it, history at quadrature points (section 11), four
steps of 1 um:

| w | mean D on the rows of nodes at 0 / 7 / 14 / 21 um from the notch | slope at 1 um (kN/um) |
|---|---|---|
| 2 um (lc / 2, Borden) | 1.00 / 0.30 / 0.04 / 0.005 | 0.0546 |
| 4 um (lc) | 1.00 / 0.30 / 0.04 / 0.005 | 0.0546 |
| 8 um | 1.00 / 0.99 / 0.13 / 0.016 | 0.0445 |
| 12 um | 0.99 / 1.00 / 0.80 / 0.10 | 0.0437 |
| Ambati | - | about 0.043 |

With `w` up to `lc` the history reaches only the Gauss points nearer the
notch line, the first row of nodes off the line keeps `D = 0.3`, and the
specimen is as stiff as with the line notch (0.055). Ambati, on cells of the
same size, gets about 0.043, so his initial field must reach the first row of
nodes; the paper does not give its width. From `w = 8 um`, about one cell,
the first row is broken and the slope is within 4 % of Ambati's, the same as
the band notch of section 12. A full run with `w = 8 um` is in progress.

### 15. Latching the hybrid mask within a step

The option of this section was removed from the code after section 18:
with `staggered: lagged` the mask is computed once per step and cannot
cycle, and with `coupled` the latch leaves 34 non-converged steps after the
peak. The sections below record what it did.

Option A of step 2, available as `damage.hybrid_latch: true` (default
false, see below). A history point is latched off
by the hybrid test (`psi- > psi+`) only when it switches off for the second
time in a step (off, on, off). From then on it stays off until the end of the
step, and every point is tested afresh at the next step (`_hybrid_offs` and
`_hybrid_latched` in `solver.py`, applied in `update_history`). A point that
switches once keeps the unlatched behaviour, so only a cycling point can be
latched and the period-2 cycle of section 6 ends after one period. The log
reports, at the end of each step, how many points are latched.

A first version latched a point at its first switch-off. It removed the
cycle in this case (section 16 uses it), but it also kept off the points
switched off by the early iterates of a step. In `pellet_quench_2D_xy` it
changed the iteration path (step 1: 9 iterations instead of 18) and left
step 6 non-converged, where the committed code converges. The second
switch-off rule replaced it.

First 46 steps (to 8.25 um) of the case of section 14:

| | Without latch (run of section 14) | Latch, second switch-off |
|---|---|---|
| non-converged steps | 11 | 0 |
| staggered iterations per step | 7 to 11, 200 when not converged | 7 to 12 |
| latched points per step | - | 0 to 2 |

Other cases with the second switch-off rule, against the committed code
(266964e):

| case | result |
|---|---|
| `pellet_quench_2D_xy`, 100 steps | same iterations at every step, same results to the last digit, gold PASS; one point latched, at step 97 |
| `two_elliptical_cavities_2D` | gold PASS |
| `sen_shear` as committed (line notch, triangle mesh), first 14 steps | 1 non-converged step against 6 |
| `sen_tension`, 701 steps | same iterations at every step, 3 non-converged steps as without, gold PASS (energies identical); 37 points latched, at one step |
| `plate_thermal_shock_2D`, 50 steps | gold FAIL: fracture energy -2.2 %, elastic energy +1.1 %, damaged node fraction +1.9 %, mean D -0.5 %; total energy, max D, crack counts and temperatures PASS. 0 non-converged steps with and without. About 9300 points latched at step 4, 50 to 350 per step after |

In `plate_thermal_shock_2D` the mask does cycle, and the latch moves the
converged solution to a different path: the same cracks, with about 2 % of
the energy moved from fracture to elastic. The latch is therefore an option
(`hybrid_latch`, default false): the default behaviour and the golds are
those of the committed code. With the option on, the 46-step run above is
reproduced exactly (same iterations, same latched points, same energy).

### 16. Full runs with the initial history field, with and without the latch

Case of section 14 (`w = 8 um`, history at quadrature points), 796 steps.
The latch here is the first version (first switch-off, section 15).

| | Without latch | With latch | Band notch (section 12) | Ambati |
|---|---|---|---|---|
| run time | 8 h 30 min | 5 h 3 min | 9 h 23 min | - |
| non-converged steps up to 12 um | 89 | 0 | - | - |
| non-converged steps after 12 um | 133 | 79 | - | - |
| force at 10 um | 0.436 kN | 0.436 kN | 0.439 kN | about 0.43 kN |
| peak | 0.573 kN at 14.1 um | 0.573 kN at 14.1 um | 0.567 kN at 13.3 um | about 0.49 kN at 12.5 um |
| force at 20 um | 0.377 kN | 0.380 kN | 0.433 kN | about 0.36 kN |
| crack tip at 20 um | (0.78, 0.15) mm | (0.78, 0.15) mm | (0.71, 0.10) mm at 21 um | - |
| crack reaches the bottom edge | about 24 um, x = 0.99 mm | about 24 um, x = 0.98 mm | about 25 um, x = 0.76 mm | near the lower-right corner |
| force at 30 um | 0.005 kN | 0.178 kN, rising from 0.168 at 28 um | 0.001 kN | about 0.47 kN, rising |

The latch removes every non-converged step before the peak and does not
change the force up to the peak. After the peak 79 steps still do not
converge. Their damage residual does not oscillate: it stalls at a level
between 1e-7 and 1e-3 and decreases slowly, while the latch keeps adding
points (up to about 700 per step) as the crack advances. This is the slow
convergence of the alternate minimisation during crack growth, a different
problem from the cycle of section 6.

With the latch the specimen keeps a residual force after the crack reaches
the bottom edge and the force rises again, as in Ambati, but at a lower
level (0.18 against about 0.47 kN at 30 um); without the latch it drops to
zero. The crack runs further right than with the band notch and than in
Ambati, and reaches the bottom edge about 6 um earlier. The peak is still
16 % above Ambati's.

### 17. The stored run predates the coupling of damage into the mechanics

The stored run and the gold of this case (8 August 2026) were produced with
z3st 0.2.0. Commit 8c24cf9 (15 September 2026) made the damage iterate
`self.D` itself, so that each mechanical solve within a step sees the current
damage. Before it, the mechanics saw `D^n` for the whole step: `u`, hence
`psi+` and `psi-`, did not change with the damage iterate, the hybrid mask
could not cycle and the loop converged in a few iterations. In effect this
is one staggered iteration per step (u with `D^n`, then d), the one-pass
scheme of Miehe et al. (2010). Ambati does not use it: he iterates u and d
within the step and stops on a normalised energy slope (his Sect. 3.4,
Eqs. 41-43, Tol about 10 deg), and his Figs. 10-11 show that one staggered
iteration per step gives inaccurate results, more so with larger increments.
The coupled scheme of HEAD is the one described in the docs
(`staggered_theory.rst`) and in the paper.

The case as committed, on HEAD (266964e), first steps:

| | Iterations per step | Energies at steps 1 to 8 |
|---|---|---|
| stored run (0.2.0) | 12 7 8 8 8 8 8 8 | - |
| HEAD | 12 7 8 X 14 25 26 X | - |
| HEAD with `D_new = self.D` reverted to a separate function | 12 7 8 8 8 8 8 8 8 8 8 | identical to the stored run in all printed digits |

Option `damage.staggered` (default `coupled`). With `lagged` the history
is computed once per step from the step-start displacement,
`H = max(H^n, psi+(u^n))`, and kept fixed while the loop solves d and then u
at that d: Box 1 and Eq. (54) of Miehe, Hofacker and Welschinger, CMAME 199
(2010) 2765-2778. It is not the 0.2.0 behaviour, which solved u with `D^n`
and then d with the history of that u (the two subproblems in the opposite
order). On this case, first 13 steps with `lagged`: iterations
12 2 7 8 8 8 8 8 8 8 8 8 8, no step non-converged; the fracture energy at
step n equals that of the 0.2.0 run at step n - 1 to all printed digits
(the crack lags by one increment), and the elastic energy is 0.04 to 0.4 %
higher. `coupled` keeps `regression/two_elliptical_cavities_2D` at PASS.

The difference between the stored run and every run of sections 3 to 16 is
the coupling, not the mesh, the code version otherwise, or round-off. All
runs from section 3 on use the coupled scheme; the stored run they are
compared with uses the 0.2.0 one. The cycle of section 6 and the slow
convergence after the peak belong to the coupled scheme.

### 18. Coupled with latch against lagged, and Ambati's stopping criterion

Case of section 16 (Ambati's mesh, history at quadrature points,
`notch_history` with `w = 8 um`, star-convex `gamma_star: 5`, damage
tolerance 1e-7), 700 steps (14.8 um), current latch (section 15). Total
energy printed at every staggered iteration in a scratch copy of the code.

| | `coupled` + `hybrid_latch` | `lagged` (Miehe) |
|---|---|---|
| run time | 3 h 51 min | 45 min |
| non-converged steps | 34, all after step 558 (12.6 um) | 0 |
| mean iterations per converged step | 18 up to step 550, 102 after | 6.5 |
| force at 5 / 10 / 12 um | 0.2215 / 0.4358 / 0.5171 kN | 0.2219 / 0.4358 / 0.5171 kN |
| force at 14 um | 0.570 kN | 0.583 kN |
| peak | 0.573 kN at 14.1 um | 0.587 kN at 14.6 um |

The two schemes give the same force up to 13 um. The lagged peak is 2.4 %
higher and 0.5 um later, as expected from the one-increment lag of the
history (Miehe 2010, Sect. 4.1.5). The coupled peak equals that of section 16
(first-version latch). Both are 17 to 20 % above Ambati's (about 0.49 kN at
12.5 um): the staggered scheme does not explain the gap.

Ambati's criterion (`gamma = atan(N (E^{N-1} - E^N) / (E^1 - E^N))`,
stop at `|gamma| <= 10 deg`) applied to the 34 non-converged steps: it stops
26 of them, at iteration 4 to 199, with the energy at that iteration within
5e-9 to 7e-5 (relative) of the energy at iteration 200; it never stops the
other 8. In these steps the total energy rises within the step (by 2e-6 to
2.5e-5 relative from the first iterate) and over the last 50 iterations
oscillates by 5e-8 to 1.2e-5: the cycle is negligible in energy, but the
criterion, which assumes a monotone decrease, does not detect its end.

### 19. Anderson acceleration on the committed case

`solver_settings.anderson_depth` (default 0, off) applies Anderson mixing (Storvik
et al., CMAME 381 (2021) 113822, Alg. 1 and 3) to the damage iterate: mixing
while the norm of the damage increment decreases, relaxation (`relax_D`)
once it rises, back to a restarted mixing after `anderson_switch` (default 5)
decreasing norms. Storvik switches on the mechanical residual of a
variational model; here the indicator is the damage increment, the quantity
that the convergence test uses.

Committed case (fine mesh, line notch, `coupled`, damage tolerance 1e-7),
split `amor`, `anderson_depth: 5`, all 796 steps, 8 h 7 min:

| | committed, star-convex + latch (step 2) | `amor` + Anderson |
|---|---|---|
| steps run | 524 (stopped at 13 um) | 796 (30 um) |
| non-converged steps | 163 | 41 |
| mean iterations per converged step | - | 7.9 to the peak, 18 after |
| first peak | 0.455 kN at 10.2 um | 0.460 kN at 10.5 um |
| force at 12 / 15 / 20 / 24.8 um | - | 0.417 / 0.393 / 0.402 / 0.363 kN |
| force at 30 um | - | 0.004 kN (bottom ligament broken at 25.0 um) |

Ambati (Fig. 13): first peak about 0.49 kN at 12.5 um, about 0.35 kN near
16 um, about 0.47 kN at 30 um. The split changes the peak by 1 %; the peak
is 6 % below Ambati's and 2 um earlier. After the peak the force stays at
0.36 to 0.42 kN without the minimum and the rise of Ambati's curve.

The non-converged steps are 11 and 12, three at the peak (267 to 269), and
36 spread over the propagation (294 to 769). Anderson does not remove them.
In these steps the damage increment falls to about 8e-7 and jumps to about
2.5e-6 when Gauss points switch side of the hybrid mask, so the residual
never reaches 1e-7; with star-convex the same steps (4, 8, 11 on a test run)
fail with and without a safeguarded restart of the mixing.

At step 770 (25.0 um) a band along the bottom edge, x from 0.93 to 0.99 mm,
breaks and the force falls from 0.363 to 0.003 kN, as at 28.2 um in the
stored run (section 12). The step converged in 132 iterations: the break is
not an artefact of a non-converged step (step 1f). At step 769 the main
crack has reached about (0.71, 0.10) mm.

## Next steps

1. Done (section 8): on a uniform 7 um mesh the crack reaches the bottom
   edge; initiation is late on that mesh.
1b. Done (section 9): on Ambati's mesh the crack does not start from the
   notch and the edges are damaged.
1c. Done (section 11): history at quadrature points does not change the
   initiation on Ambati's mesh.
1d. Done in part (sections 12 and 14): with the notch as a band of fully
   damaged cells the crack starts at about Ambati's load and follows his
   path. The initial history field is implemented (`notch_history`); with
   Borden's width it is too narrow for this mesh, with `w = 8 um` its slope
   matches the band. Full run with `w = 8 um` in progress.
1e. Done (section 13): plane strain, boundary conditions, constants and
   force are as in Ambati. The initial slope is set by the notch: a sharp
   slit gives 0.046 kN/um, Ambati about 0.043, the band notch and the
   initial history field with `w = 8 um` 0.0445 (section 14).
1f. Find why the bottom ligament breaks at 28.2 um. With `amor` and
   Anderson (section 19) it breaks at 25.0 um in a step that converged, so
   convergence does not remove it.
2. Done (sections 15 and 16). Option A (`damage.hybrid_latch`, off by
   default; a point latched off at its second switch-off in a step) removes
   the cycle before the peak, and the log reports the latched points per
   step. Suite checked (section 15): the latch changes
   `plate_thermal_shock_2D` only (about 2 % of the energy), hence the
   default. Still open: after the peak many steps converge slowly (79 with
   the first-version latch, section 16). The committed `sen_shear` with the
   latch was stopped at step 524 of 796 (13 um, 163 non-converged steps,
   peak 0.455 kN at 10.2 um), since the case is to be rebuilt (step 5).
   Issue to open on the cycle. The latch was later removed (section 15). On Ambati's mesh (section 18)
   `staggered: lagged` converges every step in a fifth of the time, with
   the same force up to 13 um and a 2.4 % higher peak. Anderson
   acceleration (section 19) cuts the non-converged steps of the committed
   case to 41 of 796 with `amor`, but does not remove the cycle.
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
