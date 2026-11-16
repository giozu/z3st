# The three 1D beam codes: implementation comparison

This note compares how the three 1D codes of the hexagonal duct are built: the
model, discretisation, section integrals, input, boundary conditions, output
and verification. On VP1 they give the same answer (12.2500 mm at the tip,
`beam_timoshenko_fx_NOTES.md` section 5). The differences are in what each
code can do beyond VP1, and how.

| Name | Path | Author | Size |
|---|---|---|---|
| `beam1d` | `U_hex_beam/beam1d/` (package, 9 modules) | Bianca | ~930 lines |
| `beam_fx.py` | `U_hex_beam/VP1/beam_fx.py` (reuses beam1d) | Bianca | 244 lines |
| `beam_timoshenko_fx.py` | `U_hex_shell/bowing/corner-corner/IAEA/` | ours | 331 lines |

Read from the code at `8fd37ec` (16 Nov 2026).

## 1. At a glance

| | beam1d | beam_fx | beam_timoshenko_fx |
|---|---|---|---|
| Beam theory | Euler–Bernoulli | Euler–Bernoulli | Timoshenko (shear GA), = EB in free bow |
| Unknowns | per plane (u, φ = u′) + rod u_z | per plane (u, M) + rod u_z | (u_z, w, θ) in one coupled system |
| Bending planes | x and y, separately | x and y, separately | one: the hot-to-cold direction |
| Element | cubic Hermite (C¹) + linear rod | Lagrange P3 (P2 option), mixed | Lagrange P2/P2/P2, mixed |
| Assembly | numpy, by hand, dense K | dolfinx/UFL, PETSc | dolfinx/UFL, PETSc |
| Linear solver | `np.linalg.solve` | MUMPS (saddle point) | PETSc LU via `LinearProblem` |
| Mesh | uniform nz + forced nodes at the stations | gmsh 1D `mesh.geo`, stations as physical points | `create_interval`, uniform n (multiple of 8 for the VP kinks) |
| E over the section | uniform only (`NotImplementedError`) | uniform only | E(T): S couples axial and bending |
| α | number or secant table α(T) | from beam1d | card hook (`thermal_strain`, `expansion_law`) |
| Section integrals | exact annulus, trapezoids, Gauss along the flat and through the wall | from beam1d | thin wall: mid-wall line, 2400 midpoints |
| Temperature input | z3st expression string or key-point table (Ma 2019, 12 points) | from beam1d | the case's Python temperature module (`T_vp1`, `T_vp2`) |
| Supports and loads | YAML: displacement, rotation, force, moment, pads, gravity | displacement, end rotations, forces, gravity | clamp at z = 0 (hard-coded), free top |
| Pads / contact | unilateral, active set | no | no |
| Output | `fields.npz`: u, u′, u″, M, N, fibre σ, reactions, pads | `fields_fx.npz`: u, u′, M, N | figure + `_summary.json`: w, u_z, θ vs closed form and 3D |
| Verification | `non-regression.py`, 31 checks at 1e-8 | vs beam1d | vs closed form (hex_post / vp2_post) and the 3D fields / gold |
| Cases | VP1 (own case folder) | VP1 | VP1, VP1 T-dependent, VP2 (reads the 3D case folders) |

## 2. Beam model and unknowns

**beam1d.** Euler–Bernoulli. Per bending plane a = x, y, the unknowns are the
centroid displacement u_a and its slope φ = u_a′ at each node. The plane is
loaded by the free curvature k_a = −κ*_a through the energy
½∫EI (u″ − k)² dz (`element.py`). The axial rod is solved separately with
½∫EA (u_z′ − ε₀*)² dz. Three independent linear problems (x, y, z), which is
valid because E is uniform (B1 = 0, Ma Eq. 17).

**beam_fx.** The same physics. The bending plane is in mixed form with the
moment as the second unknown, M = EI (u″ − k):

```
∫ M τ / EI + ∫ u′ τ′ + ∫ k τ = Σ_ends θ τ n          (constitutive, weak)
∫ M′ v′ + ∫ q v + Σ F_i v(z_i) = 0                   (equilibrium, weak)
```

A prescribed end rotation θ is natural (right-hand side). At a free end
M = 0 is imposed in the M space. The system is a saddle point, so it is
solved with MUMPS. The output contains M directly, without differentiating
u twice.

**ours.** Timoshenko with stretching. One coupled system in
(u_z, w, θ), with ε = u_z′ − c θ′ on the section:

```
∫ EA u′du′ − S(θ′du′ + u′dθ′) + EI θ′dθ′ + GA(w′ − θ)(dw′ − dθ) = ∫ N* du′ + M* dθ′
```

Only first derivatives, so C⁰ P2 suffices. In free bow the shear force is
zero, so w′ = θ and the result is EB. The single coupled system is what makes
E(T) possible: S ≠ 0 couples u_z and θ.

**Difference that matters.** Hers solves two bending planes. Ours solves only
the plane of `PHI_HOT` (hot-to-cold corner). A load with a component in the
other plane (flat-flat gradients, pads off-axis) needs a second plane in ours.
Ours has a shear stiffness. Hers is pure EB (no GA), which is exact for free
bowing but stiffer under concentrated pad forces.

## 3. Discretisation and mesh

| | beam1d | beam_fx | ours |
|---|---|---|---|
| Shape functions | Hermite N(ξ, L) with the physical L in the slope functions | Basix Lagrange P3 | Basix Lagrange P2 |
| Continuity of w | C¹ | C⁰ (u), C⁰ (M) | C⁰ (w, θ) |
| Non-uniform mesh | yes (slope functions carry L) | yes | yes (mixed); the Hermite option only uniform (dolfinx slope-dof bug) |
| Nodes at the stations | forced (`solver.mesh(H, nz, stations)`) | gmsh physical points, checked against geometry.yaml | by choosing n a multiple of 8 |
| Element integrals | Gauss in z, `n_gauss` (default 3) | UFL quadrature | UFL quadrature |

beam1d assembles a dense n × n matrix with numpy. For nz = 80 that is a
162 × 162 matrix, trivial. For thousands of elements it would be slow (dense
O(n²) memory, O(n³) solve), but such meshes are not needed in 1D.

## 4. Loads: from the section to the beam

All three integrate ε* = α(T − T_ref) over the section at every height.

| | beam1d / beam_fx | ours |
|---|---|---|
| Section | the exact wall (annulus between the outer and inner hexagons), trapezoids split at the quarter and mid points of each flat | the mid-wall line, dA = t dp |
| Quadrature | Gauss–Legendre along the flat and through the wall | 2400 midpoints along the line, none through the wall |
| Exact for | polynomials of low degree, T piecewise linear along the contour with kinks at the split points | linear T for κ (the thin-wall error cancels in M*/EI); O((t/a)²) error otherwise |
| Quantities | ε₀* = (1/A)∫ε* dA, κ* = (1/I)∫ε*(x, y) dA (vectors) | resultants N* = ∫Eε* dA, M* = −∫Eε* c dA, and EA, S, EI, GA |
| Along z | beam1d: evaluated directly at the element Gauss points. beam_fx: interpolated into DG_p (exact only with T kinks on nodes) | on 4001 levels, then `np.interp` to P1 coefficient functions at the mesh vertices (kinks must be nodes) |

**Consequences.**
- For constant E and linear T (VP1) both are exact.
- For a thick wall with a non-linear section T, beam1d's section is the more
  accurate.
- For E(T), ours interpolates the section sums with P1 between the vertices.
  That is the source of its 2e-5 mm error on VP1 T-dependent. Evaluating at the
  quadrature points, as beam1d does, would remove it.

## 5. Input and case structure

**beam1d.** A case folder laid out like a z3st case: `input.yaml`,
`geometry.yaml` (D, t, H, orientation, nz, stations),
`boundary_conditions.yaml` (thermal and mechanical blocks), material card,
`Allrun`, `Allclean`, `non-regression.py`. The temperature is a z3st-style
expression string (evaluated with `eval`) or a key-point CSV table (12 points
on the contour, Ma et al. 2019, `contour_linear` or `ma12` section models).
This is the closest to how z3st itself works, and the natural starting point
for a z3st beam regime.

**beam_fx.** Uses beam1d's case reader, plus a gmsh 1D mesh with the stations
as physical points. Boundary conditions are found by region name (gmsh tag),
not by height.

**ours.** No case folder of its own. It is a driver over the existing 3D case
folders: it reads `post.yaml`, the material card, the case's temperature
module and `geometry.yaml`. The temperature is the same Python function the
3D run uses as its eigenstrain, so 1D and 3D are guaranteed to see the same T.
The supports are hard-coded (clamped base, free top).

## 6. Boundary conditions and loads

| | beam1d | beam_fx | ours |
|---|---|---|---|
| Displacement at a station | yes | yes | base only |
| Rotation | at a station (Dirichlet on φ) | at the ends only (natural) | base only (Dirichlet on θ) |
| Point force | yes | yes (added to the vector) | no |
| Point moment | yes | no | no |
| Pads (unilateral, gap) | yes, active set | no | no |
| Gravity | yes (ρ A g) | yes | no |

For pads or springs, ours would add a term to the bilinear form (a node at the
pad height, k added to that dof, or a `dS` measure on that vertex). beam1d
already has the full contact logic.

## 7. Output and verification

**beam1d.** `fields.npz` with displacements, slopes, curvatures, M, N, fibre
stress σ(x, y, z) and its resultants, reactions and pad forces.
`non-regression.py` runs 31 checks at 1e-8: section properties against closed
forms and against the 3D mesh (`mesh3d.py`), free strains, stations, M, N, σ,
reactions, element matrices against Oechsner (K_e, rod), rigid modes, the
Kronecker property of the shape functions, a tip-load cantilever and a
3-element mesh. `fdm.py` gives a finite-difference cross-check (Ma Sect. 3.2).
`vp1_vs_z3st.py` compares with z3st 3D and splits the gap into thermal shear
and the rest.

**beam_fx.** `fields_fx.npz` (u, u′, M, N) and the comparison with beam1d.
`beam_fx_NOTES.txt` lists checks with loads (tip force, distributed load,
simply supported beam, imposed rotation), exact to the printed digits.

**ours.** One figure and one `_summary.json` per case: w, u_z and θ′ against
the closed form and the 3D fields and gold, with the station table. Three
cases (VP1, VP1 T-dependent, VP2). No element-level checks and no loaded
cases.

## 8. What each could take from the others

**Ours from beam1d.**
- The exact section with Gauss through the wall (removes the thin-wall error).
- Evaluating the loads at the quadrature points instead of interpolating
  (removes the E(T) interpolation error).
- The YAML case structure, stations and the boundary-condition reader.
- Pads with the active set.
- A second bending plane.

**beam1d / beam_fx from ours.**
- The S coupling (one coupled axial–bending system) for E(T) cards: beam1d
  currently stops on a non-uniform E.
- A shear stiffness GA, if pad forces become large or the duct short.
- The T-dependent and VP2 comparisons.

**Between beam_fx and ours (both dolfinx).**
- beam_fx's (u, M) mixed form gives M directly and imposes end rotations
  naturally.
- Ours' (w, θ) form gives the rotation directly, takes E(T) and has a shear
  stiffness.
- A dolfinx beam regime in z3st could start from either. The YAML layer
  already exists in beam1d.
