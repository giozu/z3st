# thermal_bending_plate_3D: free flat plate under a linear temperature

A flat steel plate, free of loads and restraints, carries a temperature that
varies linearly through its thickness. The hot face expands more than the cold
one and the plate bends into a shallow dome. Because nothing stops it bending,
it carries no stress. The case checks the temperature, the membrane strain, the
curvature, the whole displacement field and the vanishing stress against the
closed-form solution.

The derivation is in `thermal_stress_flat_plate.ipynb` (same folder). For a
free plate with an arbitrary T(x) it gives

    sigma_yy = sigma_zz = alpha E / (1 - nu) * [ T_mean
                          + 12 (x - L/2) / L^3 * int_0^L T (x - L/2) dx - T(x) ]

and the bracket vanishes when T(x) is linear. The companion case
`verification/thermal/thick_slab_non_adiabatic_3D` holds the plate flat
(`Clamp_x` on `xmin`), so it checks the membrane term alone; this case checks
the bending.

In the local suite (carries a gold). Not in `cases_ci.txt`.

## Geometry and mesh

x runs through the thickness, y and z lie in the plane of the plate.

| | value |
|---|---|
| thickness Lx | 0.2 m |
| half-spans Ly, Lz | 1.0 m (a quarter of a 2 m x 2 m plate) |
| mesh | structured hexahedra, 21 x 41 x 41 nodes, 32000 cells |
| element size | 10 mm through the thickness, 25 mm in the plane |
| elements | linear Lagrange for T and u |

Labels: `xmin` and `xmax` are the plate faces, `ymin` and `zmin` the symmetry
planes, `ymax` and `zmax` the outer edges.

## Material

`vessel_steel_0.yaml` (shared card): E = 177 GPa, nu = 0.30,
alpha = 1.7e-5 1/K, k = 48.1 W/m K, T_ref = 300 K, no gamma heating.

## Boundary conditions

Thermal:

- T = 600 K on `xmin`, T = 400 K on `xmax`;
- all other faces adiabatic (natural condition).

With constant k and no source the exact temperature is T(x) = 600 - 1000 x K.

Mechanical:

- `Clamp_y` on `ymin` and `Clamp_z` on `zmin`: symmetry planes. The exact
  solution has u_y = 0 on y = 0 and u_z = 0 on z = 0, so they add no spurious
  restraint. They block the translations in y and z and all three rotations.
- all other faces traction-free; no mechanical load.

The translation in x is left free. A Dirichlet condition on u_x over any face
or edge would conflict with the exact u_x, which varies everywhere, and Z3ST
has no point constraint. `remove_rigid_nullspace: true` finds the free mode and
attaches it to the matrix as its nullspace. The matrix is then singular, so the
mechanical solve uses CG with GAMG (`iterative_amg`) rather than LU. The
displacement is determined up to a constant in x, and `non-regression.py`
removes the mean shift before comparing.

The singular system is solvable because a thermal load does no work on a rigid
motion: f_i = int eps(phi_i) : D eps_th dV, and eps(r) = 0 for a rigid motion r.

## Solution procedure

One steady step. The staggered loop runs with relaxation 1: T does not depend
on u, so the second iteration only confirms convergence.

## Analytic solution

With g = (To - Ti) / Lx = -1000 K/m the total strain equals the thermal strain,
eps_ij = alpha (T - T_ref) delta_ij, and the displacement is, up to an
x-translation,

    u_x = alpha [ (Ti - T_ref) x + g x^2 / 2 - g (y^2 + z^2) / 2 ]
    u_y = alpha (T(x) - T_ref) y
    u_z = alpha (T(x) - T_ref) z

so that

- eps_0 = alpha (T_mean - T_ref) = 3.4e-3,
- kappa = alpha g = -0.017 1/m (radius about 59 m),
- sigma = 0 everywhere, edges included.

## Checks (tolerance 2e-2)

| metric | how | rel. error |
|---|---|---|
| `L2_error_T` | nodal T against 600 - 1000 x | 2e-15 |
| `eps_0` | u_y / Ly at x = Lx/2, y = Ly, z = 0 | 5e-9 |
| `kappa` | parabola fitted to u_x(y) on x = Lx/2, z = 0; kappa = -2 x (y^2 coefficient) | 3.2e-3 |
| `L2_error_u` | whole displacement field against the exact one, mean x-shift removed | 9.2e-4 |
| `rms_sigma_xx/yy/zz` | rms over all cells, divided by alpha E / (1 - nu) abs(To - Ti) | 9.4e-4 |

The stress scale in the last row is the stress the same gradient would cause if
the plate were held flat (about 430 MPa at the faces). Without the bending term
the analytic stress would differ from zero by that amount.

## Mesh

The curvature carries the largest error. It comes from the in-plane element
size: linear hexahedra cannot represent the parabola u_x ~ y^2 within an
element.

| nodes (x, y, z) | element (mm) | kappa error | L2 u error | run time |
|---|---|---|---|---|
| 21, 21, 21 | 10 x 50 x 50 | 1.56e-2 | 4.6e-3 | 34 s |
| 11, 21, 21 | 20 x 50 x 50 | 1.26e-2 | 3.8e-3 | 16 s |
| 21, 41, 41 (committed) | 10 x 25 x 25 | 3.20e-3 | 9.2e-4 | 40 s |

Halving the in-plane size divides the curvature error by about 5; the
thickness division matters little.

## Figures (`plots.py`, run by `Allrun`)

- `mesh.png`: the undeformed mesh.
- `deformed.png`: the plate with the displacement magnified 10 times, coloured
  by u_x, the undeformed outline behind it.
- `profiles.png`: T, total strain and stress (per cell) across the thickness at
  the plate centre, against the analytic lines. The stress panel also shows the
  held-flat stress for scale.

## Run

    conda activate z3st
    ./Allrun      # gmsh -> z3st -> non-regression -> plots
    ./Allclean
