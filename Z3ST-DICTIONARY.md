# Z3ST dictionary

Commands and keys you'll use most when you write a z3st case, or a `.py` script that
reads a case, drives it, or extends it. Names were checked against z3st 0.3.2 and
dolfinx 0.10/0.11.

- **Part A: FEniCSx.** The libraries z3st is built on (dolfinx, UFL, PETSc, MPI, gmsh).
  You need them for custom material callables, `diagnostics.py`, and standalone scripts.
- **Part B: z3st native.** The CLI, the case files and their YAML keys, the `Spine` object,
  and the post-processing utilities.
- **Part C: shortcuts.** VS Code, Jupyter notebooks and the bash terminal (WSL).

---

# Part A: FEniCSx

## A.1 Imports

```python
import numpy as np
import ufl
import dolfinx
from dolfinx import fem, mesh, io, plot
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from mpi4py import MPI
from petsc4py import PETSc
```

## A.2 Mesh

| Command | What it does |
|---|---|
| `mesh.create_interval(MPI.COMM_WORLD, n, [x0, x1])` | 1D mesh |
| `mesh.create_rectangle(comm, [[x0,y0],[x1,y1]], [nx,ny], cell_type=mesh.CellType.triangle)` | 2D rectangle (`quadrilateral` also works) |
| `mesh.create_box(comm, [[x0,y0,z0],[x1,y1,z1]], [nx,ny,nz], mesh.CellType.hexahedron)` | 3D box |
| `io.gmsh.read_from_msh("mesh.msh", comm, rank=0, gdim=3)` | Reads a Gmsh file. Returns `.mesh`, `.cell_tags`, `.facet_tags`. z3st calls this for you. |
| `msh.topology.dim` / `msh.geometry.dim` | Topological / geometric dimension |
| `fdim = tdim - 1` | Facet dimension |
| `msh.topology.create_connectivity(fdim, tdim)` | Required before facet→cell lookups |
| `msh.geometry.x` | Node coordinates `(n_nodes, 3)` |
| `msh.topology.index_map(tdim).size_local` | Number of cells owned by this rank |
| `mesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], 0.0))` | Boundary facets matching a geometric marker |
| `mesh.locate_entities(msh, tdim, marker)` | Cells matching a marker |
| `mesh.meshtags(msh, fdim, indices, values)` | Builds a tag object (indices must be sorted) |
| `facet_tags.find(tag_id)` | Facet indices that carry a given physical tag |
| `cell_tags.values`, `cell_tags.indices` | Raw tag arrays |

## A.3 Function spaces and functions

| Command | What it does |
|---|---|
| `fem.functionspace(msh, ("Lagrange", 1))` | Scalar P1 space (temperature) |
| `fem.functionspace(msh, ("Lagrange", 1, (msh.geometry.dim,)))` | Vector P1 space (displacement) |
| `fem.functionspace(msh, ("DG", 0))` | Cell-wise constant space (stress, history variables) |
| `fem.functionspace(msh, ("DG", 0, (3, 3)))` | Cell-wise tensor space |
| `V.sub(i)` / `V.sub(i).collapse()` | Component subspace (`collapse` returns `(V_i, dof_map)`) |
| `fem.Function(V, name="Temperature")` | Discrete field. `name` is what ParaView shows. |
| `fem.Constant(msh, PETSc.ScalarType(300.0))` | Constant you can update later through `.value` without recompiling the form |
| `u.x.array` | Local dof values as a numpy view (owned dofs plus ghosts) |
| `u.x.scatter_forward()` | Syncs ghost values after you edit `u.x.array` by hand |
| `u.x.array[:] = v.x.array` | Copy values (same space) |
| `u.interpolate(lambda x: 300 + 10*x[0])` | Fills from a Python function of `x` (shape `(3, N)`) |
| `u.interpolate(fem.Expression(expr, V.element.interpolation_points))` | Fills from a UFL expression (for example, project σ onto DG0) |
| `V.dofmap.index_map.size_local * V.dofmap.index_map_bs` | Number of owned dofs (use it to slice `x.array` before MPI reductions) |
| `V.tabulate_dof_coordinates()` | Coordinates of each dof |

## A.4 Boundary conditions (raw dolfinx)

| Command | What it does |
|---|---|
| `dofs = fem.locate_dofs_topological(V, fdim, facets)` | Dofs that sit on the given facets |
| `dofs = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[1], 0))` | Dofs located by coordinate |
| `fem.dirichletbc(PETSc.ScalarType(600.0), dofs, V)` | Dirichlet BC with a scalar value |
| `fem.dirichletbc(u_bc_function, dofs)` | Dirichlet BC from a Function |
| `fem.dirichletbc(0.0, fem.locate_dofs_topological(V.sub(2), fdim, f), V.sub(2))` | Fixes one component only (a z3st `Clamp_z`) |
| Neumann / Robin | No object needed: add `ds(tag)` terms to the weak form |

## A.5 UFL (weak forms)

| Symbol | Meaning |
|---|---|
| `u, v = ufl.TrialFunction(V), ufl.TestFunction(V)` | Unknown and test function |
| `ufl.grad`, `ufl.div`, `ufl.nabla_grad` | Differential operators |
| `ufl.inner(a, b)`, `ufl.dot(a, b)` | Full contraction / single contraction |
| `ufl.sym(ufl.grad(u))` | Small-strain tensor ε(u) |
| `ufl.tr(A)`, `ufl.dev(A)`, `ufl.det(A)`, `ufl.Identity(d)` | Tensor algebra |
| `ufl.SpatialCoordinate(msh)` | `x` inside forms (`x[0]` is r in axisymmetric problems) |
| `ufl.FacetNormal(msh)` | Outward normal `n` |
| `ufl.CellDiameter(msh)` | Element size `h` |
| `ufl.exp`, `ufl.ln`, `ufl.sqrt`, `ufl.sin`, `ufl.cos`, `ufl.tanh` | Math (use these, not `np.*`, inside forms) |
| `ufl.conditional(ufl.gt(T, 1000), a, b)` | if/else (`lt`, `ge`, `le`, `eq`) |
| `ufl.max_value(a, b)`, `ufl.min_value(a, b)` | Pointwise max / min |
| `ufl.as_vector([...])`, `ufl.as_matrix([[...]])`, `ufl.as_tensor` | Builds a vector or tensor |
| `ufl.derivative(F, u, du)` | Gâteaux derivative (Jacobian for Newton) |
| `ufl.dx`, `ufl.ds`, `ufl.dS` | Cell, exterior-facet, interior-facet measures |
| `dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)` → `dx(1)` | Integral over material tag 1 |
| `ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags)` → `ds(12)` | Integral over facet tag 12 |
| `ufl.Measure("dx", ..., metadata={"quadrature_degree": 4})` | Sets the quadrature degree |

Typical steady heat equation with k(T), a flux and a convective BC:

```python
a = ufl.inner(k * ufl.grad(u), ufl.grad(v)) * dx + h * u * v * ds(13)
L = q3 * v * dx - q_flux * v * ds(12) + h * T_ext * v * ds(13)
```

Linear thermo-elasticity:

```python
eps   = lambda w: ufl.sym(ufl.grad(w))
sigma = lambda w: 2*mu*eps(w) + lmbda*ufl.tr(eps(w))*ufl.Identity(d)
a = ufl.inner(sigma(u), eps(v)) * dx
L = ufl.inner((3*lmbda + 2*mu) * alpha * (T - T_ref) * ufl.Identity(d), eps(v)) * dx
```

## A.6 Solving

| Command | What it does |
|---|---|
| `LinearProblem(a, L, bcs=[bc], petsc_options={...}, petsc_options_prefix="th_")` | Linear solve (the prefix is mandatory in dolfinx ≥ 0.10) |
| `uh = problem.solve()` | Returns the solution Function |
| `NonlinearProblem(F, u, bcs=[bc], petsc_options={...}, petsc_options_prefix="mech_")` then `problem.solve()` | Newton via SNES. `F` is the residual in `u` (no TrialFunction). |
| `{"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}` | Direct solver (z3st `direct_mumps`) |
| `{"ksp_type": "cg", "pc_type": "hypre", "pc_hypre_type": "boomeramg"}` | Iterative solver (z3st `iterative_hypre`) |
| `{"ksp_type": "gmres", "pc_type": "gamg"}` | Iterative solver (z3st `iterative_amg`) |
| `{"snes_type": "newtonls", "snes_linesearch_type": "bt", "snes_rtol": 1e-8}` | Newton options |
| `{"ksp_monitor": None, "snes_monitor": None}` | Prints residuals |

## A.7 Integrals and MPI-safe reductions

| Command | What it does |
|---|---|
| `val = fem.assemble_scalar(fem.form(T * dx))` | Partial integral **on this rank only** |
| `msh.comm.allreduce(val, op=MPI.SUM)` | Global value (always do this, even in serial) |
| `msh.comm.allreduce(local_max, op=MPI.MAX)` | Global max (for example T_max) |
| `msh.comm.rank == 0` | Guard for prints and file writes |
| `mpirun -n 4 python -m z3st` | Parallel run |

## A.8 I/O and visualisation

| Command | What it does |
|---|---|
| `with io.XDMFFile(comm, "out.xdmf", "w") as f: f.write_mesh(msh); f.write_function(u, t)` | XDMF time series (P1 only) |
| `io.VTXWriter(comm, "out.bp", [u], engine="BP4")` | ADIOS2 output (any order) |
| `io.XDMFFile(comm, "mesh.xdmf", "r").read_mesh()` | Reads a mesh |
| `topology, cells, geom = plot.vtk_mesh(V)` → `pyvista.UnstructuredGrid(...)` | Quick PyVista plot |

## A.9 gmsh (`.geo`)

| Command | What it does |
|---|---|
| `gmsh mesh.geo -2` / `-3` | Meshes in 2D / 3D and writes `mesh.msh` |
| `gmsh -setnumber D 0.25 -setstring orientation vertical mesh.geo -3` | Overrides parameters from the command line |
| `If (!Exists(D)) D = 0.200; EndIf` | Default that `-setnumber` can override |
| `Physical Volume("steel", 1) = {...};` | Cell tag (must match `labels:` in `geometry.yaml`) |
| `Physical Surface("inner", 12) = {...};` | Facet tag (3D). In 2D: `Physical Curve` for facets, `Physical Surface` for cells. |
| `Transfinite Curve {..} = n+1;` / `Transfinite Surface` / `Recombine Surface` | Structured quad or hex mesh |
| `Extrude {0,0,H} { Surface{..}; Layers{nz}; Recombine; }` | Extrudes 2D to 3D |
| `Mesh.MshFileVersion = 4.1;` | Format dolfinx reads |

---

# Part B: z3st native

## B.1 Command line

| Command | What it does |
|---|---|
| `./Allrun` | gmsh → `python -m z3st` → `non-regression.py` → convergence plot |
| `./Allclean` | Removes generated output |
| `python -m z3st > log_z3st.md` | Runs the case in the current folder (reads `input.yaml`). Redirected output becomes Markdown. |
| `python -m z3st --mesh_plot` | Previews facet tags with PyVista before solving |
| `python -m z3st --debug` | Verbose output and HeatFlux diagnostics |
| `Z3ST_PLAIN_LOG=1 python -m z3st > log.txt` | Plain log, no Markdown rewrite |
| `python non-regression.py` | Analytic check and gold check → `output/non-regression.json` |
| `python -m z3st.utils.plot_convergence log_z3st.md` | Staggered-residual plot (`convergence.png`) |
| `gmsh $(python3 -m z3st.utils.geo_args) mesh.geo -3` | Passes the scalars from `geometry.yaml` to gmsh as `-setnumber` / `-setstring` |
| `z3st/cases/non-regression_local.sh` (`--list`) | Local suite: every case that has an `Allrun` and a gold file (`sandbox/` is skipped) |

Minimal `Allrun`:

```bash
#!/bin/bash
DIM=3                      # gmsh dimension
# POST="plots.py"          # optional extra scripts
. "$(python3 -c 'import z3st,os;print(os.path.dirname(z3st.__file__))')/utils/allrun.sh"
```

## B.2 Case folder

```
my_case/
├── Allrun / Allclean
├── input.yaml                  # what to solve and how
├── geometry.yaml               # dimensions + label ↔ tag map
├── boundary_conditions.yaml    # thermal / mechanical / damage BCs
├── mesh.geo  (→ mesh.msh)
├── my_material.yaml            # optional case-local material card
├── k_mymat.py                  # optional case-local callable (k(T), ...)
├── diagnostics.py              # optional per_step(problem, step, t) hook
├── non-regression.py
└── output/  fields.xdmf/.h5  or  fields_NNNN.vtu,  non-regression.json,  non-regression_gold.json
```

## B.3 `input.yaml`

| Key | Values / meaning |
|---|---|
| `mesh_path`, `geometry_path`, `boundary_conditions_path` | Paths relative to the case folder |
| `materials: {name: path.yaml}` | `name` must match the cell label in `geometry.yaml` and the keys in the BC file |
| `regime` | `1d` \| `2d` (plane strain) \| `3d` \| `axisymmetric` (x = r, y = z) |
| `models.thermal` / `mechanical` / `damage` / `plasticity` / `cluster` / `porosity` | `true` / `false` |
| `models.gap_conductance` | `{type: Fixed, value: 5000}` or `{type: Gas, ...}` |
| `models.contact` | `{surface_a, surface_b, penalty_stiffness, initial_gap}` (PCMI) |
| `solver_settings.max_iters` | Maximum staggered iterations per step |
| `solver_settings.relax_T` / `relax_u` / `relax_D` | Under-relaxation factors (1.0 means none) |
| `solver_settings.relax_adaptive`, `relax_growth`, `relax_shrink`, `relax_min`, `relax_max` | Adaptive relaxation |
| `solver_settings.relax_aitken` | Aitken Δ² relaxation on u |
| `thermal.analysis` | `stationary` \| `transient` |
| `<physics>.solver` | `linear` \| `nonlinear` (SNES) |
| `<physics>.linear_solver` | `direct_mumps` \| `iterative_amg` \| `iterative_hypre` |
| `<physics>.rtol`, `stag_tol`, `convergence` | KSP tolerance, staggered tolerance, `rel_norm` \| `norm` |
| `mechanical.order` | FE order of u (1 or 2) |
| `mechanical.remove_rigid_nullspace` | `true` when only thermal load acts and no Dirichlet BC on u (free body) |
| `damage.type`, `lc`, `split`, `hybrid_constraint` | `AT1`\|`AT2`, length scale, `amor`\|`miehe`\|`star_convex` |
| `time: [t0, t1, ...]`, `lhr: [q0, q1, ...]` | Piecewise-linear power history (s, W/m) |
| `n_steps` | Total number of steps (int) or a per-segment list |
| `output.format` | `xdmf` (single time series) \| `vtu` (one file per step) |

Can be edited **while the case runs** (applied at the next step): `*.stag_tol`, `*.rtol`,
`solver_settings.relax_*`, `max_iters`, `damage.hybrid_constraint`/`gamma_star`.
Changes to anything else need a restart.

## B.4 `geometry.yaml`

```yaml
name: hex_shell_3d
D: 0.200          # top-level scalars → geometry parameters (and geo_args → gmsh)
Lz: 1.000         # z3st reads Lx / Ly / Lz / Ri / Ro when they are present
labels:           # name → Physical tag in mesh.geo
  steel: 1        # cell tag  = material name in input.yaml
  inner: 12       # facet tag = region in boundary_conditions.yaml
```

## B.5 `boundary_conditions.yaml`

Structure: `physics → material → list of BCs`. Every `region` must be a key in `labels`.

| Physics | `type` | Keys |
|---|---|---|
| thermal | `Dirichlet` | `temperature` (K, scalar or list of length `n_steps`) |
| thermal | `Neumann` | `flux` (W/m², **positive = out of the body**, so heat entering is negative) |
| thermal | `Robin` | `h_conv` (W/m²K) + `T_ext` (K) for convection, **or** `pair: <region>` for gap coupling |
| mechanical | `Clamp_x` / `Clamp_y` / `Clamp_z` | Fixes that component to 0 (no `Clamp_z` in 2D or axisymmetric: use `Clamp_y`) |
| mechanical | `Dirichlet_x/y/z` | `displacement` (or `value`), m, scalar or list |
| mechanical | `Dirichlet` | `displacement: [ux, uy(, uz)]` |
| mechanical | `Neumann` | `traction` (Pa, along the normal; negative = pressure pushing in) |
| mechanical | `Slip_x/y/z` | Symmetry plane |
| damage | `Dirichlet` | Fixed `D` |

```yaml
thermal:
  steel:
    - { type: Neumann,   region: inner, flux: -5000.0 }
    - { type: Dirichlet, region: outer, temperature: 600.0 }
mechanical:
  steel:
    - { type: Clamp_z, region: bottom }
```

## B.6 Material card (`*.yaml`)

| Key | Units / meaning |
|---|---|
| `E`, `nu` | Pa, – (z3st derives `lmbda`, `G`, `bulk_modulus`) |
| `alpha`, `T_ref`, `T_initial` | 1/K, K (stress-free reference temperature), K (initial condition) |
| `k` | W/mK, **or** `"module.func"` → callable `k(T)` (a case-local `.py` file works) |
| `cp`, `rho` | J/kgK, kg/m³ |
| `constitutive` | `lame` \| `hyperelastic` \| `plasticity` \| `custom` (+ `stress_function: "pkg.mod.func"`) |
| `yield_strength`, `hardening_modulus` | J2 plasticity |
| `sigma_c` or `Gc` | Phase field (z3st computes one from the other using `lc`) |
| `fissile`, `heavy_metal_fraction`, `radial_profile`, `axial_profile` | Fuel power source |
| `swelling` / `eigenstrain` (+ `swelling_rate`) | Constant or burnup-driven swelling |
| `creep: norton`, `creep_A0`, `creep_n`, `creep_Q`, `creep_irr_B`, `fast_flux` | Creep |
| `gamma_heating`, `mu_gamma`, `gamma_inner_radius` | γ heating |
| `cracking: isotropic` | Fuel cracking (Barani) |

Case-local callable (as in `U_hex_shell/k_1515.py`). Write it so it works on both UFL
expressions and numpy arrays, so `non-regression.py` can reuse it:

```python
def k(T):
    return 13.95 + 0.01163 * (T - 273.15)   # only + - * / and ufl-compatible math
```

## B.7 Python API: driving z3st from a script

The same sequence `__main__.py` runs (needs dolfinx; run it from the case folder):

```python
from z3st.core.spine import Spine
from z3st.utils.utils_load import load, generate_power_history

inp  = load("input.yaml")
geom = load(inp["geometry_path"])
mats = {name: load(path) for name, path in inp["materials"].items()}

problem = Spine(input_file=inp, mesh_file=inp["mesh_path"], geometry=geom)
problem.load_materials(**mats)
times, lhrs, _ = generate_power_history(inp["time"], inp["lhr"], n_steps=inp["n_steps"] - 1, filename=None)
problem.n_steps = len(times)
problem.parameters(lhr=lhrs[0]); problem.initialize_fields(); problem.set_boundary_conditions()

for step in range(1, len(times)):
    dt = times[step] - times[step - 1]
    problem.current_step = step
    problem.parameters(lhr=lhrs[step]); problem.set_power(); problem.update_state(dt)
    problem.solve(max_iters=inp["solver_settings"]["max_iters"], dt=dt)
    problem.get_results()
```

| Attribute / method | What it is |
|---|---|
| `problem.mesh`, `cell_tags`, `facet_tags`, `label_map` | dolfinx mesh, tags, `{name: id}` from `geometry.yaml` |
| `problem.mgr` | `MeshManager`: `locate_domain_dofs(tag, V)`, `locate_facets_dofs(tag, V)`, `summary()` |
| `problem.V_t`, `V_m`, `V_d`, `Q` | Spaces: T (P1), u (vector P-order), D (P1), DG0 |
| `problem.T`, `u`, `D`, `burnup`, `q_third` | Solution and state `Function`s (use `.x.array` for values) |
| `problem.materials[name]` | Resolved material dict (includes `lmbda`, `G`, callables) |
| `problem.dx_tags[tag]`, `ds_tags[tag]` | Per-tag UFL measures (built inside `solve`) |
| `problem.weight` | Integration weight: `2πr` in axisymmetric, 1 otherwise. Multiply integrands by it. |
| `problem.on` | `{"thermal": True, ...}` model switches |
| `problem.regime`, `tdim`, `fdim` | Regime and dimensions |
| `problem.current_step`, `n_steps` | Step counters |
| `problem.get_results()` | Builds the UFL `stress`, `strain`, `stress_th`, energy dicts per material |
| `problem.solve(max_iters, dt)` | One staggered step. Returns `converged`. |
| `problem.heat_flux(problem.T)` | Heat-flux diagnostic |
| `problem.compute_energy_balance(problem.u)` | `(E_el, E_frac)` with damage enabled |
| `problem.snapshot_state()` / `restore_state(snap)` | Save and roll back (used by time adaptivity) |

## B.8 `diagnostics.py`: a case-local hook

If this file sits in the case folder, `python -m z3st` loads it and calls `per_step` after every step:

```python
import numpy as np
from mpi4py import MPI

def per_step(problem, step, t):
    comm = problem.mesh.comm
    imap = problem.V_t.dofmap.index_map
    T = problem.T.x.array[: imap.size_local]                 # owned dofs only
    T_max = comm.allreduce(T.max() if T.size else -np.inf, op=MPI.MAX)
    if comm.rank == 0:
        with open("output/history.csv", "a") as f:
            f.write(f"{step},{t},{T_max}\n")
```

Examples: `cases/regression/pwr_rod_2D/diagnostics.py`, `cases/verification/fuel/creep_shrink_fit_2D/diagnostics.py`.

## B.9 Reading results (no dolfinx required)

Output fields: `Temperature`, `Displacement` (nodal); `Stress`, `Strain`, `HeatFlux`,
`StrainEnergyDensity`, `VonMises`, `Hydrostatic`, `Damage`, `ContactPressure` (DG0 = per cell);
`Burnup`, `Porosity`, `ClusterDensity`.

| Function (`z3st.utils.…`) | What it does |
|---|---|
| `utils_extract_xdmf.list_fields_xdmf(xdmf)` | Field names in `fields.h5` |
| `utils_extract_xdmf.list_steps_xdmf(xdmf, field)` | Saved steps |
| `utils_extract_xdmf.extract_field_xdmf(xdmf, field, step_index=-1)` | `x, y, z, data`. Cell fields come back at the cell centres. For tensors, `data.reshape(-1, 3, 3)`. |
| `utils_extract_vtu.list_fields(vtu)` | Same for VTU |
| `utils_extract_vtu.extract_field(vtu, field)` | `x, y, z, data` |
| `utils_extract_vtu.extract_temperature` / `extract_displacement` / `extract_stress(vtu, "xx")` | Shortcuts |
| `utils_extract_vtu.extract_principal_stresses(...)`, `extract_cylindrical_field(...)` | Principal stresses, r-θ-z components |
| `utils_extract_vtu.average_section(...)`, `average_section_radial(...)` | Section averages |
| `utils_plot.plot_field_along_r_xyz(...)`, `plotter_sigma_temperature_cylinder/slab(...)` | Standard plots |
| `utils_load.load(path)` | YAML → dict |

## B.10 `non-regression.py` building blocks (`z3st.utils.non_regression`)

| Function | What it does |
|---|---|
| `case_paths(__file__)` | `CASE_DIR, VTU_FILE, OUT_JSON` |
| `load_yaml(CASE_DIR, "input.yaml")` | Parses one of the case's YAML files |
| `load_case(CASE_DIR, material=None)` | `(geometry, input, material)` |
| `line(vtu, "Temperature", y=0.5, z=0.0, tol=1e-5)` | Profile along the free coordinate, sorted |
| `metric(num, ref)` | Comparison against an analytic value (abs + rel error) |
| `error_metric(err)` | The value is already an error (L2 norm, non-uniformity) |
| `tracked(val)` | No reference: never fails the analytic check, still compared with the gold |
| `finish(errors, TOL, OUT_JSON, CASE_DIR)` | Analytic pass/fail and gold regression → `non-regression.json` |

Skeleton:

```python
import os
from z3st.utils.non_regression import finish, load_yaml, metric
from z3st.utils.utils_extract_xdmf import extract_field_xdmf

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
XDMF     = os.path.join(CASE_DIR, "output", "fields.xdmf")
OUT_JSON = os.path.join(CASE_DIR, "output", "non-regression.json")

x, y, z, T = extract_field_xdmf(XDMF, "Temperature")
errors = {"T_max_K": metric(T.max(), 612.3)}
finish(errors, 1e-2, OUT_JSON, CASE_DIR)
```

To bless a gold: once the case is right, copy `output/non-regression.json` to
`output/non-regression_gold.json`. From then on the case is part of the local suite.

---

# Part C: Shortcuts

VS Code on Windows connected to WSL. On a Mac, replace `Ctrl` with `Cmd`.

## C.1 Editing

| Shortcut | What it does |
|---|---|
| `Ctrl+/` | Comments or uncomments the selected lines (`#` in Python/YAML, `//` in C-like files) |
| `Shift+Alt+A` | Block comment around the selection (`"""…"""` in Python). YAML has no block comment: use `Ctrl+/`. |
| `Ctrl+]` / `Ctrl+[` | Indents / outdents the selected lines |
| `Alt+↑` / `Alt+↓` | Moves the line up / down |
| `Shift+Alt+↓` | Duplicates the line |
| `Ctrl+Shift+K` | Deletes the line |
| `Ctrl+D` | Selects the next occurrence of the word (repeat to add more cursors) |
| `Ctrl+Shift+L` | Selects all occurrences of the word |
| `Alt+Click` | Adds a cursor |
| `Shift+Alt+drag` | Column (box) selection, useful for YAML tables and BC lists |
| `F2` | Renames a symbol everywhere (Python) |
| `Ctrl+Space` | Triggers autocompletion |
| `Shift+Alt+F` | Formats the document |
| `Ctrl+Shift+[` / `Ctrl+Shift+]` | Folds / unfolds the current block |
| `Ctrl+K Ctrl+0` / `Ctrl+K Ctrl+J` | Folds / unfolds everything |
| `Ctrl+Z` / `Ctrl+Shift+Z` | Undo / redo |

`.geo` files have no comment support until you install a Gmsh syntax extension.
Until then, type `//` by hand, or use `Shift+Alt+drag` to put it on several lines at once.

## C.2 Navigating

| Shortcut | What it does |
|---|---|
| `Ctrl+P` | Opens a file by name (type `input.yaml`, `non-regression`, ...) |
| `Ctrl+Shift+P` | Command palette: every VS Code command, searchable |
| `Ctrl+G` | Goes to a line number |
| `Ctrl+Shift+O` | Jumps to a function or class in the current file |
| `F12` / `Ctrl+Click` | Goes to definition (for example, from `extract_field_xdmf` into z3st's code) |
| `Alt+F12` | Shows the definition in a pop-up, without leaving the file |
| `Shift+F12` | Finds all references |
| `Alt+←` / `Alt+→` | Goes back / forward after a jump |
| `Ctrl+F` / `Ctrl+H` | Find / replace in the file |
| `Ctrl+Shift+F` / `Ctrl+Shift+H` | Find / replace across the project (for example, every case that uses `Clamp_z`) |
| `Ctrl+Tab` | Switches between open files |
| `Ctrl+W` | Closes the tab |
| `Ctrl+\` | Splits the editor (for example, `input.yaml` next to `boundary_conditions.yaml`) |
| `Ctrl+B` | Shows / hides the sidebar |
| `Ctrl+Shift+E` | File explorer |
| `Ctrl+Shift+G` | Source control (git changes) |

## C.3 Markdown

| Shortcut | What it does |
|---|---|
| `Ctrl+Shift+V` | Opens the rendered preview in a new tab |
| `Ctrl+K V` | Opens the preview **side by side** (scrolls with the editor) |
| Right-click the file in the explorer → *Open Preview* | Same, from the explorer |

`log_z3st.md` is written in Markdown, so `Ctrl+K V` gives you a readable run log, with one heading per step and iteration.

## C.4 Terminal

| Shortcut / action | What it does |
|---|---|
| ``Ctrl+` `` | Shows / hides the integrated terminal |
| ``Ctrl+Shift+` `` | Opens a new terminal |
| Right-click a folder in the explorer → *Open in Integrated Terminal* | Terminal that already sits **in that folder** |
| `Ctrl+Shift+5` | Splits the terminal |
| `Ctrl+↑` / `Ctrl+↓` | Jumps to the previous / next command in the terminal output |
| `Ctrl+Shift+C` | Opens an external terminal in the workspace folder |

In bash:

| Command / key | What it does |
|---|---|
| `Tab` (twice to list) | Completes file and folder names |
| `↑` / `↓` | Previous / next command |
| `Ctrl+R` then type | Searches the command history (for example, `Ctrl+R` `gmsh`) |
| `Ctrl+C` | Stops the running command (for example, a z3st run) |
| `Ctrl+L` or `clear` | Clears the screen |
| `Ctrl+A` / `Ctrl+E` | Moves to the start / end of the line |
| `Ctrl+U` / `Ctrl+K` | Deletes up to the start / end of the line |
| `cd -` | Goes back to the previous folder |
| `cd ..`, `cd ~`, `pwd`, `ls -la` | Basic navigation |
| `code .` | Opens the current folder in VS Code |
| `code file.py` | Opens a file in the current window |
| `explorer.exe .` | Opens the current WSL folder in Windows Explorer |
| `conda activate z3st` | Activates the z3st environment |
| `tail -f log_z3st.md` | Follows the log while the case runs (`Ctrl+C` to stop following) |
| `python -m z3st > log_z3st.md 2>&1 &` | Runs in the background (`jobs`, `fg` to bring it back) |
| `grep -rn "Clamp_z" z3st/cases` | Searches text in files |
| `history \| grep gmsh` | Finds a previous command |

## C.5 Python and Jupyter

| Shortcut | What it does |
|---|---|
| `Ctrl+Shift+P` → *Python: Select Interpreter* | Picks the `z3st` conda env (do it once per workspace) |
| `F5` / `Ctrl+F5` | Runs the file with / without the debugger |
| `F9` | Toggles a breakpoint |
| `F10` / `F11` / `Shift+F11` | Step over / into / out (debugger) |
| `Shift+Enter` (in a `.py` file) | Sends the selected lines to the Python terminal |

Notebook (`.ipynb`). `Esc` enters command mode (blue bar). `Enter` goes back to editing.

| Key | What it does |
|---|---|
| `Shift+Enter` | Runs the cell and moves to the next one |
| `Ctrl+Enter` | Runs the cell and stays on it |
| `Alt+Enter` | Runs the cell and inserts a new one below |
| `A` / `B` | Inserts a cell above / below (command mode) |
| `D D` | Deletes the cell (command mode) |
| `Ctrl+Z` | Undoes the cell deletion (command mode) |
| `M` / `Y` | Turns the cell into Markdown / code (command mode) |
| `Shift+↑` / `Shift+↓` | Selects several cells (command mode) |
| Right-click the cell → *Join With Next Cell* | Merges cells |
| `Ctrl+Shift+\` | Splits the cell at the cursor (edit mode) |
| `L` | Shows / hides line numbers |
| `O` | Shows / hides the cell output |
| `Ctrl+/` | Comments lines inside a cell |

Restart the kernel from the notebook toolbar after editing a `.py` that the notebook imports
(for example, `k_1515.py`), or add `%load_ext autoreload` + `%autoreload 2` in the first cell.

---

## Common pitfalls

- A region name that is not in `labels` stops the run with `Region 'x' not found in label_map`.
- Neumann `flux`: a negative value means heat **enters** the body.
- Free body with only a thermal load: set `mechanical.remove_rigid_nullspace: true`, or add enough `Clamp_*` BCs.
- In 2D and axisymmetric runs the axial direction is **y**: use `Clamp_y`, not `Clamp_z`.
- In a script, sum `assemble_scalar` results with `comm.allreduce`. Write files only on `rank == 0`.
- `k(T)` callable: inside forms, T is a UFL object. Use `ufl.exp`, not `np.exp` (or keep it polynomial).
- Stress and strain are DG0. Compare them at cell centres, not at nodes.
