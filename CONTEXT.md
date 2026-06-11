# Z3ST — Repository Context

**Z3ST** is an open-source **FEniCSx-based finite-element framework** for coupled thermo-mechanical material analysis, written in Python. It is developed by **Giovanni Zullo** (Politecnico di Milano), licensed under **Apache 2.0**, version **0.1.0 (2025)**.

- **Repository:** https://github.com/giozu/z3st
- **DOI:** 10.5281/zenodo.17748028
- **Citation:** `CITATION.cff` + BibTeX entry in README
- **Continuous integration:** GitHub Actions (`.github/workflows/ci.yml`, `static.yml`)
- **Docs:** Sphinx sources under `docs/source/`, built by GitHub Actions
- **Python:** ≥ 3.10 (pyproject.toml)
- **Runtime dependencies:** `numpy`, `scipy`, `sympy`, `pandas`, `matplotlib`, `pyvista≥0.42`, `meshio`, `pyyaml`, `gmsh` + (external) **dolfinx / FEniCSx / basix / UFL / PETSc / MPI** provided by the Conda env `z3st_env.yml`.

---

## 1. Purpose and scope

Z3ST provides a **modular and extensible** environment to simulate:

- steady-state and transient **heat conduction** in arbitrary multi-material domains;
- **linear / hyperelastic / elasto-plastic** mechanical response;
- **phase-field fracture** (damage) with AT1 and AT2 models;
- **gap conductance** between domains (fixed or gas-type);
- **1D cluster-dynamics** advection–diffusion problems;
- arbitrary, spatially-dependent **internal heat sources** (γ-heating, user-defined);
- **fuel-performance behaviours carried at the material level**: burnup accumulation, radial power shaping `f(r, bu)`, burnup-driven swelling, and penalty **pellet-clad contact (PCMI)** with contact-coupled gap conductance;
- fully **YAML-driven** configuration for reproducibility, with the ability to plug in user **Python material modules** at runtime.

All physical models can be enabled/disabled independently via the `models:` block of `input.yaml`, and coupled physics are solved with a **staggered scheme** with adaptive relaxation.

---

## 2. High-level architecture

The code is a single installable package `z3st/` exposed as a CLI (`python -m z3st`) plus a set of post-processing utilities. Heavy `dolfinx`-dependent modules are lazily imported through `z3st/__init__.py::__getattr__`, so lightweight utilities can be imported without MPI/dolfinx installed.

### 2.1 Directory structure

```
z3st/
├── LICENSE                           Apache 2.0
├── README.md
├── CITATION.cff
├── GIT-COMMANDS.md                   internal git cheat-sheet
├── pyproject.toml / setup.py         installable Python package
├── z3st_env.yml                      Conda env recipe (FEniCSx + deps)
├── clean.sh                          cleanup helper
├── docs/                             Sphinx documentation
│   └── source/
│       ├── index.rst, api.rst, installation.rst, usage.rst,
│       │ physics_models.rst, examples.rst, quick_reference.rst,
│       │ differentiable_features.rst, troubleshooting.rst, ...
│       └── images/
├── .github/workflows/
│   ├── ci.yml                        non-regression CI
│   └── static.yml                    Sphinx docs build
└── z3st/                             main Python package
    ├── __init__.py                   lazy-import facade, exposes Spine/Solver/etc
    ├── __main__.py                   CLI entry point (runs cases from input.yaml)
    ├── core/                         FEM core infrastructure
    │   ├── config.py                 Config mixin (parses input.yaml)
    │   ├── finite_element_setup.py   FE function spaces (V_t, V_m, V_d, V_c, V_pl, Q)
    │   ├── solver.py                 staggered solver, PETSc options
    │   ├── spine.py                  Spine — top-level driver (multi-inheritance of mixins)
    │   ├── diagnostic.py             structured logger
    │   └── mesh/                     mesh sub-package
    │       ├── reader.py             Gmsh → dolfinx loader
    │       ├── manager.py            MeshManager (tags, geometry metrics, normals)
    │       └── plotter.py            PyVista mesh/tag preview
    ├── models/                       physical model mixins (plugged into Spine)
    │   ├── thermal_model.py          heat conduction + Dirichlet/Neumann/Robin BCs
    │   ├── mechanical_model.py       linear, Voigt, hyperelastic, plane_stress, custom
    │   ├── damage_model.py           phase-field AT1/AT2, Miehe/Amor splits
    │   ├── plasticity_model.py       J2 (return mapping) + custom crystal plasticity hook
    │   ├── gap_model.py              Fixed / Gas gap conductance (+ contact-coupled)
    │   ├── contact_model.py          penalty pellet-clad mechanical contact (gap closure)
    │   └── cluster_dynamic_model.py  1D advection-diffusion (DG+SIPG+upwind)
    ├── materials/                    YAML material cards + Python modules
    │   ├── steel.yaml, austenitic_steel.yaml, martensitic_steel.yaml,
    │   │ high_carbon_steel.yaml, T91.yaml, 15_15Ti.yaml,
    │   │ vessel_steel.yaml, vessel_steel_0.yaml
    │   ├── uo2.yaml, zircaloy.yaml
    │   ├── ceramic.yaml, oxide.yaml, plastic.yaml, lead.yaml, h2o.yaml
    │   ├── ceramic.py, oxide.py                    Python-side property callables (k(T), Gc(mesh))
    │   └── fuel_profiles.py, fuel_swelling.py      fuel-behaviour callables:
    │                                               radial power f(r,bu) + burnup-driven swelling(bu)
    ├── utils/                        post-processing + helpers
    │   ├── writer.py                 unified VTU / XDMF OutputWriter (pre-compiled, single-pass)
    │   ├── mesh_builder.py           build simple meshes from YAML geometry
    │   ├── plot_convergence.py       staggered-residual plots
    │   ├── utils_extract_vtu.py      scalar/vector/tensor extraction from VTU
    │   ├── utils_extract_xdmf.py     same for XDMF
    │   ├── utils_load.py             YAML loader + power-history generator
    │   ├── utils_plot.py             1D/radial plots (T(r), σ_rr(r), …)
    │   ├── utils_verification.py     reference analytical benchmarks
    │   ├── output.py                 stdout / JSON helpers
    │   ├── z-gui.py                  interactive PyVista 3D viewer
    │   ├── interactive_gui.ipynb     notebook front-end
    │   └── geo_files/                reusable Gmsh .geo templates
    │       ├── annulus_3D.geo, coaxial_cylinders_2D.geo,
    │       │ full_cylinder.geo,
    │       │ lenticular_bubble.geo, lenticular_bubble_3D.geo,
    │       │ mesh_1bubble.geo
    ├── examples/                     minimal didactic setups
    │   ├── thin_slab/
    │   └── cylindrical_shell/
    └── cases/                        verification / validation / demo cases
        ├── non-regression.sh, non-regression_github.sh
        ├── non-regression_github.py, non-regression_summary.txt
        └── <one folder per case — see §6>
```

### 2.2 The `Spine` driver

`z3st/core/spine.py::Spine` is the central orchestrator. It inherits from **all** physics and infrastructure mixins and exposes the high-level workflow invoked by `python -m z3st`:

```python
class Spine(Config,
            FiniteElementSetup,
            Solver,
            ThermalModel, MechanicalModel,
            GapModel,
            DamageModel,
            ClusterDynamicsModel,
            PlasticityModel):
    ...
```

Key responsibilities:

1. Load mesh (`load_mesh`) and wrap it in a `MeshManager` (cell tags, facet tags, geometry metrics: area/perimeter/inner_radius, facet normals, topological/facet dimensions).
2. Initialize each enabled model (from `self.on = {thermal, mechanical, damage, cluster, plasticity}`).
3. `load_materials(**materials)` — resolve material dictionaries, compute Lamé parameters from `E, ν`, bulk modulus, and automatically convert between `Gc ↔ σ_c` according to the active damage model (AT1 or AT2):
   - **AT2:** `Gc = (256/27) · lc · σ_c² / E`
   - **AT1:** `Gc = (8/3) · lc · σ_c² / E`
   Also recognises `k` or `Gc` given as **symbolic Python callables** (`module.function`) and resolves them via `importlib`.
4. `initialize_fields()` — instantiate `T`, `u`, `D`, `H`, `c`, … with initial conditions from material cards.
5. `set_boundary_conditions()` — read `boundary_conditions.yaml` and dispatch to thermal/mechanical/damage BC setters.
6. `set_power()` — build `q_third` (W/m³) per material: fissile (LHR/area), γ-heating with exponential decay in rect geometry or `K₀(μr)/K₀(μRᵢ)` in cylindrical geometry, spherical decay, etc. The cylindrical/spherical γ-attenuation profile is normalised at a **per-material reference surface** `gamma_inner_radius` (defaults to the geometry `inner_radius`, so existing cases are unchanged); a layer that sits inboard of the geometry reference — e.g. a thermal shield ahead of the vessel — sets it so its `K₀` profile is normalised at its own inner surface rather than the vessel's. **Radial power form factor (source bus):** a fissile material may carry a `radial_profile` callable (resolved like `k`/`Gc`), evaluated on the fuel dofs and area-normalised to mean 1, so the volumetric source is shaped `q''' = (LHR/area)·f(r, bu)` while preserving the integral. A built-in rim-peaking profile lives in `materials/fuel_profiles.py`; a mechanistic TUBRNP-style profile drops in behind the same hook.
7. `update_state(dt)` — advance each material's own history once per step (**state bus**, the `material.update_state(dt, fields)` channel). Currently: a fissile material accumulates its local **burnup** into the `self.burnup` field (MWd/kgU) from the deposited power, `Δbu = q_third·dt/(ρ·HM_frac·8.64e10)` (`heavy_metal_fraction` from the card). Called *before* the solve so fields at `t_k` are consistent with the state at `t_k` (behaviours that consume burnup — swelling, fuel-k — see the end-of-step value).
8. `solve(dt, max_iters)` — dispatch to `solve_staggered`.
9. `get_results()` — build symbolic UFL strain / stress / stress_mech / stress_th / energy_density dictionaries per material. The eigenstress `stress_th = −ℂ:ε*` carries the **total inelastic eigenstrain** ε* (thermal + material-carried swelling/creep — see §4.2), so fuel swelling enters equilibrium through the same channel as thermal expansion.
10. `compute_energy_balance()` (damage only) — assembles elastic + fracture energy (`E_el`, `E_frac`) for diagnostics written to `energies.txt`.

### 2.3 Execution flow (`z3st/__main__.py`)

```
read input.yaml
└── load geometry.yaml, material YAMLs
    └── Spine(input, mesh, geometry)
        ├── Config.__init__
        ├── FiniteElementSetup.__init__          # build V_t, V_m, V_d, V_c, V_pl, Q
        ├── Solver.__init__                       # relaxation factors, adaptive options
        └── Model.__init__ (thermal/mech/...)
    ├── load_materials
    ├── generate_power_history(times, lhrs, n_steps)
    ├── parameters(lhr), initialize_fields(), set_boundary_conditions()
    └── for each (t, lhr):
        ├── parameters(lhr); set_power()   # source bus: q''' = (LHR/area)·f(r, bu)
        ├── update_state(dt)               # state bus: accumulate burnup (before solve)
        ├── solve(max_iters, dt)           # staggered
        ├── get_results()                  # symbolic σ, ε, ψ, eigenstress
        ├── compute_energy_balance         # if damage
        ├── writer.write(t, step)          # VTU per-step files OR single XDMF time series
        └── diagnostics.per_step(...)      # optional case-local CSV hook (if diagnostics.py present)
```

**Output (`utils/writer.py::OutputWriter`).** The unified writer pre-allocates Function targets and compiles all field Expressions once at construction, so each `write()` is just interpolate + I/O (no UFL JIT in the loop). It emits a single set of **merged, domain-wide** fields — `Stress`, `VonMises`, `Hydrostatic`, `StrainEnergyDensity`, `HeatFlux` (plus `Temperature`, `Displacement`, `Strain`, `Damage`, `Burnup`, …) — each cell filled from its own material's expression, rather than one field per material. Format is selected by `output.format`: **`vtu`** writes per-step `fields_NNNN.vtu`; **`xdmf`** writes a single `fields.xdmf` + `.h5` time series for the whole run (the right choice for long transients — one file, opened directly in ParaView). Note: dolfinx XDMF is for ParaView / dolfinx read-back, not generic pyvista/meshio parsing. For scripted post-processing of long runs, a case-local `diagnostics.py` exposing `per_step(problem, step, t)` can stream scalar trajectories (e.g. burnup, gap, contact pressure, T_max) to a CSV — format-independent and step-count-independent.

CLI flags:
- `--debug`       enable verbose debugging
- `--mesh_plot`   preview surface tags with PyVista (via `core.mesh.plotter.MeshPlotter`)

**Markdown log filter** (`__main__.py`, top of file). When `sys.stdout` is *not* a TTY (i.e., the user redirected output, e.g. `python -m z3st > log.md`), a line-by-line stdout shim rewrites the existing decorated output into Markdown:
- Morse-code dividers `--.. ..- .-.. .-.. --- ...` → `***` (horizontal rule; `***` rather than `---` to avoid setext-heading ambiguity).
- `[STEP NN/MM] rest` → `## Step NN/MM: rest`.
- `--- Staggering iteration N/M ---` → `#### Iteration N/M`.
- `--. foo --..` (spine section markers) → `### foo`.
- `__Foo__` (initializer headers) → `### Foo`.
- `[DESCRIPTION]` → `## Description`.
- `[INFO|WARNING|ERROR|SUCCESS] body` → `**[TAG]** body`.

Interactive runs are pass-through. Set `Z3ST_PLAIN_LOG=1` to force the raw, unfiltered output even when redirected. None of the solver / model / config modules were touched; the filter is a single point of intercept in `__main__.py`.

**Hot-reload of input.yaml parameters** (`__main__.py`, `_reload_hot_params`). At the start of each time step, `__main__.py` re-reads `input.yaml` and propagates allow-listed parameter changes in-place into the in-memory config dicts (which are shared by reference with `problem.dmg_cfg` / `mech_cfg` / `th_cfg`, so the changes are immediately visible to the solver on the next step). The user can edit `input.yaml` mid-run and changes apply at the next step boundary (latency ≤ one step's wall-time). A one-line notice is printed when a value actually changes; silent on no-change steps. The reload is robust to mid-edit reads (transient `yaml.YAMLError` / `FileNotFoundError` → silent skip; the previous values stay in effect).

Allow-listed (hot-reloadable):
- `damage.{stag_tol, rtol, hybrid_constraint, gamma_star}`
- `mechanical.{stag_tol, rtol}`
- `thermal.{stag_tol, rtol}`
- `solver_settings.{max_iters, relax_T, relax_u, relax_D, relax_adaptive, relax_growth, relax_shrink, relax_min, relax_max}`

NOT hot-reloadable (intentionally — would invalidate pre-allocated FE structures or pre-compiled UFL Expressions): mesh / geometry / BC paths, materials, regime, `models.*` toggles, `damage.{type, lc, split}`, `mechanical.constitutive`, `plasticity.mode`, time history / `n_steps`. Edits to these are silently ignored mid-run; restart the simulation to apply them.

---

## 3. Core modules

### 3.1 `core/config.py — Config`

Parses the user YAML and fills:

- `self.on`: `{thermal, mechanical, damage, cluster, plasticity}` flags from `models:`
- `self.gap_model`, `self.h_gap_value` from `models.gap_conductance`
- paths to geometry, mesh, boundary conditions
- `self.n_steps`
- `self.regime ∈ {1d, 2d, 3d, axisymmetric, plane_stress}` (validated up front since 2026-06-10; `1d` is used by the teaching cases)

### 3.2 `core/finite_element_setup.py — FiniteElementSetup`

Allocates the FE function spaces on `self.mesh`:

| Space  | Element                               | Purpose                                 |
|--------|----------------------------------------|-----------------------------------------|
| `V_t`  | Lagrange P1 (scalar)                   | Temperature                             |
| `V_m`  | Lagrange P(order) vector (dim = tdim)  | Displacement (order from `mechanical.order`) |
| `V_d`  | Lagrange P1 (scalar)                   | Damage (when `damage` on)               |
| `Q`    | DG0 (scalar)                           | History / driving force H, per-cell     |
| `V_c`  | DG1 (scalar)                           | Cluster density (when `cluster` on)     |
| `V_pl` | quadrature (3×3 tensor), degree `q_degree` | Plastic strain tensor ε_p (when `plasticity` on) |
| `Q_pl` | quadrature (scalar), degree `q_degree`  | Cumulative plastic strain p             |

### 3.3 `core/solver.py — Solver`

Implements the staggered solver and PETSc options.

**PETSc options** via `get_solver_options(physics, solver_type, rtol)`:
- `direct_mumps`  → LU + MUMPS (`preonly`)
- `iterative_amg` → CG/GMRES + GAMG preconditioner
- `iterative_hypre` → CG/GMRES + HYPRE BoomerAMG

KSP selection: **CG** for thermal and damage (SPD), **GMRES** for mechanical.

**Integration measures** (`_build_measures`):
- axisymmetric ⇒ `w = 2π r` (spatial coordinate)
- 2D / 3D ⇒ `w = 1`
- optional `quadrature_degree` (required for plasticity because of quadrature element spaces).
- builds `self.dx_tags[tag]` and `self.ds_tags[id]` keyed per cell/facet tag.

**`_thermal_step`** (backward Euler if `thermal.analysis = transient`):
```
a_t(u,v) = ∫ w k ∇u·∇v dx   (+ ∫ w (ρ cp/dt) u v dx if transient)
         + ∫ w h u v ds                (Robin: gap or convective)
L_t(v)   = ∫ w q''' v dx
         + ∫ w (ρ cp/dt) Tⁿ v dx       (transient)
         + ∫ w (−flux) v ds            (Neumann)
         + ∫ w h T_ext v ds            (Robin)
```
- Supports `dt = 0`: preserves IC, only applies BCs.
- Post-solve relaxation: `T ← α_T · T_new + (1 − α_T) · T_old`.
- Convergence on `‖ΔT‖` or `‖ΔT‖/‖T‖` (L2) depending on `thermal.convergence`.
- Adaptive relaxation (EMA residual) scales `α_T` between `relax_min` and `relax_max`.

**`_mechanical_step`** (linear or SNES Newton for non-linear/hyperelastic):
- Updates step-dependent Dirichlet displacements and tractions (`raw` list indexed by `current_step`).
- Body force `(0, -ρg)` (2D/axisym) or `(0,0,-ρg)` (3D).
- Linear form:
  ```
  a_m = ∫ w σ_mech(u):ε(v) dx     L_m = ∫ w b·v dx + ∫ w t·v ds  − ∫ w σ_th(T):ε(v) dx
  ```
- Non-linear `F_m` residual assembled with `hyperelastic_residual` when `constitutive = hyperelastic`, otherwise via `sigma_mech(u, mat)` at the current displacement.
- SNES options: `newtonls`, line-search `basic` or `bt`, inner solver MUMPS / AMG / HYPRE.

**`_damage_step`** — linear elliptic problem per iteration (AT1 or AT2):

- **AT2:**
  `a = ((H + 1) u v + lc² ∇u·∇v) w dx`,  `L = H v w dx`.
- **AT1:**
  `a = (2 H + diag_shift) u v w dx + (3 Gc lc / 4) ∇u·∇v dx`,
  `L = (2 H − 3 Gc/(8 lc)) v w dx`.
  Produces the sharp elastic threshold (`σ > σ_c`) and requires a post-solve clipping `D ∈ [0,1]`.

In both cases, post-solve **irreversibility** is enforced pointwise:
`D_new ← max(D_new, D_old)` and clipped to `[0,1]`.

**`_cluster_step`** (DG upwind + SIPG):

Solves, per time step (implicit Euler) on a 1D advection–diffusion:
```
∂c/∂t + v ∂c/∂n − D ∂²c/∂n² = 0
```
with upwind interior-facet flux for advection, SIPG for diffusion, and a **mass-conservation rescaling** of the form `C_tot = ∫ c·n dn = const`. Logs Péclet number for diagnostics.

**`solve_staggered(max_iter, dt, rtol_*, stag_tol_*)`** — outer loop:
1. `_build_measures`
2. allocate local copies `T_new/old`, `u_new/old`, `D_new/old`, `c_new` if the corresponding model is active.
3. for each iteration: thermal → mechanical → (`update_history(u)` →) damage → cluster
4. check `conv_th ∧ conv_mech ∧ conv_damage`; on success push local solutions back into `self.T, self.u, self.D, self.c` and trigger `update_plastic_history(u)` when plasticity is active.
5. if `max_iter` exceeded, keep last iterate and warn.

### 3.4 `core/mesh/*`

- `reader.py` — loads Gmsh `.msh` files into `dolfinx.mesh` with `cell_tags` and `facet_tags`.
- `manager.py` — `MeshManager` stores geometry metadata (`geometry_type ∈ {rect, cyl/cylinder, sphere}`, `Lx, Ly, Lz, Ri, Ro`), the `label_map` from `geometry.yaml`, computes area/perimeter/inner_radius and exposes `locate_domain_dofs`, `locate_facets_dofs`, facet normals, dimensions.
- `plotter.py` — optional PyVista viewer for surface tags.

### 3.5 `core/diagnostic.py`

Custom structured logger (`log.info`, `log.warning`, …) reused across modules.

---

## 4. Physical models (`z3st/models/`)

### 4.1 Thermal (`thermal_model.py`)

- Supported BCs: `Dirichlet`, `Neumann` (flux), `Robin` (two modes — gap by `pair:` or convective by `h_conv + T_ext`).
- `heat_flux(T)` diagnostic: per-material average `|q|` + per-component flux (dimension-aware, `r`/`z` labels in axisymmetric; supports symbolic `k(T)` cards). Printed after each step's `get_results()` **only under `--debug`** — the writer's `HeatFlux` field is the default output channel.
- Config keys read from `input.yaml::thermal`:
  `solver`, `linear_solver`, `rtol`, `stag_tol`, `convergence`, `analysis` (`stationary | transient`).

### 4.2 Mechanical (`mechanical_model.py`)

Constitutive modes (via `material.constitutive` — see §5):

| mode            | Behaviour                                                            |
|-----------------|----------------------------------------------------------------------|
| `lame` (default)| Isotropic small-strain `σ = λ tr(ε) I + 2 G ε`                       |
| `voigt`         | Uses user-provided 6×6 `C_matrix`, or isotropic fallback             |
| `hyperelastic`  | Neo-Hookean `ψ = μ/2(I_C − 3) − μ ln J + λ/2 (ln J)²`, σ = (1/J) P F⊤ |
| `plasticity`    | J2 return-mapping with linear isotropic hardening                     |
| `custom`        | Loads `material.stress_function = "pkg.mod.func"`, `σ = f(u, T, material, model)` |

Regime handling in `epsilon(u)`:
- **axisymmetric:** strain in cylindrical coords `(r, θ, z)`; `ε_rr = ∂u_r/∂r`, `ε_θθ = u_r/r`, `ε_zz = ∂u_z/∂z`, `ε_rz = ½(∂u_r/∂z + ∂u_z/∂r)`.
- **2D:** 3×3 tensor with zero z-components (plane strain).
- **3D:** `sym(∇u)`.
- **plane_stress:** replaces `λ` with `λ_ps = 2 G λ / (λ + 2G)` and forces `σ_zz = σ_xz = σ_yz = 0`.

BC types (`set_mechanical_boundary_conditions`): `Dirichlet`, `Dirichlet_x/y/z`, `Neumann` (scalar traction along facet normal), `Clamp_x/y/z`, `Slip_x/y/z`. Step-dependent BCs supported via a list of values of length `n_steps`.

**Eigenstrain bus (thermal + material inelastic strains).** Equilibrium puts an eigenstress on the RHS, `σ_th = −ℂ:ε*`, where `ε*` is the *total* inelastic eigenstrain returned by `MechanicalModel.eigenstrain(T, material)`:
- thermal expansion `α (T − T_ref) I` (when a thermal field is active);
- a constant volumetric swelling `(swelling/3) I` from a scalar `swelling` card field;
- a state-dependent material eigenstrain via an `eigenstrain` callable (`"pkg.mod.func"`, resolved like `k`/`Gc`) — e.g. burnup-driven fuel swelling `(swelling_rate·bu/3) I` in `materials/fuel_swelling.py`, which reads the `burnup` field.

For a purely thermal eigenstrain this reduces exactly to `σ_th = −(3λ + 2G) α (T − T_ref) I`. Because ε* is a UFL tensor the Newton tangent stays automatic, and fuel swelling/creep need no change to the momentum balance — *"fuel is a material"*: each region's inelastic behaviour travels with its own material, applied wherever the material's thermal block or own eigenstrain is active.

Damage coupling: when damage is active, `σ ← g(D) σ` with `g(D) = (1−D)² + K` (K = 1e-6 regularization). The eigenstress is degraded by the same `g(D)` so a fully-damaged cell recovers the traction-free crack-face limit.

### 4.3 Damage (`damage_model.py`)

Phase-field fracture with two variational models:

- **AT2**: quadratic local dissipation `w(D)=D²`, no elastic threshold.
  `H = (2 lc / Gc) · ψ⁺`  (non-dimensional history).
  Miehe spectral split (positive eigenvalues of `ε`).
- **AT1**: linear local dissipation `w(D)=D`, analytical elastic threshold `σ_c`.
  `H = ψ⁺` (physical, J/m³).
  Amor (volumetric/deviatoric) split.

Elastic-energy splits (selectable via `damage.split: amor | miehe | star_convex`; if absent, defaults to Amor for AT1 and Miehe for AT2 — the historical pairing):
- `psi_miehe_spectral(u, mat)` — 2D closed-form + 3D via Cardano's formula with smooth clamping.
- `psi_amor_split(u, mat)` — `ψ⁺ = ½ λ ⟨tr ε⟩₊² + G dev(ε):dev(ε)`; `ψ⁻ = ½ λ ⟨tr ε⟩₋²`.
- `psi_star_convex(u, mat)` — Vicentini, Zolesi, Carrara, Maurini, De Lorenzis 2024 (Int. J. Fract. 247:291-317). One-parameter generalisation of Amor controlled by `dmg_cfg["gamma_star"]` (a *model* parameter in the `damage:` block, default 0 → reduces to Amor; intentionally not a per-material property). `ψ⁺ = G|dev ε|² + (λ/2)[⟨tr ε⟩₊² − γ⋆⟨tr ε⟩₋²]`, `ψ⁻ = (1+γ⋆)(λ/2)⟨tr ε⟩₋²`. Satisfies all five criteria in the Vicentini 2024 Table 2 (the other splits do not). `γ⋆ > 0` raises the compressive-vs-tensile critical-stress ratio.

`update_history(u)` — vectorised, per-material update of the history field `H` on DG0. Supports **Ambati-Gerasimov-De Lorenzis hybrid constraint** (`dmg_cfg.hybrid_constraint`, default `True`): where `ψ⁻ > ψ⁺` locally, contribution to H is set to 0 to suppress crack growth in compression.

`compute_energy_balance(u)` returns `(E_el, E_frac)`:
- AT2: `γ = Gc/2 · (D²/lc + lc |∇D|²)`
- AT1: `γ = 3 Gc/8 · (D /lc + lc |∇D|²)`

BCs: only `Dirichlet` on D (e.g. `D=0` on healthy boundary).

### 4.4 Plasticity (`plasticity_model.py`)

- **J2 plasticity** with isotropic linear hardening; constitutive update computed symbolically through `ufl.conditional` (returns stress `σ = σ_tr − 3μ Δp n`, yield stress `σ_y = σ_y0 + H p`).
- History fields `ep`, `p`, `ep_n`, `p_n` on quadrature spaces.
- `mode: custom` in `plasticity:` block hooks a user crystal-plasticity function `get_cp_internal_variables` inside the material's module — used by the `demo_CP_single_grain` case.
- `update_plastic_history(u)` refreshes `ep_n`, `p_n` after the staggered iteration converges.

### 4.5 Gap conductance (`gap_model.py`)

Two modes (selected via `models.gap_conductance.type`):

- **Fixed** — constant `h_gap` (W/m²K) read from `value:`.
- **Gas** — `k_gas = value · 1e-4 · T_gap^0.79`, then `h_gap = k_gas / gap_size` where `gap_size` is computed as the mean distance between two paired labelled facet groups (via SciPy cKDTree on facet centroids), and `T_gap = ½ (T_inner + T_outer)`.

Invoked inside `_thermal_step` when a Robin BC is defined with `pair:` to another subdomain.

**Contact-coupled conductance.** When `gap_conductance.contact_coupling.enabled` is set, a solid-contact term is added on gap closure (Todreas & Kazimi, *Nuclear Systems I*, 3rd ed., Eqs. 8.141/8.142): the emergent contact pressure (from `contact_model`) raises `h_gap` above the open-gap gas value, so closing the gap cools the fuel. Parameters: `meyer_hardness` (Pa), `gas_thickness` (m, roughness-based residual gas space).

### 4.5bis Creep (`creep_model.py`)

Implicit Norton creep (`ε̇_eq = A0·exp(−Q/RT)·σ_eq^n`) for a material carrying `creep: norton` + `creep_A0/n/Q` on its card — the dissipative extension of the energy-first design (incremental variational principle, Ortiz–Stainier). The cell-local minimisation over Δε_cr condenses to the scalar radial-return equation per point; a DG0 **predictor** Δγ₀ holds its exact root (vectorised numpy Newton, refreshed before every mechanical solve, consistency gated in the staggered convergence test), and the UFL stress carries **one symbolic Newton step** from the predictor — so `ufl.derivative` yields exactly the implicit-function-theorem consistent tangent through a trivially small expression tree (a fully unrolled symbolic Newton explodes FFCx). The accumulated `ε_cr` is a per-material DG0 tensor state advanced once per converged step; it enters the trial through the eigenstrain channel and the output stress via `creep_output_stress`. Mechanical steps auto-promote to SNES when creep is active. v1 scope: isotropic Lamé, no damage/plasticity combination, regimes with 3×3 strain tensors. Verified by `V_creep_verification` (1e-14) and `V_creep_relaxation_verification` (4e-15 vs the BE recursion; O(dt) defect pinned).

### 4.6 Cluster dynamics (`cluster_dynamic_model.py`)

1D advection–diffusion solver for defect-cluster size distributions `c(n,t)`:
- `∂c/∂t = −v ∂c/∂n + D ∂²c/∂n²`
- initial conditions: `constant` (on labelled region) or `gaussian`.
- DG1 space with upwind advection and SIPG diffusion.
- Mass conservation: `∫ c·n dn` is rescaled to the initial target every step.

### 4.7 Mechanical contact (`contact_model.py`)

Penalty pellet-clad mechanical contact (gap closure / PCMI), enabled via the `models.contact` block. An explicit fixed-point scheme integrated into the staggered loop:

- **Gap measurement** — each mechanical iteration measures the current normal gap from the displacement iterate as `gap = g0 + ⟨u_r⟩_b − ⟨u_r⟩_a`, the boundary-integral mean radial displacement of the two paired facing surfaces (`surface_a` = pellet outer, `surface_b` = clad inner).
- **Penalty traction** — on penetration (`gap < 0`) a pressure `p = k_pen · ⟨−gap⟩₊` is applied as `t = −p·n` on both facing surfaces (UFL/AD supplies the tangent). Config: `penalty_stiffness` (Pa/m), `initial_gap` (m).
- Verified against the analytical plane-stress Lamé interference-fit solution (`cases/V_coaxial_contact_verification`, ~3.5 %). The emergent contact pressure also feeds the contact-coupled gap conductance (§4.5), so thermal + mechanical PCMI are two-way coupled.

---

## 5. Material database (`z3st/materials/`)

Materials are plain YAML cards. Common fields:

| Key                 | Units      | Meaning                                                |
|---------------------|------------|--------------------------------------------------------|
| `name`              | —          | Human-readable name                                    |
| `E`, `nu`           | Pa, —      | Young's modulus, Poisson's ratio                        |
| `k`                 | W/(m·K) or `"mod.func"` | Thermal conductivity (scalar OR symbolic)    |
| `cp`, `rho`         | J/(kg·K), kg/m³ | Specific heat, density                             |
| `alpha`             | 1/K        | Thermal expansion coefficient                          |
| `T_ref`             | K          | Stress-free / reference temperature                    |
| `T_initial`         | K          | Initial temperature in the field                        |
| `mu_gamma`          | 1/m        | γ-ray attenuation coefficient                          |
| `gamma_heating`     | W/m³       | Volumetric γ heating magnitude                          |
| `gamma_inner_radius`| m          | Per-material reference radius for the cylindrical/spherical γ-attenuation profile (optional; defaults to geometry `inner_radius`; set it for an inboard layer such as a shield) |
| `fissile`           | bool       | If true, `q''' = LHR / area` in the pellet              |
| `heavy_metal_fraction` | —       | M_U / M_compound (e.g. 0.8815 for UO2); burnup accumulated per kg heavy metal |
| `radial_profile`    | string     | `"pkg.mod.func"` radial power form factor `f(r, bu)` (source bus); see `fuel_profiles.py` |
| `radial_peak_amplitude` / `radial_peak_exponent` | — | parameters of the built-in `rim_peaking` profile `f = 1 + A(r/Ro)^p` |
| `axial_profile`     | string     | `"pkg.mod.func"` axial power form factor `f(z)` (source bus, composed with the radial one); built-ins in `fuel_profiles.py`: `chopped_cosine`, `tabulated_axial` |
| `axial_extrapolated_length` | m  | extrapolated length L′ of the chopped cosine `f(z) = cos(π(z−z_mid)/L′)` (default 1.1·L; peaking factor = 1/[(2L′/πL)·sin(πL/2L′)]) |
| `axial_table_z` / `axial_table_f` | m, — | elevation/factor lists for `tabulated_axial` (piecewise-linear, end values held outside the range); only the *shape* matters — the mean-1 normalisation makes the absolute scale irrelevant |
| `swelling`          | —          | constant volumetric swelling ΔV/V (isotropic eigenstrain `(ΔV/V)/3 · I`) |
| `eigenstrain`       | string     | `"pkg.mod.func"` state-dependent eigenstrain callable (e.g. swelling(bu)); see `fuel_swelling.py` |
| `swelling_rate`     | (MWd/kgU)⁻¹ | ΔV/V per unit burnup for `fuel_swelling.solid_swelling`  |
| `sigma_c` / `Gc`    | Pa, J/m²   | Phase-field critical stress OR fracture energy (one is derived from the other given `lc`) |
| `yield_strength`    | Pa         | Initial yield stress (J2 plasticity)                    |
| `hardening_modulus` | Pa         | Linear isotropic hardening modulus                      |
| `constitutive`      | string     | `lame`, `voigt`, `hyperelastic`, `plasticity`, `custom` |
| `C_matrix`          | 6×6        | User-provided elasticity matrix (for `voigt` mode)      |
| `stress_function`   | string     | `"pkg.mod.func"` for `constitutive: custom`             |

The framework auto-fills `lmbda`, `G`, `bulk_modulus` from `(E, ν)` at material load time.

Available cards (non-exhaustive):

- **Steels:** `steel.yaml`, `austenitic_steel.yaml`, `martensitic_steel.yaml`, `high_carbon_steel.yaml`, `T91.yaml`, `15_15Ti.yaml`, `vessel_steel.yaml`, `vessel_steel_0.yaml`
- **Materials:** `uo2.yaml` (E = 358 GPa, k = 5 W/mK, α = 1e-5/K, Gc = 15 kJ/m², T_initial = 1023 K), `zircaloy.yaml` (Zircaloy-4)
- **Ceramics / oxides:** `ceramic.yaml` (+ Python `ceramic.py::k(T)`), `oxide.yaml` (+ Python `oxide.py::k(T), Gc(mesh)` with grain-boundary heterogeneity via `tanh(|y|/half_width)`)
- **Other:** `plastic.yaml`, `lead.yaml`, `h2o.yaml`
- **Fuel-behaviour callables:** `fuel_profiles.py::rim_peaking` (radial power form factor `f(r, bu)` — source bus), `fuel_swelling.py::solid_swelling` (burnup-driven swelling eigenstrain — eigenstrain bus). These realise the *"fuel is a material"* design: state-dependent fuel physics (radial power, swelling; later densification, fuel-k(bu), creep) live in the material and travel with its region, rather than as global solver toggles.

---

## 6. Simulation cases (`z3st/cases/`)

Each case folder is self-contained:

```
<case>/
├── Allrun        shell driver (gmsh → python -m z3st → non-regression.py → plot)
├── Allclean      cleanup helper
├── input.yaml    simulation configuration
├── geometry.yaml geometry definition + label ↔ tag map
├── boundary_conditions.yaml  thermal / mechanical / damage BCs
├── mesh.geo      Gmsh input
├── mesh.msh      generated mesh (NOT tracked in git per .gitignore — Allrun regenerates it with gmsh on each run)
├── non-regression.py         post-run comparison vs. output/non-regression_gold.json
└── output/                   auto-generated VTU/XDMF + plots
```

The suite is driven by `z3st/cases/non-regression.sh` (local) and `non-regression_github.sh` (CI) and summarised in `non-regression_summary.txt`. Each case's `non-regression.json` carries two verdicts (since 2026-06-10): `"summary"` (analytic-tolerance check) and `"regression"` (vs the blessed `non-regression_gold.json`); the local summary reports both per case, and CI fails when either is FAIL (previously only Allrun crashes failed CI — numerical regressions were invisible).

### 6.1 Catalogue of cases

**00 — Tutorial / example**
- `00_example` — uniaxial traction of a rectangular steel block (3D, linear elasticity, `Neumann` traction + `Clamp` BCs on 3 faces).

**1–6 — Slabs and thin cylindrical shells (thermal / thermo-mechanical benchmarks)**
- `1_thin_slab_2D`, `1_thin_slab_neumann_2D`, `1_thin_slab_neumann_3D`, `1_thin_slab_non_linear`
- `2_thin_cylindrical_shell_2D`
- `3_thick_slab_adiabatic_3D`, `3_thin_slab_adiabatic_2D`
- `4_thick_cylindrical_shell_adiabatic_2D`, `4_thin_cylindrical_shell_adiabatic_2D`
- `5_thick_slab_non_adiabatic_3D`, `5_thin_slab_non_adiabatic_2D`
- `6_thick_cylindrical_shell_non_adiabatic_2D`, `6_thin_cylindrical_shell_non_adiabatic_2D`

**7 — Box heated**
- `7_box_heated`

**8–9 — Thick cylinders with regime variations**
- `8_thick_cylindrical_shell_plane_strain_2D`
- `9_thick_cylindrical_shell_GPS_2D`, `9_thick_cylindrical_shell_GPS_3D` (generalized plane strain)

**11–14 — Cylinders (Mariotte, gradients, annular, full cylinder, fracture)**
- `11_thin_cylindrical_shell_Mariotte`
- `12_cylindrical_shell_thermal_gradient_2D`, `12_cylindrical_shell_thermal_gradient_3D`
- `13_annular_cylinder`
- `14_full_cylinder`
- `14_full_cylinder_cracking` (3D gold-standard McClenny reproducer), `14_full_cylinder_cracking_2D_xy` (plane-strain McClenny Fig. 8 reproducer — the workhorse for the paper's case-14 chapter), `14_full_cylinder_thermal_2D_rz` (axisymmetric thermal verification, damage off). See §9c for the full three-variant rationale.

**15 — Cavities and pressurised bodies**
- `15_single_elliptical_cavity_2D`, `15_two_elliptical_cavities_2D`
- `15_spherical_pressurized_cavity`

**16 — Multi-body coupling**
- `16_coaxial_cylinders_3D`

**17 — Stress–strain curves**
- `17_stress_strain_curve_displacement`, `17_stress_strain_curve_stress`
- `17_stress_strain_curve_double_crack`, `17_stress_strain_curve_knotch`

**18 — 2D fracture benchmarks**
- `18_box_crack_2D`, `18_box_knotch_2D`

**19 — Single-edge notched (classical phase-field benchmarks)**
- `19_single-edge_notched_shear_test`
- `19_single-edge_notched_tension_test`

**20 — Plasticity**
- `20_plasticity_2D`

**I, II — Utilities**
- `I_mesh_sensitivity_2D` — mesh convergence study
- `II_attenuation_map` — γ attenuation in materials

**V_* — Analytical-verification cases** (closed-form checks; `V_` = verification, renamed from `U_*` on 2026-06-09)
- `V_swelling_verification` — constant volumetric swelling eigenstrain (free expansion → σ ≈ 0, exact `u`).
- `V_fuel_swelling_verification` — burnup-driven swelling reading the `burnup` field (the eigenstrain bus consuming the state bus).
- `V_burnup_verification` — burnup accumulation + radial-power source bus on an axisymmetric pellet (closed-form mean burnup; rim/core ratio = 1 + A).
- `V_axial_power_verification` — axial-power source bus (chopped-cosine `f(z)`, Todreas & Kazimi 1-D axial problem) on a tall axisymmetric fuel column: closed-form mean burnup (machine precision); axial peaking factor = 1/[(2L′/πL)·sin(πL/2L′)]; end/peak = cos(πL/2L′).
- `V_axial_table_verification` — tabulated axial profile (`tabulated_axial`, piecewise-linear node-wise peaking factors — the standard core-physics input): closed-form mean burnup (exact); table-node ratio f₃/f₁ (machine precision); peak/mean = max f / trapezoid mean.
- `V_creep_verification` — implicit Norton creep, constant-stress uniaxial bar (backward Euler exact): total/creep/radial strain vs closed form at 1e-14 (radial pins the deviatoric −½ flow).
- `V_creep_relaxation_verification` — stress relaxation at held strain: Z3ST ≡ the scalar backward-Euler recursion at 4e-15; deviation from the exact `σ(t)` equals the predicted O(dt) defect (2.11% at 50 steps, pinned to 2e-13).
- `V_coaxial_contact_verification` — penalty contact pressure vs the analytical plane-stress Lamé interference-fit.

**U_* — Extended / demo cases**
- `U_coaxial_contact_2D` — 2D-rz PCMI penalty-contact demo (oxide pellet + steel clad, power ramp → gap closure, contact pressure).
- `U_pwr_rod_2D` — generic-PWR fuel-rod segment (4.5 mm pellet, 65 µm gap, Zircaloy clad): coupled thermal + mechanical + gap conductance + penalty contact + burnup + swelling → **burnup-driven gap closure and PCMI** over a multi-year power history. BCs include the 15.5 MPa coolant pressure on the clad outer and (since 2026-06-10) a 2 MPa He fill-gas pressure on all gap/plenum-facing surfaces (`lateral_1`, `top_1`, `inner_2`); coolant is still a Dirichlet 580 K (no film drop — see punch list CODE-FEATURE-3 for the coolant module). Gold-protected since 2026-06-10: `non-regression.py` reads end-state PCMI scalars from `output/history.csv` (PCMI onset 30.4 MWd/kgU; final gap −0.78 µm, p = 38.8 MPa, bu = 44.5 MWd/kgU) with the mean burnup checked against the closed form; in the local suite, not in CI (~3 min run).
- `U_pressure_vessel_2D`, `U_cluster_dynamics_test`, `U_quarter_block`, `U_spherical_shell`. (`U_box_knotch_3D`, `U_slab_contact`, `U_thick_cylindrical_shell_plane_stress` were removed in commit f1bb70b; note this leaves the `plane_stress` regime with no exercising case.)

**demo_CP_single_grain** — crystal-plasticity single-grain demo using the `custom` constitutive + `plasticity.mode: custom` hook.

### 6.2 Typical `input.yaml` structure

```yaml
mesh_path: mesh.msh
geometry_path: geometry.yaml
boundary_conditions_path: boundary_conditions.yaml

materials:
  uo2: ../../materials/uo2.yaml

regime: axisymmetric            # 2d | 3d | axisymmetric | plane_stress

solver_settings:
  coupling: staggered
  max_iters: 200
  relax_T: 0.8
  relax_u: 0.6
  relax_D: 0.5
  relax_adaptive: true
  relax_growth: 1.1
  relax_shrink: 0.9
  relax_min: 0.1
  relax_max: 0.95

models:
  thermal: true
  mechanical: true
  damage: true
  # plasticity: true
  # cluster: true
  # gap_conductance:
  #   type: Fixed              # or Gas
  #   value: 5000.0
  #   contact_coupling: { enabled: true, meyer_hardness: 9.65e8, gas_thickness: 4.0e-6 }
  # contact:                   # penalty pellet-clad mechanical contact (PCMI)
  #   surface_a: pellet_outer  # facet group on body A
  #   surface_b: clad_inner    # facet group on body B
  #   penalty_stiffness: 5.0e13
  #   initial_gap: 65.0e-6

# Fissile / fuel behaviour is configured on the MATERIAL cards, not here:
#   fissile, heavy_metal_fraction, radial_profile (source bus),
#   swelling / eigenstrain + swelling_rate (eigenstrain bus). See §5.

thermal:
  analysis: transient          # stationary | transient
  solver: linear
  linear_solver: iterative_hypre   # direct_mumps | iterative_amg | iterative_hypre
  rtol: 1.0e-6
  stag_tol: 1.0e-3
  convergence: rel_norm        # rel_norm | norm

mechanical:
  solver: linear               # linear | nonlinear (SNES)
  linear_solver: iterative_hypre
  rtol: 1.0e-6
  stag_tol: 1.0e-3
  convergence: rel_norm
  order: 1                     # FE order

damage:
  type: AT2                    # AT1 | AT2
  solver: linear
  linear_solver: iterative_hypre
  rtol: 1.0e-6
  stag_tol: 1.0e-3
  convergence: rel_norm
  lc: 1.0e-4
  hybrid_constraint: true
  # split: star_convex          # optional; amor | miehe | star_convex.
                                # When star_convex, optionally set gamma_star >= -1
                                # (defaults to 0 -> reduces to Amor exactly).
  # gamma_star: 1.0             # star-convex model parameter; see damage_model.py::psi_star_convex

time:  [0.0, 0.01]
lhr:   [0.0, 0.0]
n_steps: 20

output:
  format: xdmf                 # vtu | xdmf  (vtkhdf prototyped but reverted; dolfinx 0.10's vtkhdf submodule lacks per-field naming)
  filename: simulation.xdmf
```

### 6.3 Typical `boundary_conditions.yaml`

```yaml
thermal:
  uo2:
    - { type: Dirichlet, region: contact_wall, temperature: 263.15 }

mechanical:
  uo2:
    - { type: Clamp_y, region: bottom }
    - { type: Clamp_x, region: axis }
```

Thermal BC types: `Dirichlet`, `Neumann` (flux), `Robin` (either `pair:` for gap-coupled or `h_conv + T_ext` for convective).
Mechanical BC types: `Dirichlet` (vector), `Dirichlet_x/y/z`, `Neumann` (scalar pressure along facet normal), `Clamp_x/y/z`, `Slip_x/y/z`. Lists of length `n_steps` give step-dependent histories.
Damage BC types: `Dirichlet` (`D = const`).

---

## 7. Capabilities summary

| Capability                        | Status / notes                                                     |
|-----------------------------------|---------------------------------------------------------------------|
| Stationary heat conduction        | ✓ (linear)                                                          |
| Transient heat conduction         | ✓ backward-Euler                                                    |
| Linear elasticity                 | ✓ Lamé / Voigt / plane-stress                                       |
| Anisotropic elasticity            | ✓ via user-provided 6×6 `C_matrix`                                  |
| Hyperelasticity                   | ✓ Neo-Hookean (SNES Newton)                                         |
| Thermo-mechanical coupling        | ✓ staggered, adaptive relaxation                                    |
| Phase-field fracture (AT1, AT2)   | ✓ Miehe/Amor split, hybrid constraint, irreversibility              |
| J2 plasticity                     | ✓ return mapping + linear isotropic hardening (quadrature elements) |
| Crystal plasticity                | experimental — via `custom` constitutive hook (`demo_CP_single_grain`) |
| Gap conductance                   | ✓ Fixed or Gas (k_gas = f(T_gap), gap_size from facet centroids)    |
| Cluster dynamics (1D)             | ✓ DG upwind + SIPG, mass-conservation renormalisation                |
| Axisymmetric / 2D / 3D / plane-stress regimes | ✓ all consistent with the integration weight `w`      |
| Volumetric heating                | ✓ fissile (LHR/area), γ-heating (rect / cyl / sphere analytic decay), user `q'''` |
| Burnup accumulation               | ✓ per-fissile-material `burnup` field via `update_state(dt)` (state bus)            |
| Radial power shaping              | ✓ `radial_profile` form factor `f(r, bu)` (source bus); built-in rim-peaking        |
| Axial power shaping               | ✓ `axial_profile` form factor `f(z)` (source bus, composed `f_r·f_z`, single mean-1 normalisation); built-ins: chopped cosine (T&K), tabulated (node-wise peaking factors) |
| Cladding creep (implicit, AD)     | ✓ Norton + Arrhenius via the incremental variational principle (`models/creep_model.py`): condensed radial return on the displacement space, DG0 predictor + one symbolic Newton step → exact IFT consistent tangent by `ufl.derivative`; per-material `ε_cr` DG0 state; card keys `creep: norton`, `creep_A0/n/Q`; verified to 1e-14 (constant stress) and 4e-15 vs the BE recursion (relaxation) |
| Integrated-power diagnostic       | ✓ `set_power` prints the exact FE integral of the fissile source per material per step (regime-weighted, MPI-reduced); note the mean-1 normalisation is *nodal*, so a radially peaked profile integrates to LHR·Lz·⟨f⟩_area/⟨f⟩_nodal (= 1.2·LHR·Lz for rim-peaking A=3, p=8) — pinned by the `total_power` checks in the burnup-family `V_` cases |
| Fuel swelling                     | ✓ constant ΔV/V or burnup-driven eigenstrain (eigenstrain bus)                       |
| Pellet-clad contact (PCMI)        | ✓ penalty contact + contact-coupled gap conductance (verified vs analytical Lamé)   |
| Python material callables         | ✓ `k(T)`, `Gc(mesh)`, `radial_profile(r,bu)`, `eigenstrain(bu)` loaded via `importlib` |
| Post-processing                   | unified `OutputWriter` (merged domain-wide fields; per-step VTU **or** single-file XDMF time series), case-local diagnostics CSV, radial/1D plots, PyVista viewer, notebook GUI |
| Mesh IO                           | Gmsh (`.msh`), YAML-based mesh builder                              |
| Parallelism                       | MPI via `MPI.COMM_WORLD`; PETSc (MUMPS / GAMG / HYPRE BoomerAMG)     |
| CI / non-regression               | per-case `non-regression.py` vs. JSON gold, summarised in `non-regression_summary.txt` |

---

## 8. Running a case

```bash
git clone https://github.com/giozu/z3st.git
conda env create -f z3st_env.yml && conda activate z3st
pip install -e .

cd z3st/cases/00_example/
gmsh -3 mesh.geo                    # or: ./Allrun
python -m z3st > log_z3st.md
python non-regression.py
python ../../utils/plot_convergence.py
```

Optional flags: `--debug`, `--mesh_plot`.

---

## 9. Extensibility points

- **Python material modules**: `material.k = "z3st.materials.ceramic.k"` or `material.Gc = "z3st.materials.oxide.Gc"` — any importable callable that returns a UFL expression.
- **Custom constitutive laws**: `material.constitutive: custom` + `material.stress_function: "pkg.mod.func"`; signature `f(u, T, material, model) -> σ`.
- **Custom crystal plasticity**: add `get_cp_internal_variables(u, T, material, model)` next to the stress function.
- **Planned integrations** (README roadmap): `dolfinx_mpc` for multi-point constraints, `dolfinx_materials` for standardised material libraries, **Merope** (microstructure generation), rate-theory solvers, Monte Carlo workflows, contact mechanics.

---

## 9bis. Currently missing capabilities

The following capabilities are **not present** in Z3ST v0.1.0 and would need to be implemented to reproduce polycrystalline-RVE studies such as Aydiner et al. (2024) on dual-phase steels:

- **Cohesive zone model (e.g. Park–Paulino–Roesler / PPR)** — no `models/cohesive_model.py` exists. Intergranular decohesion (F/F, F/M) requires zero-thickness interface elements with a traction–separation law, which FEniCSx does not support natively as ABAQUS UEL does. Would need either a custom implementation (mesh duplication along grain boundaries + paired surface elements with a traction–separation potential) or a reformulation via discontinuous Galerkin / penalty-based surface terms.

- **Uncoupled ductile damage indicators (e.g. Bao–Wierzbicki, Modified Mohr–Coulomb)** — no triaxiality-dependent accumulated-damage model is implemented. Would be a straightforward addition on the existing quadrature spaces (`Q_pl`), driven by the plasticity state, as a post-processing indicator without two-way coupling to the stress response.

- **Multi-point constraints / periodic BCs with master-node coupling** — `dolfinx_mpc` is listed as a planned integration but is not yet wired into the BC infrastructure. This is required for enforcing constant stress triaxiality on an RVE (e.g. Eq. 22 of Aydiner et al.) and for true periodic boundary conditions on polycrystalline cells.

---

## 9c. Active work-in-progress: case 14 thermal-shock fracture (UO2)

The case-14 family is being calibrated against McClenny et al., JNM 565 (2022) 153719, in preparation for figures in the companion paper at `~/research-manuscripts/z3st_paper/`. State as of 2026-05-15.

**Three case-14 variants, each with a distinct role (verdict reached 2026-05-07 after extensive iteration and an independent NotebookLM consult):**

| Case directory | Role | Why |
|---|---|---|
| `14_full_cylinder_cracking/` (3D) | **Gold-standard** McClenny reproducer. | Faithfully captures the 60° azimuthal wedge AND the radial-only heat transfer (top/bottom faces zero-flux, matching McClenny's alumina-spacer design). |
| `14_full_cylinder_cracking_2D_xy/` (2D plane strain) | **McClenny Fig. 8 reproducer** at lower compute cost. | Plane strain (no axial gradients) is consistent with the alumina spacer's role; the 60° contact arc on a transverse cross-section is exactly McClenny's 2D representation (their Fig. 8 top, Fig. A.13, Fig. A.14). Modeled as upper-half disc with mirror symmetry on y=0 (= 30° contact in the upper half). |
| `14_full_cylinder_thermal_2D_rz/` (2D axisymmetric) | **Verification only** (thermal + linear-elastic). Damage disabled. | Axisymmetric mode mathematically prohibits azimuthal variation, so it cannot represent the 60° contact wedge. Any axisymmetric idealization (full-circumference cooling, single-face quench, etc.) either contradicts McClenny's experimental design or produces an unphysical annular damage band. Useful only to verify z3st's axisymmetric thermal solver against the analytic Bessel-series solution. |

The **alumina spacer detail** is critical and easy to overlook: McClenny p.3 notes "Insulation is placed on one side of the fuel pellet to ... eliminate axial thermal contact between the bottom of the capsule and the UO2 so that conductive radial heat transfer was the primary method of heat transfer to occur. This was intentionally designed to form a stress concentration on the contact region to induce fracturing." This means the experiment is actively radial-only by design — plane strain (2D-xy) is the *correct* dimensional reduction, and any axisymmetric variant with axial gradients (e.g. cooling only the top face) would *contradict* the experiment.

**Methodological framing (the paper's contribution).** Z3ST's damage block implements the **Ambati et al. (2015) hybrid (isotropic-anisotropic) phase-field formulation** (Comput Mech 55:383-405, Eq. 27). McClenny et al. instead use the **Miehe anisotropic formulation with viscous Allen-Cahn evolution** (their Eq. 10, with viscosity `eta = 1e-8 s/mm`). The case-14 chapter is therefore not "same problem, same model, different code"; it is a benchmark that the hybrid model — whose mechanical block is **linear** and which has no viscosity-tuned kinetics — captures the same crack topology at a per-iteration cost roughly an order of magnitude lower (Ambati §3.1). The Ambati paper is checked into `z3st/cases/14_full_cylinder_cracking/`.

**Implementation correspondence:** verified against the Ambati paper on 2026-05-05.
- Eq. (27a) `sigma = (1-D)^2 dPsi0/de` ↔ `damage_model.py:30-37` `g(D) = (1-D)^2 + K`, applied to the full stress in the linear mechanical block of `solver.py::_mechanical_step`.
- Eq. (27b) `-l^2 Lap d + d = (2l/Gc)(1-d) H+` ↔ `solver.py:528-530` AT2 weak form (`(H+1) u v + l^2 grad u . grad v = H v`), with `H = (2l/Gc) Psi+` from `damage_model.py:54` and irreversibility `H = max(H_old, H_new)` at line 252.
- Eq. (27c) `Psi+ < Psi- => d := 0` ↔ `damage_model.py:233-245`, **softened**: Z3ST sets `H -> 0` in compression cells rather than `D -> 0`. Equivalent under monotonic loading (the thermal-shock case here); more physical than the literal Ambati projection under cyclic loading because it preserves accumulated damage. This deviation is intentional and documented in the docstring of `update_history`.

**Calibration choice (uo2.yaml, 2026-05-15):**
- The material card declares `sigma_c: 1.0e+9 Pa` directly; `Gc` is auto-derived in `spine.py:147-149` (AT1: `Gc = (8/3)·lc·σc²/E ≈ 372 J/m²` at `lc = 5e-5 m`). The AT1 elastic threshold is `ψc = σc²/(2E) = 3·Gc/(16·lc) ≈ 1.4 MJ/m³`, below the cold-rim tensile-hoop-stress driving energy (order 4–5 MJ/m³ at the peak rim hoop-stress), so damage initiates and propagates along the contact arc. Note: McClenny et al.'s macro-tuned effective `Gc ≈ 80 kJ/m²` (their Table 3) is **not directly reachable in strict AT1** because of the algebraic identity above — matching their `Gc` with `σc = 1 GPa` would require `lc ~ 0.48 m`, far larger than the pellet radius. The present calibration (`σc = 1 GPa`, `lc = 50 μm`, `Gc ≈ 372 J/m²`) is the AT1 sweet spot under the rigid identity, and the AT1+Ambati-hybrid framing is a methodological alternative to McClenny's Miehe-AT2 + viscous-Allen-Cahn, **not** a direct reproduction of their effective `Gc`. (Earlier docs and uo2.yaml comment carried `Gc ≈ 93 J/m²` and `ψc ≈ 0.35 MJ/m³`, off by 4×; root cause: σc=500 MPa was inadvertently used in the arithmetic, which gives σc²/4. Corrected here on 2026-05-18.)
- An earlier `σc = 2 GPa` calibration (Gc ≈ 1490 J/m², ψc ≈ 5.6 MJ/m³) caused the prescribed seed to stall after ~50–100 μm of propagation: the singular tip concentration was the only place that crossed the threshold, and as soon as the tip cells softened (D → 1), the stress relaxation **shielded** the next-shell cells back below threshold. Diagnostic (probing the actual tip at r = 9.75 mm, θ = 15°): ψ⁺ in the next-shell ring was 0.07–0.19 MJ/m³ vs ψc = 5.6 MJ/m³. Halving σc → 1 GPa drops the threshold below the bulk tensile-rim driving force everywhere on the contact arc, breaking the shielding lock.
- The earlier worry that `σc = 1 GPa` would trigger spurious bulk damage from the plane-strain Amor-split deviatoric channel is **resolved**: the regime-aware-eigenstrain fix in `damage_model.py::_thermal_eigenstrain` (eth_zz = 0 in 2D Cartesian) suppresses the bulk pollution. The deciding diagnostic: bulk T (r < 0.5R) is uniformly at T_init = 1023 K (no cooling penetration yet at t = 0.1 s); 0 bulk nodes have D > 0.1 outside the contact wedge. Only 2 nodes at the (-R, 0) Clamp_x pin show D > 0.1 (a known FE point-constraint stress-concentration artifact, spatially contained to a single element, harmless).

**Mesh refinement (mesh.geo, 2026-05-15):**
- `lc_outer = 12.5 μm` (h/lc = 0.25, the Borden/Miehe ideal for AT1 bandwidth resolution); previously 25 μm (= lc/2) was borderline and contributed to the shielding lock.
- `Field[2].DistMax = 5.0e-3` (down from 8.0e-3) — keeps the fine zone confined to the expected crack-growth corridor; element count grows ~15% (15589 → 17965 nodes), not 4× as the lc_outer halving alone would suggest.

**Two structural code fixes applied 2026-05-06/07 (regime-aware elastic strain in the damage driver):**
1. The damage driving force `psi_pos` is evaluated on the *elastic* strain `eps_el = eps(u) − α·(T−T_ref)·I`, not on the total strain `eps(u)`. Without this, uniform thermal expansion in the bulk (where the body is unconstrained) produces a spurious `psi_pos` that drives damage everywhere, leading to a runaway `E_el` cascade. Fixed in `damage_model.py::psi_amor_split`, `psi_miehe_spectral`, `psi_split`, `crack_driving_force`, `update_history`; `solver.py::solve_staggered` now passes `T=T_new` to `update_history`.
2. The thermal eigenstrain z-component is *suppressed* in 2D Cartesian regimes (`regime: 2d` or `plane_stress`). Otherwise the plane-strain constraint `eps_zz = 0` combined with `eth_zz = α·ΔT` produces `eps_el_zz = −α·ΔT`, whose deviatoric component generates a uniform bulk `psi_pos ≈ (2/3)·G·α²·ΔT² ≈ 5.6 MJ/m³` at our parameters — exactly at the AT1 threshold for `σc = 2 GPa`, causing immediate divergence. Fixed in `damage_model.py::_thermal_eigenstrain` and mirrored in `mechanical_model.py::elastic_energy_density`. Axisymmetric and 3D regimes are unaffected (their `eps_zz` is dynamic).
3. `compute_energy_balance` now (a) uses the elastic strain for `E_el` (consistent with (1)), (b) applies the regime weight `w = 2π·r` for axisymmetric integrals.

**What still needs running (the user runs these locally — do not invoke from this assistant):**
```bash
# 3D (gold standard)
cd ~/z3st/z3st/cases/14_full_cylinder_cracking/
./Allrun

# 2D-xy (McClenny Fig. 8 reproducer)
cd ~/z3st/z3st/cases/14_full_cylinder_cracking_2D_xy/
./Allrun

# 2D-rz (verification only, damage off)
cd ~/z3st/z3st/cases/14_full_cylinder_thermal_2D_rz/
./Allrun
```

`Allrun` chains `Allclean → gmsh → python -m z3st > log_z3st.md → non-regression.py` in each case.

**Run results (2D-xy, 2026-05-15, n_steps=100, dt=1 ms, wall-clock 311 s):**
- E_el: 236 → 2372 J (down 19% vs the pre-calibration σc=2 GPa run, because cracks now release strain energy).
- E_frac: 0.56 → 3.74 J (up 52%).
- **~5 discrete radial cracks** within the 0–30° upper-half contact arc (= ~10 across the full disc by mirror symmetry), reproducing McClenny's Fig. 8 (top) "two major + fan of shorter cracks" topology. The pre-crack at θ=15° anchors the pattern; secondary cracks nucleate spontaneously at the elastic shielding length scale along the rest of the contact arc.
- Crack penetration depth: D = 1.0 from r ≈ 9.0 mm out to the rim (1 mm penetration; seed length is 250 μm so ~0.75 mm propagation beyond the seed).
- No spurious bulk damage: only 20 nodes with D > 0.1 outside the contact wedge (edge effect at the θ=30° BC discontinuity), and 2 nodes at the (-R, 0) Clamp_x pin. The bulk (r < 0.5R) is uniformly at T_init = 1023 K with D = 0.

**Methodological finding — secondary spontaneous nucleation (2026-05-15):**
With σc = 1 GPa, the AT1 threshold is crossed not just at the singular crack tip but across the entire cold contact arc once σ_θθ builds up. The first crack (anchored by the seed) sheds its rim load and pushes its immediate neighbors into a *compressive shadow* — those cells stay healthy. A few elements further along the arc, σ_θθ recovers above threshold and another crack nucleates. The result is the classical **Bažant/Bahr periodic thermal-shock cracking instability**: crack spacing settles to the elastic shielding length (≈ cold-front depth × O(1)). This is a stronger result than the previous (single-crack-from-seed) run and **invalidates the paper's prior claim** (main.tex line 365) that spontaneous nucleation requires a monolithic Newton-Krylov coupling — the staggered scheme with the regime-aware-eigenstrain fix and AT1 threshold handles it cleanly. The seed is still kept (i) for numerical robustness at step 0 (avoids the cascade as the rim cells transition abruptly from D=0 to D~1), and (ii) for CI reproducibility (anchors the pattern to a known location; without the seed, crack locations would drift between runs/machines).

**Expected outcomes:**
- 3D case: same physics, in 3D. Crack initiation around `t ≈ 10⁻² s` per McClenny p.8. Crack bands will appear wider than McClenny's because of the lc coarsening; the topology and timing are the diagnostic targets. The same σc=1 GPa / lc=50 μm calibration will likely apply; verify by running.
- 2D-rz case: damage is OFF. Reports thermal radial profile vs analytic Bessel series (target: <1% L2 error), thermo-elastic stress profiles, and energy balance. Useful for verifying the axisymmetric integration weight `2π·r`.

**Notes from prior iterations (don't repeat the failed paths):**
- The 2026-05-07 A+B+C calibration plan (n_steps→500, stag_tol→1e-5, relax_D→0.3) was made *obsolete* by the actual solution: switching the mechanical and damage linear solvers to `direct_mumps` (the divergence was AMG terminating on garbage residuals on the 7-orders-of-magnitude heterogeneous SPD that AT1's H = 2·ψ⁺ produces near the rim), plus the prescribed pre-crack to stabilize step 0. The current input.yaml has all tols back to default; don't chase the A+B+C path again.
- The 2026-05-07 fallback plan to add a viscous Allen–Cahn term `η ∂D/∂t` to AT1 is *not needed*. The current run completes 100 steps in 311 s with no divergence anywhere.
- The earlier suspicion that the adaptive controller silently bypasses `relax_min` was **investigated and ruled out**: `solver.py:623` enforces the floor correctly. The previous run observing `relax_D = 0.05` with `relax_min = 0.3` was a configuration artefact (user-supplied initial value below the floor), not a framework bug.
- The two SENT/SENS benchmarks in `cases/19_single-edge_notched_*` are correctly configured per Miehe 2010 / Ambati §4 and serve as standalone empirical verification of the hybrid-model implementation.

**Paper integration point:** `~/research-manuscripts/z3st_paper/main.tex` has subsection `\ref{sec:case14}` (lines 311–321) with the case-14 narrative and Figs. `case14_thermal_shock_results.png` and `case14_damage_field.png` / `case14_hoop_stress.png`. The chapter narrative foregrounds (i) the Ambati-hybrid-vs-Miehe-anisotropic formulation difference, (ii) the alumina-spacer rationale for plane-strain dimensional reduction, (iii) the σc=1 GPa calibration choice and the σc-vs-Gc-vs-lc coupling, and (iv) the secondary-nucleation finding. The figures in `figures/` and the numerical claims at line 317 and line 321 should match the values listed in **Run results** above.

---

## 10. File-level quick index

| File                                          | Role                                                  |
|-----------------------------------------------|-------------------------------------------------------|
| `z3st/__main__.py`                            | CLI entry point, time-stepping loop                   |
| `z3st/core/spine.py`                          | Top-level `Spine` driver (multi-inheritance)          |
| `z3st/core/config.py`                         | Parses `input.yaml` into `self.on`, paths, regime     |
| `z3st/core/finite_element_setup.py`           | Allocates `V_t, V_m, V_d, V_c, V_pl, Q`               |
| `z3st/core/solver.py`                         | Staggered solver, PETSc options, DG cluster solver    |
| `z3st/core/mesh/{reader,manager,plotter}.py`  | Mesh IO, tag management, PyVista preview              |
| `z3st/core/diagnostic.py`                     | Structured logger                                     |
| `z3st/models/thermal_model.py`                | Thermal BCs + heat-flux diagnostics                   |
| `z3st/models/mechanical_model.py`             | Strain/stress tensors, constitutive routes, mech BCs  |
| `z3st/models/damage_model.py`                 | AT1/AT2 phase-field, energy splits, history update    |
| `z3st/models/plasticity_model.py`             | J2 return mapping + custom crystal-plasticity hook    |
| `z3st/models/gap_model.py`                    | Fixed / Gas gap-conductance model (+ contact-coupled) |
| `z3st/models/contact_model.py`                | Penalty pellet-clad mechanical contact (PCMI)         |
| `z3st/models/creep_model.py`                  | Implicit Norton creep (incremental variational, IFT tangent by AD) |
| `z3st/models/cluster_dynamic_model.py`        | 1D advection–diffusion cluster dynamics (DG/SIPG)     |
| `z3st/materials/*.yaml`                       | Material cards                                        |
| `z3st/materials/{ceramic,oxide}.py`           | Python callables for `k(T)`, `Gc(mesh)`                |
| `z3st/materials/{fuel_profiles,fuel_swelling}.py` | Fuel-behaviour callables: radial power `f(r,bu)`, swelling(bu) |
| `z3st/utils/writer.py`                        | Unified `OutputWriter` (VTU / single-file XDMF; merged fields) |
| `z3st/utils/utils_load.py`                    | YAML loader + power-history generator                 |
| `z3st/utils/utils_plot.py`, `plot_convergence.py` | Plotting helpers                                  |
| `z3st/utils/utils_extract_{vtu,xdmf}.py`      | Field extraction from output files                    |
| `z3st/utils/utils_verification.py`            | Analytical benchmarks                                  |
| `z3st/utils/z-gui.py`, `interactive_gui.ipynb`| Interactive viewers                                   |
| `z3st/utils/geo_files/*.geo`                  | Reusable Gmsh templates                                |
| `z3st/cases/…`                                | ~40 verification / validation / demo cases            |
| `docs/source/*.rst`                           | Sphinx user & API documentation                        |

---

*Generated on 2026-04-16 for Z3ST v0.1.0; last updated 2026-06-10 (full four-agent re-audit folded into punch_list.md; CODE-P0-5 plane_stress solver fix + regime validation; gold-regression verdict now persisted to `non-regression.json` and gated in local summary + CI; 8 stale golds re-blessed; `heat_flux` fixed and re-wired under `--debug`; `U_pwr_rod_2D` gold-protected and wired into the local suite together with the `V_*` cases).*