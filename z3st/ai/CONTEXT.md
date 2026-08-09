# Z3ST — Repository Context

**Z3ST** is an open-source **FEniCSx-based finite-element framework** for coupled thermo-mechanical material analysis, written in Python. It is developed by **Giovanni Zullo** (Politecnico di Milano), licensed under **Apache 2.0**, version **0.3.0 (2026)**.

- **Repository:** https://github.com/giozu/z3st
- **DOI:** 10.5281/zenodo.17748028
- **Citation:** `CITATION.cff` + BibTeX entry in README
- **Continuous integration:** GitHub Actions (`.github/workflows/ci.yml`, `static.yml`)
- **Docs:** Sphinx sources under `docs/source/`, built by GitHub Actions
- **Python:** ≥ 3.10 (pyproject.toml)
- **Runtime dependencies:** `numpy≥2`, `scipy`, `matplotlib`, `pyvista≥0.42`, `pyyaml`, `gmsh`, `h5py`, `shapely` (extras: `post` = `pandas`; `docs`; `dev`; `nn` = `torch`, `dolfinx-external-operator`) + (external) **dolfinx / FEniCSx / basix / UFL / PETSc / MPI** provided by the Conda env `z3st_env.yml`. The `deps` audit check verifies this against `pyproject.toml`.

> **Students:** onboarding, working agreements, and per-student topic/branch tracking
> live in the private note `../personal_stuff/z3st_students.md` (not in this repo).

---

## 1. Purpose and scope

Z3ST provides a **modular and extensible** environment to simulate:

- steady-state and transient **heat conduction** in arbitrary multi-material domains;
- **linear / hyperelastic / elasto-plastic** mechanical response;
- **phase-field fracture** (damage) with AT1 and AT2 models;
- **gap conductance** between domains (fixed or gas-type);
- **1D cluster-dynamics** advection–diffusion problems;
- **thermal-gradient pore migration** (oxide-fuel restructuring: central void +
  columnar zone) via a stabilised porosity-advection model two-way coupled to
  the thermal solve (§4.10);
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
├── pyproject.toml                    installable Python package (PEP 621)
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
    │   ├── solver.py                 staggered loop + services (577 lines; the
    │   │                             physics steps live in models/, see §3.3)
    │   ├── spine.py                  Spine — top-level driver (multi-inheritance of mixins)
    │   └── mesh/                     mesh sub-package
    │       ├── reader.py             Gmsh → dolfinx loader
    │       ├── manager.py            MeshManager (tags, geometry metrics, normals)
    │       └── plotter.py            PyVista mesh/tag preview
    ├── models/                       physical model mixins (plugged into Spine)
    │   ├── thermal_model.py          heat conduction + Dirichlet/Neumann/Robin BCs
    │   ├── mechanical_model.py       linear Lamé, hyperelastic, plasticity, custom
    │   ├── damage_model.py           phase-field AT1/AT2, Miehe/Amor splits
    │   ├── plasticity_model.py       J2 (return mapping) + custom crystal plasticity hook
    │   ├── gap_model.py              Fixed / Gas gap conductance (+ contact-coupled)
    │   ├── contact_model.py          penalty pellet-clad mechanical contact (gap closure)
    │   ├── cluster_dynamic_model.py  1D advection-diffusion (DG+SIPG+upwind)
    │   ├── creep_model.py            implicit Norton (+irradiation) creep, AD tangent
    │   ├── cracking_model.py         Barani isotropic-softening fuel cracking
    │   ├── porosity_migration_model.py  thermal-gradient pore migration (SU/SUPG)
    │   └── nn_conductivity.py        neural-network k(T) (Picard / external-operator Newton)
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
    │   ├── logger.py                 framework-wide logger (`log`)
    │   ├── plot_convergence.py       staggered-residual plots
    │   ├── utils_extract_vtu.py      scalar/vector/tensor extraction from VTU
    │   ├── utils_extract_xdmf.py     same for XDMF
    │   ├── utils_load.py             YAML loader + power-history generator
    │   ├── utils_plot.py             1D/radial plots (T(r), σ_rr(r), …)
    │   ├── utils_verification.py     reference analytical benchmarks
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
        ├── non-regression_local.sh    discovery-based suite (Allrun + gold)
        ├── non-regression_github.sh   CI suite, reads cases_ci.txt
        ├── cases_ci.txt, suite_exclude.txt
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

**stdout filters** (`__main__.py::_install_stdout_filters`, top of file). One wrapper, installed unconditionally on every MPI rank, carrying two independent per-write flags: `markdown` and `emit`. **Both must stay unconditional** — the wrapper cannot be installed on some ranks and not others, because dolfinx's PETSc wrappers do collective work in `__del__` and a rank holding a different set of live objects has the collector order those destructors differently, deadlocking the run at exit. And it cannot be gated on `isatty()`, because under `mpirun` stdout *is* a tty even when the shell redirects to a file. **Nothing automated guards this.** Every case in both suites runs serially, so a regression here is invisible to them; parallel runs are exercised by hand.

`emit` is false on every rank but 0, so the ~260 diagnostic prints are not repeated once per process; `[WARNING]` and `[ERROR]` lines pass from all ranks. `markdown` is false under `Z3ST_PLAIN_LOG=1` or on a tty; when true (i.e., the user redirected output, e.g. `python -m z3st > log.md`), a line-by-line stdout shim rewrites the existing decorated output into Markdown:
- Morse-code dividers `--.. ..- .-.. .-.. --- ...` → `***` (horizontal rule; `***` rather than `---` to avoid setext-heading ambiguity).
- `[STEP NN/MM] rest` → `## Step NN/MM: rest`.
- `--- Staggering iteration N/M ---` → `#### Iteration N/M`.
- `--. foo --..` (spine section markers) → `### foo`.
- `__Foo__` (initializer headers) → `### Foo`.
- `[DESCRIPTION]` → `## Description`.
- `[INFO|WARNING|ERROR|SUCCESS] body` → `**[TAG]** body`.

Interactive runs and `Z3ST_PLAIN_LOG=1` are pass-through for the markdown transform, but the rank gate still applies. None of the solver / model / config modules were touched; the filter is a single point of intercept in `__main__.py`.

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
- `self.regime ∈ {1d, 2d, 3d, axisymmetric}` (`1d` is used by the teaching cases)

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
| `V_p`  | Lagrange P1 (scalar)                   | Porosity (when `porosity` on)           |

### 3.3 `core/solver.py — Solver`

Implements the staggered loop and the services the physics mixins consume.

**Where the steps live.** `solver.py` went from 1579 lines to
577: each physics mixin now owns its own staggered step, and the solver provides
services rather than containing the physics. The `_*_step` subsections below are
still the reference for what each step does; only their home moved.

| step | file |
|---|---|
| `_thermal_step`, `_thermal_step_nonlinear` | `models/thermal_model.py` (with `_thermal_conductivity_aux_operands`, `_build_gap_pair_aux`, `_refresh_gap_pair`) |
| `_mechanical_step` | `models/mechanical_model.py` |
| `_damage_step` | `models/damage_model.py` |
| `_cluster_step` | `models/cluster_dynamic_model.py` |
| `_porosity_step` | `models/porosity_migration_model.py` (was always there) |

Still in `solver.py`: `solve_staggered`, `get_solver_options`, `_stagger_residual`,
`_adapt_relax`, `_bc_objects`, `_value_at_step`, `_build_measures`, `_global_tags`,
`invalidate_dt_caches`, and at module level `as_bool`, `aitken_omega` and the three
rigid-body nullspace builders. `Spine` multiply-inherits everything, so every call
site is unchanged and `self._thermal_step(...)` still resolves.

**Aitken Δ² dynamic relaxation**: the displacement relaxation factor is recomputed every staggered iteration from the last two raw residuals R_k = ũ_k − u_k as ω_{k+1} = −ω_k·(R_{k−1}·ΔR)/|ΔR|², clamped to [relax_min, relax_max] — the standard cure for slowly-contracting interface coupling (Küttler & Wall). Supersedes the heuristic grow/shrink controller for u (thermal keeps the EMA controller). Two safeguards, both load-bearing: ω restarts from the configured `relax_u` at every time step (the recursion scales each new ω from the previous one, so a clamped-at-the-floor ω from a degenerate step — e.g. the zero-power initial step — would otherwise poison every later estimate), and the update is skipped when ‖ΔR‖ < 1e-8·‖R‖ (converged/zero-load steps produce a garbage quotient). Motivation: the PCMI contact phase of `regression/pwr_rod_2D` ran ~130-190 staggered iterations/step with the heuristic controller collapsed to relax_min.

**Per-step form caching**: `_thermal_step` and `_mechanical_step` build their UFL forms and Linear/Nonlinear problem objects **once per time step** (cache keyed on `current_step` and the solution-Function identities) instead of once per staggered iteration; `LinearProblem.solve()` reassembles A and b from the stored forms, and every iteration-varying input is consumed by reference — the gap conductance is a persistent `Constant` owned by the gap model, the paired-surface temperatures persistent Functions refreshed per iteration, the contact pressure was already a `Constant`, and creep predictor/state, burnup, and T live in Functions. Per-step rebuild is forced anyway by dt (creep) and the cracking rescale, which both bake floats. Effect: a 109-step PWR-rod run drops from ~2 500 UFL form constructions to ~109. Validated bit-for-bit against the gold suite (7 cases PASS/PASS spanning thermal cache, gap-pair refresh, contact, SNES creep).

**Shared step helpers**. Four physics steps had written
out the same staggered bookkeeping:
- `_stagger_residual(new, old, cfg, tol, label)` — scatter both, `copy`, `axpy(-1)`,
  the two norms, the guarded division, and the verdict. Returns
  `(converged, norm_d, rel_norm_d, residual)`. **The printed strings are parsed by
  `utils/plot_convergence.py`**, which greps the solver log for exactly
  `||ΔX||/||X|| = <float>`: changing the label or the spacing empties every
  convergence plot, with nothing to catch it.
- `_adapt_relax(name, residual, prev)` — the EMA grow/shrink controller, previously
  written three times (T, u, D) identical bar the attribute name.
- `aitken_omega(r_k, r_prev, omega, comm, n_owned, lo, hi)` — module level, shared
  with `porosity_migration_model`. The caller owns the clamp and the starting ω
  because the two uses genuinely differ: the displacement loop clamps to
  `[relax_min, relax_max]` and carries ω across time steps, porosity clamps to a
  hard-coded `[0.05, 1]` and restarts from `aitken_omega0`.
- `_regime_normal()` lives on `MechanicalModel` and is called from the solver's
  traction update, which had a copy of the dispatch.

The `coupling` key is gone: all cases said `staggered`, the other
branch raised, so `solve()` calls `solve_staggered` directly.

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

### 3.5 `utils/logger.py`

Framework-wide `logging` logger named `z3st` (`log.info`, `log.warning`, …), reused
across modules. Two things about it are load-bearing :

- **It writes to stdout, not stderr.** The case `Allrun` redirects stdout only, so on
  stderr its output never reached `log_z3st.md` — the mesh diagnostics were absent
  from the very file CI dumps when a case fails. The handler targets stdout through a
  late-bound proxy, so it does not matter whether `__main__` wrapped the stream before
  this module was imported. The `[LEVEL] message` format is what the markdown filter
  already renders as `**[LEVEL]** message`.
- **A rank filter**: INFO from rank 0 only, WARNING and ERROR from every rank. It is
  the only gate when z3st is imported as a library; under `python -m z3st` the stdout
  wrapper gates the same way (see §2.3).

It keeps its own handler with `propagate = False`. An earlier `logging.basicConfig`
here configured the *root* logger, which silently set the level and format for
matplotlib, h5py and anything else that logs.

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
| `hyperelastic`  | Neo-Hookean `ψ = μ/2(I_C − 3) − μ ln J + λ/2 (ln J)²`, σ = (1/J) P F⊤ |
| `plasticity`    | J2 return-mapping with linear isotropic hardening                     |
| `custom`        | Loads `material.stress_function = "pkg.mod.func"`, `σ = f(u, T, material, model)` |

Regime handling in `epsilon(u)`:
- **axisymmetric:** strain in cylindrical coords `(r, θ, z)`; `ε_rr = ∂u_r/∂r`, `ε_θθ = u_r/r`, `ε_zz = ∂u_z/∂z`, `ε_rz = ½(∂u_r/∂z + ∂u_z/∂r)`.
- **2D:** 3×3 tensor with zero z-components (plane strain).
- **3D:** `sym(∇u)`.

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
- `psi_amor_split(u, mat)` — `ψ⁺ = ½ K_n ⟨tr ε⟩₊² + G dev(ε):dev(ε)`; `ψ⁻ = ½ K_n ⟨tr ε⟩₋²`, with `K_n = λ + 2G/n` Amor's n-dimensional bulk modulus (n = strain-tensor dimension) so ψ⁺+ψ⁻ = ψ_el exactly.
- `psi_star_convex(u, mat)` — Vicentini, Zolesi, Carrara, Maurini, De Lorenzis 2024 (Int. J. Fract. 247:291-317). One-parameter generalisation of Amor controlled by `dmg_cfg["gamma_star"]` (a *model* parameter in the `damage:` block, default 0 → reduces to Amor; intentionally not a per-material property). `ψ⁺ = G|dev ε|² + (K_n/2)[⟨tr ε⟩₊² − γ⋆⟨tr ε⟩₋²]`, `ψ⁻ = (1+γ⋆)(K_n/2)⟨tr ε⟩₋²` (same `K_n` as Amor, matching Vicentini's κ). Satisfies all five criteria in the Vicentini 2024 Table 2 (the other splits do not). `γ⋆ > 0` raises the compressive-vs-tensile critical-stress ratio.

`update_history(u)` — vectorised, per-material update of the history field `H` on DG0. Supports **Ambati-Gerasimov-De Lorenzis hybrid constraint** (`dmg_cfg.hybrid_constraint`, default `True`): where `ψ⁻ > ψ⁺` locally, contribution to H is set to 0 to suppress crack growth in compression.

`compute_energy_balance(u)` returns `(E_el, E_frac)`:
- AT2: `γ = Gc/2 · (D²/lc + lc |∇D|²)`
- AT1: `γ = 3 Gc/8 · (D /lc + lc |∇D|²)`

BCs: only `Dirichlet` on D (e.g. `D=0` on healthy boundary).

### 4.4 Plasticity (`plasticity_model.py`)

- **J2 plasticity** with isotropic linear hardening; constitutive update computed symbolically through `ufl.conditional` (returns stress `σ = σ_tr − 3μ Δp n`, yield stress `σ_y = σ_y0 + H p`).
- History fields `ep`, `p`, `ep_n`, `p_n` on quadrature spaces.
- `mode: custom` in `plasticity:` block hooks a user crystal-plasticity function `get_cp_internal_variables` inside the material's module — used by the `verification/plasticity/crystal_single_grain` case.
- `update_plastic_history(u)` refreshes `ep_n`, `p_n` after the staggered iteration converges.

### 4.5 Gap conductance (`gap_model.py`)

Two modes (selected via `models.gap_conductance.type`):

- **Fixed** — constant `h_gap` (W/m²K) read from `value:`.
- **Gas** — `k_gas = value · 1e-4 · T_gap^0.79`, then `h_gap = k_gas / gap_size` where `gap_size` is computed as the mean distance between two paired labelled facet groups (via SciPy cKDTree on facet centroids), and `T_gap = ½ (T_inner + T_outer)`.

Invoked inside `_thermal_step` when a Robin BC is defined with `pair:` to another subdomain.

**Contact-coupled conductance.** When `gap_conductance.contact_coupling.enabled` is set, a solid-contact term is added on gap closure (Todreas & Kazimi, *Nuclear Systems I*, 3rd ed., Eqs. 8.141/8.142): the emergent contact pressure (from `contact_model`) raises `h_gap` above the open-gap gas value, so closing the gap cools the fuel. Parameters: `meyer_hardness` (Pa), `gas_thickness` (m, roughness-based residual gas space). The Ross-Stoute harmonic mean accepts symbolic k(T) cards by evaluating them at the current mean gap temperature (`_k_at_gap`; UFL folds constants, so `k_func(float)` is a plain number) — previously a symbolic fuel conductivity silently zeroed `h_contact` and disabled the contact-cooling feedback.

**Conductance under-relaxation** (`models.gap_conductance.relax`, default 1.0 = off): h is damped between staggered iterations, `h ← ω·h_new + (1−ω)·h_prev`, with the memory reset every time step. The contact-pressure → conductance → temperature → expansion → pressure feedback is the loop that chatters on gap closure; damping h attacks it at the source. `h_gap` is returned as a persistent `Constant` (updated in place) so the cached thermal form needs no rebuild.

### 4.6 Creep (`creep_model.py`)

Implicit Norton creep (`ε̇_eq = A0·exp(−Q/RT)·σ_eq^n`) for a material carrying `creep: norton` + `creep_A0/n/Q` on its card — the dissipative extension of the energy-first design (incremental variational principle, Ortiz–Stainier). The cell-local minimisation over Δε_cr condenses to the scalar radial-return equation per point; a DG0 **predictor** Δγ₀ holds its exact root (vectorised numpy Newton, refreshed before every mechanical solve, consistency gated in the staggered convergence test), and the UFL stress carries **one symbolic Newton step** from the predictor — so `ufl.derivative` yields exactly the implicit-function-theorem consistent tangent through a trivially small expression tree (a fully unrolled symbolic Newton explodes FFCx). The accumulated `ε_cr` is a per-material DG0 tensor state advanced once per converged step; it enters the trial through the eigenstrain channel and the output stress via `creep_output_stress`. Mechanical steps auto-promote to SNES when creep is active. v1 scope: isotropic Lamé, no damage/plasticity combination, regimes with 3×3 strain tensors. Verified by `verification/fuel/creep` (1e-14) and `verification/fuel/creep_relaxation` (4e-15 vs the BE recursion; O(dt) defect pinned).

**Irradiation creep**: a linear-in-stress in-pile term ε̇_irr = B·φ·σ_eq is added inside the same radial return, g(Δγ) = Δγ − Δt·[A(T)·base^n + B·φ·base] = 0 with base = σ_eq − 3GΔγ — the linear term preserves the monotone/concave structure, so the Newton convergence argument is unchanged. Card keys `creep_irr_B` (Pa⁻¹ per n/m²) and `fast_flux` (n/(m²·s)), both required together (spine validates); absent → exact no-op (both creep verification golds unchanged at machine precision). First user: the Zircaloy clad of `regression/pwr_rod_2D` (B = 2.0e-36, φ = 7.0e17 → ε̇ ≈ 1e-10 s⁻¹ at 80 MPa, ~1 % over three years — in-pile creep-down at 580 K, where thermal Norton alone is negligible).

### 4.7 Cluster dynamics (`cluster_dynamic_model.py`)

1D advection–diffusion solver for defect-cluster size distributions `c(n,t)`:
- `∂c/∂t = −v ∂c/∂n + D ∂²c/∂n²`
- initial conditions: `constant` (on labelled region) or `gaussian`.
- DG1 space with upwind advection and SIPG diffusion.
- Mass conservation: `∫ c·n dn` is rescaled to the initial target every step.

### 4.8 Mechanical contact (`contact_model.py`)

Penalty pellet-clad mechanical contact (gap closure / PCMI), enabled via the `models.contact` block. An explicit fixed-point scheme integrated into the staggered loop:

- **Gap measurement** — each mechanical iteration measures the current normal gap from the displacement iterate as `gap = g0 + ⟨u_r⟩_b − ⟨u_r⟩_a`, the boundary-integral mean radial displacement of the two paired facing surfaces (`surface_a` = pellet outer, `surface_b` = clad inner).
- **Penalty traction** — on penetration (`gap < 0`) a pressure `p = k_pen · ⟨−gap⟩₊` is applied as `t = −p·n` on both facing surfaces (UFL/AD supplies the tangent). Config: `penalty_stiffness` (Pa/m), `initial_gap` (m).
- Verified against the analytical plane-stress Lamé interference-fit solution (`cases/verification/fuel/shrink_fit`, 1.0 %). The emergent contact pressure also feeds the contact-coupled gap conductance (§4.5), so thermal + mechanical PCMI are two-way coupled.

### 4.9 Fuel cracking — isotropic softening (`cracking_model.py`)

The Barani et al. (NED 342, 2019) isotropic-softening model for fuel cracking: the cracked pellet is represented by globally rescaled elastic constants, conserving principal strains and minimising the squared principal-stress deviation between the cracked (anisotropic) and equivalent isotropic descriptions. Scaling from the VIRGIN constants (kept as `E_virgin`/`nu_virgin` on the card):

```
f(ν)     = (2/3)·(2−ν)/(2+ν)·1/(1−ν)
E_iso(n) = f(ν)^n · E            ν_iso(n) = ν / (2^n + (2^n − 1)·ν)
```

The number of cracks follows the paper's empirical correlation on the rod-average LHR, n = n₀ + (n∞ − n₀)(1 − exp(−(LHR − LHR₀)/τ)) above LHR₀, with the fitted constants LHR₀ = 5 kW/m, n₀ = 1, n∞ = 12, τ = 21 kW/m (Oguma 1983 / Walton & Husser 1983 data; all overridable via `cracking_lhr0/n0/n_inf/tau`). No healing: n is driven by the maximum LHR seen in the history (irreversible, `_lhr_max` on the card). Opt-in per material card with `cracking: isotropic` (unknown values rejected at load like `creep:`); the rescale runs once per step from `spine.parameters()`, and since the mechanical form rebuilds per step (form cache), the softened lmbda/G are consumed with no extra plumbing. At 20 kW/m: n ≈ 6.6, E_iso/E ≈ 0.11, ν_iso ≈ 0.003 — order-of-magnitude lower fuel stresses, the paper's headline effect. Unit-checked against the paper (n = 1 at 5 kW/m with E_iso/E = f(ν); Fig. 3 curve at 10/20/40 kW/m; irreversibility). Scope: elastic softening only — Jankus-Weeks cracked-fuel creep correction and healing deliberately excluded, as in the paper.

### 4.10 Porosity migration (`porosity_migration_model.py`)

Thermal-gradient pore migration for oxide-fuel
restructuring (central void + columnar zone), in the spirit of Barani et al.
(JNM 2022) / Novascone et al. (2018). The porosity field `p(x,t)` (P1 on `V_p`)
obeys the conservative advection law `∂p/∂t + ∇·(v p) = 0` with the pore
velocity **up the thermal gradient**:

```
v = v0 · (c1 + c2·T + c3·T² + c4·T³) · T^-2.5 · exp(-Hs/RT) · ∇T   [m/s]
```

(Sens vapour-transport correlation; `v0 = 1.303427e8` reproduces Sens'
4.2e-11 m/s benchmark; `Hs = 5.98e5` J/mol; all overridable in the
`porosity:` block.)

- **Discretisation:** CG P1 + streamline stabilisation, backward Euler,
  point-wise clamp to `[0,1]` after each solve. `porosity.stabilisation:
  su` (default, artificial diffusion K = (h/2|v|) v⊗v) or `supg` (consistent:
  test perturbed by τ v·∇w against the full strong residual, τ =
  ((2/dt)² + (2|v|/h)²)^-1/2 — sharper front). Note the clamp costs strict
  mass conservation (~2.5 % on the reference transient); the conservative
  DG+limiter variant lives on the `feature/DG` research branch.
- **Boundary:** natural (zero-flux) by default — the rim keeps its fabrication
  porosity as a *prediction*. `porosity.rim_inflow_porosity` (+ optional
  `porosity.rim_label`, default `outer`) switches to a weak inflow/outflow
  split (prescribed value only where `v·n < 0`).
- **Two-way thermal coupling:** `q'''` is rescaled by `1/(1-p)` (LHR
  conserved) and a material with `thermal_conductivity_model: kato_porosity`
  gets `k(T,p)` = Kato matrix conductivity × Maxwell-Eucken porosity factor,
  refreshed each staggered iteration (`update_porosity_dependent_properties`).
- **Staggered contract:** `_porosity_step(p_new, p_n, dt, ...)` keeps `p_n`
  frozen at tⁿ across staggered iterations (the pattern the cluster solver now
  follows too); under-relaxation fixed factor `porosity.relax` or Aitken
  (`porosity.aitken: true`); convergence on the mixed rel/abs criterion.
- **State/rollback:** `porosity`/`porosity_n` are in `_SNAPSHOT_FIELDS`, so the
  adaptive grid rolls them back on a failed step.
- **Verification:** `verification/fuel/porosity_migration` (Barani (U,Pu)O2
  sector: centre p → 1, void radius 0.18 r/Ro vs 0.20, rim 0.15 preserved,
  centre T 2924 K gold-tracked).

---

## 5. Material database (`z3st/materials/`)

Materials are plain YAML cards. Common fields:

| Key                 | Units      | Meaning                                                |
|---------------------|------------|--------------------------------------------------------|
| `name`              | —          | Human-readable name                                    |
| `E`, `nu`           | Pa, —      | Young's modulus, Poisson's ratio                        |
| `k`                 | W/(m·K) or `"mod.func"` | Thermal conductivity (scalar OR symbolic; `materials.fuel_thermal.k` = Fink 2000 UO2 k(T), 95 % TD) |
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
| `swelling_rate`     | (MWd/kgU)⁻¹ | ΔV/V per unit burnup for `fuel_swelling.solid_gas_densification` |
| `gas_swelling_rate` / `gas_T_onset` / `gas_T_width` | (MWd/kgU)⁻¹, K, K | gaseous-swelling amplitude and thermal-activation sigmoid of `fuel_swelling.solid_gas_densification` (defaults 4.0e-4, 1200, 150) |
| `densification_dv` / `densification_bu` | —, MWd/kgU | in-pile densification amplitude and burnup constant of `solid_gas_densification` (defaults 0.010, 2.0); set `densification_dv: 0` to switch the term off |
| `cracking`          | —          | `isotropic` enables the isotropic-softening fuel-cracking model (§4.9); `cracking_lhr0/n0/n_inf/tau` override the correlation constants |
| `creep_irr_B` / `fast_flux` | Pa⁻¹ per n/m², n/(m²·s) | irradiation-creep coefficient and fast flux for the optional ε̇ = B·φ·σ term (§4.6); both or neither |
| `sigma_c` / `Gc`    | Pa, J/m²   | Phase-field critical stress OR fracture energy (one is derived from the other given `lc`) |
| `yield_strength`    | Pa         | Initial yield stress (J2 plasticity)                    |
| `hardening_modulus` | Pa         | Linear isotropic hardening modulus                      |
| `constitutive`      | string     | `lame`, `hyperelastic`, `plasticity`, `custom` |
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

The suite is driven by `z3st/cases/non-regression_local.sh` (local) and `non-regression_github.sh` (CI) and summarised in `non-regression_summary.txt`. The local suite is discovery-based: every directory under `cases/` with both an `Allrun` and a blessed `output/non-regression_gold.json` is a member (`sandbox/` is never scanned); exceptions live in `cases/suite_exclude.txt` with a reason per line, and `--list` prints the discovered set. A case is protected if and only if it has a gold. The CI list is curated separately in `cases/cases_ci.txt`, which `non-regression_github.sh` reads. The is chosen for **coverage** against a stated time budget (22 cases, 14 min 47 s, from measured per-case times) rather than purely for speed; its header records the two models no case reaches at all, so nobody looks for them there. Each case's `non-regression.json` carries two verdicts: `"summary"` (analytic-tolerance check) and `"regression"` (vs the blessed `non-regression_gold.json`); the local summary reports both per case, and CI fails when either is FAIL.

### 6.1 Catalogue of cases

**00 — Tutorial / example**
- `verification/mechanics/uniaxial_tension` — uniaxial traction of a rectangular steel block (3D, linear elasticity, `Neumann` traction + `Clamp` BCs on 3 faces).

**1–6 — Slabs and thin cylindrical shells (thermal / thermo-mechanical benchmarks)**
- `verification/thermal/thin_slab_dirichlet_2D`, `verification/thermal/thin_slab_neumann_2D`, `verification/thermal/thin_slab_neumann_3D`, `verification/mechanics/uniaxial_tension_nonlinear`
- `verification/thermal/thin_cylindrical_shell_dirichlet_2D`
- `verification/thermal/thick_slab_adiabatic_3D`, `verification/thermal/thin_slab_adiabatic_2D`
- `verification/thermal/thick_cylindrical_shell_adiabatic_2D`, `verification/thermal/thin_cylindrical_shell_adiabatic_2D`
- `verification/thermal/thick_slab_non_adiabatic_3D`, `verification/thermal/thin_slab_non_adiabatic_2D`
- `verification/thermal/thick_cylindrical_shell_non_adiabatic_2D`, `verification/thermal/thin_cylindrical_shell_non_adiabatic_2D`

**7 — Box heated**
- `verification/thermal/box_heated`

**8–9 — Thick cylinders with regime variations**
- `verification/mechanics/lame_plane_strain_2D`
- `verification/mechanics/lame_gps_2D`, `verification/mechanics/lame_gps_3D` (generalized plane strain)

**11–14 — Cylinders (Mariotte, gradients, annular, full cylinder, fracture)**
- `verification/mechanics/mariotte_thin_shell`
- `verification/mechanics/thermal_gradient_2D`, `verification/mechanics/thermal_gradient_3D`
- `verification/mechanics/annular_cylinder`
- `verification/mechanics/full_cylinder`
- `benchmarks/damage/pellet_quench_2D_xy` (plane-strain McClenny Fig. 8 reproducer — the primary case for the paper's case-14 chapter). See §11 for the variant rationale. Axisymmetric transient-cooling verification is exercised by no case.

**15 — Cavities and pressurised bodies**
- `verification/mechanics/elliptical_cavity_2D`, `regression/two_elliptical_cavities_2D`
- `verification/mechanics/spherical_cavity`

**16 — Multi-body coupling**
- `verification/thermal/coaxial_gap_3D`

**17 — Stress–strain curves**
- `verification/mechanics/stress_strain_displacement`, `verification/mechanics/stress_strain_stress`
- `verification/damage/double_crack_2D`, `verification/damage/notched_plate_2D`

**18 — 2D fracture benchmarks**
- `verification/damage/box_crack_2D`, `verification/damage/box_notch_2D`

**19 — Single-edge notched (classical phase-field benchmarks)**
- `benchmarks/damage/sen_shear`
- `benchmarks/damage/sen_tension`

**20 — Plasticity**
- `verification/plasticity/j2_hardening_2D`

**I, II — Utilities**
- `studies/mesh_sensitivity_2D` — mesh convergence study
- `studies/attenuation_map` — γ attenuation in materials. A non-suite study (custom `run_map.py`/`plot_map.py`, no `Allrun`/`non-regression.py`).

**V_* — Analytical-verification cases** (closed-form checks)
- `verification/fuel/swelling` — constant volumetric swelling eigenstrain (free expansion → σ ≈ 0, exact `u`).
- `verification/fuel/fuel_swelling` — burnup-driven swelling reading the `burnup` field (the eigenstrain bus consuming the state bus).
- `verification/fuel/burnup` — burnup accumulation + radial-power source bus on an axisymmetric pellet (closed-form mean burnup; rim/core ratio = 1 + A).
- `verification/fuel/axial_power` — axial-power source bus (chopped-cosine `f(z)`, Todreas & Kazimi 1-D axial problem) on a tall axisymmetric fuel column: closed-form mean burnup (machine precision); axial peaking factor = 1/[(2L′/πL)·sin(πL/2L′)]; end/peak = cos(πL/2L′).
- `verification/fuel/axial_table` — tabulated axial profile (`tabulated_axial`, piecewise-linear node-wise peaking factors — the standard core-physics input): closed-form mean burnup (exact); table-node ratio f₃/f₁ (machine precision); peak/mean = max f / trapezoid mean.
- `verification/fuel/creep` — implicit Norton creep, constant-stress uniaxial bar (backward Euler exact): total/creep/radial strain vs closed form at 1e-14 (radial pins the deviatoric −½ flow).
- `verification/fuel/creep_relaxation` — stress relaxation at held strain: Z3ST ≡ the scalar backward-Euler recursion at 4e-15; deviation from the exact `σ(t)` equals the predicted O(dt) defect (2.11% at 50 steps, pinned to 2e-13).
- `verification/fuel/creep_law_discovery` — the framework's first **inverse / constitutive-identification** case (EUCLID-style, independent implementation — the published EUCLID codes are GPL-3.0 and are not used). Forward problem: the relaxation case re-run with **500** implicit steps (own `input.yaml`; ~36 min — data defect ~0.2%, below the noise; observations cached to `output/fem_stress_history.csv`). Inverse: 51 noisy (2%) observations of mean σ_zz; 5-mechanism library `{S, S², S³, S⁵, sinh S}` (S = σ/σ_ref); self-contained dual-number forward AD through the implicit BE integrator (different grid from the data → no inverse crime); damped Gauss-Newton on log-coefficients; backward elimination by strain share + **one-standard-error rule** (plain BIC could NOT separate {S²,S³} from {S³}; with the original coarse 50-step data a spurious S² absorbed the time-discretisation defect in 6/10 seeds — hence the 500-step forward run). Result: cubic Norton selected alone **10/10 seeds**, coefficient to 1.6%. Gold-blessed (`discover.py` ~16 s given the CSV); feeds the z3st paper's identification section (figure `creep_law_discovery.png`). In the local suite only, not in CI: the forward run is too long.
- `verification/fuel/shrink_fit` — penalty contact pressure vs the analytical plane-stress Lamé interference-fit. The drive is a uniform Dirichlet ramp (300 → 1500 K, lhr = 0): a uniform pellet temperature expands stress-free, so the Lamé reference is exact rather than approximate. Siblings `shrink_fit_disk` (2D plane strain) and `shrink_fit_disk_3d`.

**U_* — Extended / demo cases**
- `regression/pwr_rod_2D` — generic-PWR fuel-rod segment (4.5 mm pellet, 65 µm cold gap, Zircaloy clad), the framework's integral fuel-performance case. Physics: Fink (2000) UO2 k(T) (`materials/fuel_thermal.py`), Robin coolant film (h = 3.5e4 W/(m²·K), T = 580 K), burnup + rim-peaking radial power, solid + gaseous swelling with early-life densification (`fuel_swelling.solid_gas_densification`), Barani isotropic-softening fuel cracking (§4.9: n ≈ 6.6 at power, E_iso/E ≈ 0.11), Zircaloy thermal Norton + irradiation creep (§4.6), gas gap conductance with contact coupling (+ relax 0.5 damping), penalty contact, 15.5 MPa coolant and 2 MPa He fill-gas pressures. History: ramp to 20 kW/m in 20 d, 1800 d hold; weighted time grid `n_steps: [8, 60, 40]` — fine through the gap-closure window at ~330 d, strided across the creep plateau. Solver: MUMPS both blocks, Aitken Δ² relaxation, mech stag_tol 5e-4. Gold state: T_max peaks at 1210.8 K at day 65, PCMI onset at day 315 / 9.95 MWd/kgU average burnup, the contact-conductance feedback cools the pellet to 1149.7 K (−61 K on closure), contact pressure plateaus at 22.35 MPa; end of life at 1800 d with 58.23 MWd/kgU average (139.1 peak rim) and gap −0.45 µm. Mean burnup matches the closed form to machine precision. Wall-clock ~105 min for 109 steps. Gold-protected (end-state PCMI scalars from `output/history.csv`, mean burnup against the closed form), excluded from the routine local suite via `cases/suite_exclude.txt`; not in CI.
- `verification/thermal/spherical_shell` — gold-protected (semi-analytic checks).
- `U_pressure_vessel_2D` (`cases/sandbox/`): its `non-regression.py` only extracts CSV/plots (no asserts, never writes `non-regression.json`), and the `non-regression_gold.json` on disk is orphaned — it holds Lamé-style L2 errors the current script cannot produce (likely inherited from a deleted case).
- `U_cluster_dynamics_test`, `U_quarter_block` — unvalidated sandboxes under `cases/sandbox/`.

**verification/plasticity/crystal_single_grain** — crystal-plasticity single-grain demo using the `custom` constitutive + `plasticity.mode: custom` hook.

### 6.2 Typical `input.yaml` structure

```yaml
mesh_path: mesh.msh
geometry_path: geometry.yaml
boundary_conditions_path: boundary_conditions.yaml

materials:
  uo2: ../../materials/uo2.yaml

regime: axisymmetric            # 1d | 2d | 3d | axisymmetric

solver_settings:
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
| Linear elasticity                 | ✓ isotropic Lamé (1d / 2d plane strain / 3d / axisymmetric)         |
| Anisotropic elasticity            | ✗ not implemented |
| Hyperelasticity                   | ✓ Neo-Hookean (SNES Newton)                                         |
| Thermo-mechanical coupling        | ✓ staggered; adaptive (EMA grow/shrink) or Aitken Δ² dynamic relaxation (`relax_aitken`); per-step form caching; gap-conductance under-relaxation (`gap_conductance.relax`) |
| Phase-field fracture (AT1, AT2)   | ✓ Miehe/Amor split, hybrid constraint, irreversibility              |
| J2 plasticity                     | ✓ return mapping + linear isotropic hardening (quadrature elements) |
| Crystal plasticity                | experimental — via `custom` constitutive hook (`verification/plasticity/crystal_single_grain`) |
| Gap conductance                   | ✓ Fixed or Gas (k_gas = f(T_gap), gap_size from facet centroids)    |
| Cluster dynamics (1D)             | ✓ DG upwind + SIPG, mass-conservation renormalisation                |
| Axisymmetric / 2D / 3D / plane-stress regimes | ✓ all consistent with the integration weight `w`      |
| Volumetric heating                | ✓ fissile (LHR/area), γ-heating (rect / cyl / sphere analytic decay), user `q'''` |
| Burnup accumulation               | ✓ per-fissile-material `burnup` field via `update_state(dt)` (state bus)            |
| Radial power shaping              | ✓ `radial_profile` form factor `f(r, bu)` (source bus); built-in rim-peaking        |
| Axial power shaping               | ✓ `axial_profile` form factor `f(z)` (source bus, composed `f_r·f_z`, single mean-1 normalisation); built-ins: chopped cosine (T&K), tabulated (node-wise peaking factors) |
| Cladding creep (implicit, AD)     | ✓ Norton + Arrhenius via the incremental variational principle (`models/creep_model.py`): condensed radial return on the displacement space, DG0 predictor + one symbolic Newton step → exact IFT consistent tangent by `ufl.derivative`; per-material `ε_cr` DG0 state; card keys `creep: norton`, `creep_A0/n/Q`; optional irradiation creep ε̇ = B·φ·σ in the same radial return (`creep_irr_B` + `fast_flux`); verified to 1e-14 (constant stress) and 4e-15 vs the BE recursion (relaxation) |
| Fuel cracking (isotropic softening) | ✓ Barani et al. (2019) model (`models/cracking_model.py`): n(LHR_max) macro-cracks rescale E and ν from the virgin constants; irreversible; card `cracking: isotropic` (§4.9) |
| Constitutive-law identification   | ✓ EUCLID-style sparse mechanism selection from simulation data (`cases/verification/fuel/creep_law_discovery/discover.py`): candidate-library fit via self-contained forward-mode AD (dual numbers) through the implicit BE integrator, Gauss-Newton + backward elimination + one-SE rule; cubic Norton recovered 10/10 noise seeds; independent of the GPL-3.0 EUCLID codes |
| Integrated-power diagnostic       | ✓ `set_power` prints the exact FE integral of the fissile source per material per step (regime-weighted, MPI-reduced); note the mean-1 normalisation is *nodal*, so a radially peaked profile integrates to LHR·Lz·⟨f⟩_area/⟨f⟩_nodal (= 1.2·LHR·Lz for rim-peaking A=3, p=8) — pinned by the `total_power` checks in the burnup-family `V_` cases |
| Fuel swelling                     | ✓ constant ΔV/V or burnup-driven eigenstrain (eigenstrain bus); `solid_gas_densification` adds T-activated gaseous swelling + early-life densification |
| Time stepping                     | ✓ piecewise-linear power history; `n_steps` as an int (duration-proportional) or per-segment interval list (decouples resolution from segment length) |
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

cd z3st/cases/verification/mechanics/uniaxial_tension/
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

## 10. Missing capabilities

The following capabilities are **not present** in Z3ST v0.1.0 and would need to be implemented to reproduce polycrystalline-RVE studies such as Aydiner et al. (2024) on dual-phase steels:

- **Cohesive zone model (e.g. Park–Paulino–Roesler / PPR)** — no `models/cohesive_model.py` exists. Intergranular decohesion (F/F, F/M) requires zero-thickness interface elements with a traction–separation law, which FEniCSx does not support natively as ABAQUS UEL does. Would need either a custom implementation (mesh duplication along grain boundaries + paired surface elements with a traction–separation potential) or a reformulation via discontinuous Galerkin / penalty-based surface terms.

- **Uncoupled ductile damage indicators (e.g. Bao–Wierzbicki, Modified Mohr–Coulomb)** — no triaxiality-dependent accumulated-damage model is implemented. Would be a straightforward addition on the existing quadrature spaces (`Q_pl`), driven by the plasticity state, as a post-processing indicator without two-way coupling to the stress response.

- **Multi-point constraints / periodic BCs with master-node coupling** — `dolfinx_mpc` is listed as a planned integration but is not yet wired into the BC infrastructure. This is required for enforcing constant stress triaxiality on an RVE (e.g. Eq. 22 of Aydiner et al.) and for true periodic boundary conditions on polycrystalline cells.

### 10.1 Cluster dynamics: what the case does and does not pin

`verification/cluster/mass_conservation_1D` (9 s, in `cases_ci.txt`) exercises
`models/cluster_dynamic_model.py`. It pins the total defect mass, the per-step
rescale factor and the final peak density against a gold. All five metrics are
`tracked()`: they have no closed form at this configuration, so the gold is the
only thing constraining them.

**What the case verifies is that the rescale works, not that the scheme conserves
mass.** The rescale holds the first moment `C1 = integral of c*n dn` constant, but
the PDE solved is `dc/dt = -v dc/dn + D d2c/dn2`, whose first moment obeys

    dC1/dt = v*C0 + D*c(n_min),     C0 = integral of c dn

With `v > 0` clusters grow, so `C1` grows: it is not an invariant of this
equation. The per-step gap before rescaling is `dt/n_bar` with `n_bar = C1/C0`
(1.9-3.6 % at `n_bar ~ 1.3`, falling as the peak advances). It is not a
discretisation error -- refining the mesh tenfold leaves it unchanged, and moving
the initial peak to `n = 20` drops it to 0.26 %, matching `dt/n_bar` exactly.

The open question is a modelling one: if the constraint is physically right,
because clusters exchange monomers among themselves, then the PDE is missing the
free-monomer sink that must accompany growth, and the rescale stands in for it.

Note that the mesh starts at `n = 1`, not 0, and the upwind scheme imposes no
inflow there.

**A methodological warning.** `contact_model.py` was
initially listed here as having zero cases. It has seven, three of them with golds,

and `verification/fuel/shrink_fit_disk` is in CI. The error: contact is configured
as a *block* —

```yaml
models:
  contact:
    surface_a: lateral_1
    penalty_stiffness: 5.0e13
```

— and `Config` stores `bool(models.get("contact", False))`, which is `True` for any
non-empty dict. Grepping for `contact: true` finds none of them. When auditing which
cases reach a model, parse the yaml and evaluate the switch the way `Config` does;
a grep for a syntactic form answers a different question.

---

## 11. Case 14 — thermal-shock fracture (UO2)

`benchmarks/damage/pellet_quench_2D_xy` is the plane-strain reproducer of
McClenny et al., JNM 565 (2022) 153719, Fig. 8, and the case behind the paper's
case-14 chapter. Plane strain carries the 60 degree azimuthal contact wedge that an
axisymmetric idealisation cannot represent.

**Calibration (`uo2.yaml`).** The card sets `sigma_c: 1.0e+9` Pa; `Gc` is derived
in `spine.py` (AT1: `Gc = (8/3)*lc*sigma_c^2/E ~ 372 J/m2` at `lc = 5e-5 m`). The
AT1 elastic threshold is `psi_c = sigma_c^2/(2E) = 3*Gc/(16*lc) ~ 1.4 MJ/m3`, below
the cold-rim tensile-hoop driving energy (4-5 MJ/m3 at the peak rim hoop stress),
so damage initiates and propagates along the contact arc.

McClenny's macro-tuned effective `Gc ~ 80 kJ/m2` (their Table 3) is not reachable
in strict AT1: the identity above would require `lc ~ 0.48 m`, larger than the
pellet radius. This calibration is the best-conditioned choice under that identity,
and the AT1 + Ambati hybrid is a methodological alternative to McClenny's
Miehe-AT2 + viscous Allen-Cahn, not a reproduction of their effective `Gc`.

**Mesh.** `mesh.geo` refines the rim, where the hoop stress peaks and the damage
band localises.

**Two properties the damage driver depends on:**

1. `psi_pos` is evaluated on the elastic strain `eps_el = eps(u) - alpha*(T-T_ref)*I`,
   not on the total strain. On the total strain, uniform thermal expansion in the
   unconstrained bulk drives damage everywhere.
2. `compute_energy_balance` uses the elastic strain for `E_el`, and applies the
   regime weight `w = 2*pi*r` for axisymmetric integrals.

**Solvers.** The mechanical and damage blocks use `direct_mumps`. AMG terminates on
garbage residuals against the seven-order-of-magnitude heterogeneous SPD that AT1's
`H = 2*psi_plus` produces near the rim. A prescribed pre-crack stabilises step 0.

## 12. File-level quick index

| File                                          | Role                                                  |
|-----------------------------------------------|-------------------------------------------------------|
| `z3st/__main__.py`                            | CLI entry point, time-stepping loop                   |
| `z3st/core/spine.py`                          | Top-level `Spine` driver (multi-inheritance)          |
| `z3st/core/config.py`                         | Parses `input.yaml` into `self.on`, paths, regime     |
| `z3st/core/finite_element_setup.py`           | Allocates `V_t, V_m, V_d, V_c, V_pl, Q`               |
| `z3st/core/solver.py`                         | Staggered solver, PETSc options, DG cluster solver    |
| `z3st/core/mesh/{reader,manager,plotter}.py`  | Mesh IO, tag management, PyVista preview              |
| `z3st/utils/logger.py`                       | Framework-wide logger (`log`)                         |
| `z3st/models/thermal_model.py`                | Thermal BCs + heat-flux diagnostics                   |
| `z3st/models/mechanical_model.py`             | Strain/stress tensors, constitutive routes, mech BCs  |
| `z3st/models/damage_model.py`                 | AT1/AT2 phase-field, energy splits, history update    |
| `z3st/models/plasticity_model.py`             | J2 return mapping + custom crystal-plasticity hook    |
| `z3st/models/gap_model.py`                    | Fixed / Gas gap-conductance model (+ contact-coupled) |
| `z3st/models/contact_model.py`                | Penalty pellet-clad mechanical contact (PCMI)         |
| `z3st/models/creep_model.py`                  | Implicit Norton creep (incremental variational, IFT tangent by AD) |
| `z3st/models/cluster_dynamic_model.py`        | 1D advection–diffusion cluster dynamics (DG/SIPG)     |
| `z3st/models/cracking_model.py`               | Barani isotropic-softening fuel cracking (§4.9)       |
| `z3st/models/porosity_migration_model.py`     | Thermal-gradient pore migration, SU/SUPG (§4.10)       |
| `z3st/models/nn_conductivity.py`              | Neural-network k(T): Picard or external-operator Newton |
| `z3st/materials/*.yaml`                       | Material cards                                        |
| `z3st/materials/{ceramic,oxide}.py`           | Python callables for `k(T)`, `Gc(mesh)`                |
| `z3st/materials/{fuel_profiles,fuel_swelling}.py` | Fuel-behaviour callables: radial power `f(r,bu)`, swelling(bu) |
| `z3st/utils/writer.py`                        | Unified `OutputWriter` (VTU / single-file XDMF; merged fields) |
| `z3st/utils/utils_load.py`                    | YAML loader + power-history generator                 |
| `z3st/utils/utils_plot.py`, `plot_convergence.py` | Plotting helpers                                  |
| `z3st/utils/utils_extract_{vtu,xdmf}.py`      | Field extraction from output files                    |
| `z3st/utils/utils_verification.py`            | Analytical benchmarks                                  |
| `z3st/utils/interactive_gui.ipynb`            | Interactive viewer (notebook)                         |
| `z3st/utils/geo_files/*.geo`                  | Reusable Gmsh templates                                |
| `z3st/cases/…`                                | ~40 verification / validation / demo cases            |
| `docs/source/*.rst`                           | Sphinx user & API documentation                        |

---

*Z3ST v0.3.0 — repository context.*

*This file describes the code as it stands, not how it got there. Keep it in step
in the same commit as the change, and keep it that way by **replacing and
deleting**, never by appending. A dated note or a paragraph about what something
used to be belongs in the git history.*
