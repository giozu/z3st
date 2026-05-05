# Z3ST — Repository Context

**Z3ST** is an open-source **FEniCSx-based finite-element framework** for coupled thermo-mechanical material analysis, written in Python. It is developed by **Giovanni Zullo** (Politecnico di Milano), licensed under **Apache 2.0**, version **0.1.0 (2025)**.

- **Repository:** https://github.com/giozu/z3st
- **DOI:** 10.5281/zenodo.17748028
- **Citation:** `CITATION.cff` + BibTeX entry in README
- **Continuous integration:** GitHub Actions (`.github/workflows/ci.yml`, `static.yml`)
- **Docs:** Sphinx sources under `docs/source/`, built by GitHub Actions
- **Python:** ≥ 3.10 (pyproject.toml)
- **Runtime dependencies:** `numpy`, `scipy`, `pandas`, `matplotlib`, `pyvista≥0.42`, `meshio`, `pyyaml`, `gmsh` + (external) **dolfinx / FEniCSx / basix / UFL / PETSc / MPI** provided by the Conda env `z3st_env.yml`.

---

## 1. Purpose and scope

Z3ST provides a **modular and extensible** environment to simulate:

- steady-state and transient **heat conduction** in arbitrary multi-material domains;
- **linear / hyperelastic / elasto-plastic** mechanical response;
- **phase-field fracture** (damage) with AT1 and AT2 models;
- **gap conductance** between domains (fixed or gas-type);
- **1D cluster-dynamics** advection–diffusion problems;
- arbitrary, spatially-dependent **internal heat sources** (γ-heating, user-defined);
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
    │   ├── gap_model.py              Fixed / Gas gap conductance
    │   └── cluster_dynamic_model.py  1D advection-diffusion (DG+SIPG+upwind)
    ├── materials/                    YAML material cards + Python modules
    │   ├── steel.yaml, austenitic_steel.yaml, martensitic_steel.yaml,
    │   │ high_carbon_steel.yaml, T91.yaml, 15_15Ti.yaml,
    │   │ vessel_steel.yaml, vessel_steel_0.yaml
    │   ├── uo2.yaml, zircaloy.yaml
    │   ├── ceramic.yaml, oxide.yaml, plastic.yaml, lead.yaml, h2o.yaml
    │   └── ceramic.py, oxide.py                    Python-side callables
    │                                               (e.g. k(T), Gc(mesh))
    ├── utils/                        post-processing + helpers
    │   ├── export_vtu.py             VTU writer for T, u, D, σ, ε
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
6. `set_power()` — build `q_third` (W/m³) per material: fissile (LHR/area), γ-heating with exponential decay in rect geometry or `K₀(μr)/K₀(μRᵢ)` in cylindrical geometry, spherical decay, etc.
7. `solve(dt, max_iters)` — dispatch to `solve_staggered`.
8. `get_results()` — build symbolic UFL strain/stress/stress_mech/stress_th/energy_density dictionaries per material.
9. `compute_energy_balance()` (damage only) — assembles elastic + fracture energy (`E_el`, `E_frac`) for diagnostics written to `energies.txt`.

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
        ├── parameters(lhr); set_power()
        ├── solve(max_iters, dt)           # staggered
        ├── get_results()                  # symbolic σ, ε, ψ
        ├── compute_energy_balance         # if damage
        └── export_vtu  or  XDMF           # output/fields_*.vtu / .xdmf
```

CLI flags:
- `--debug`       enable verbose debugging
- `--mesh_plot`   preview surface tags with PyVista (via `core.mesh.plotter.MeshPlotter`)

---

## 3. Core modules

### 3.1 `core/config.py — Config`

Parses the user YAML and fills:

- `self.on`: `{thermal, mechanical, damage, cluster, plasticity}` flags from `models:`
- `self.gap_model`, `self.h_gap_value` from `models.gap_conductance`
- paths to geometry, mesh, boundary conditions
- `self.n_steps`
- `self.regime ∈ {2d, 3d, axisymmetric, plane_stress}`

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
- `heat_flux(T)` diagnostic: per-material average `|q|`, `q_x`, `q_y`, `q_z`.
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

Thermal stress coupling: `σ_th = −(3λ + 2G) α (T − T_ref) I` subtracted from the mechanical right-hand side.

Damage coupling: when damage is active, `σ ← g(D) σ` with `g(D) = (1−D)² + K` (K = 1e-6 regularization).

### 4.3 Damage (`damage_model.py`)

Phase-field fracture with two variational models:

- **AT2**: quadratic local dissipation `w(D)=D²`, no elastic threshold.
  `H = (2 lc / Gc) · ψ⁺`  (non-dimensional history).
  Miehe spectral split (positive eigenvalues of `ε`).
- **AT1**: linear local dissipation `w(D)=D`, analytical elastic threshold `σ_c`.
  `H = ψ⁺` (physical, J/m³).
  Amor (volumetric/deviatoric) split.

Elastic-energy splits:
- `psi_miehe_spectral(u, mat)` — 2D closed-form + 3D via Cardano's formula with smooth clamping.
- `psi_amor_split(u, mat)` — `ψ⁺ = ½ λ ⟨tr ε⟩₊² + G dev(ε):dev(ε)`; `ψ⁻ = ½ λ ⟨tr ε⟩₋²`.

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

### 4.6 Cluster dynamics (`cluster_dynamic_model.py`)

1D advection–diffusion solver for defect-cluster size distributions `c(n,t)`:
- `∂c/∂t = −v ∂c/∂n + D ∂²c/∂n²`
- initial conditions: `constant` (on labelled region) or `gaussian`.
- DG1 space with upwind advection and SIPG diffusion.
- Mass conservation: `∫ c·n dn` is rescaled to the initial target every step.

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
| `fissile`           | bool       | If true, `q''' = LHR / area` in the pellet              |
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
├── mesh.msh      generated mesh (checked in for CI)
├── non-regression.py         post-run comparison vs. output/non-regression_gold.json
└── output/                   auto-generated VTU/XDMF + plots
```

The suite is driven by `z3st/cases/non-regression.sh` (local) and `non-regression_github.sh` (CI) and summarised in `non-regression_summary.txt`.

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
- `14_full_cylinder_cracking` (3D), `14_full_cylinder_cracking_2D_rz` (axisymmetric UO₂ quench with AT2 phase-field)

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

**U_* — Extended/user cases**
- `U_box_knotch_3D`, `U_cluster_dynamics_test`, `U_coaxial_contact_2D`, `U_pressure_vessel_2D`, `U_quarter_block`, `U_slab_contact`, `U_spherical_shell`, `U_thick_cylindrical_shell_plane_stress`

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

time:  [0.0, 0.01]
lhr:   [0.0, 0.0]
n_steps: 20

output:
  format: xdmf                 # xdmf | vtu
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
| Python material callables         | ✓ `k(T)`, `Gc(mesh)` loaded via `importlib`                         |
| Post-processing                   | VTU / XDMF writers, radial/1D plots, PyVista viewer, notebook GUI   |
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
python -m z3st > log.z3st
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
| `z3st/models/gap_model.py`                    | Fixed / Gas gap-conductance model                     |
| `z3st/models/cluster_dynamic_model.py`        | 1D advection–diffusion cluster dynamics (DG/SIPG)     |
| `z3st/materials/*.yaml`                       | Material cards                                        |
| `z3st/materials/{ceramic,oxide}.py`           | Python callables for `k(T)`, `Gc(mesh)`                |
| `z3st/utils/export_vtu.py`                    | VTU writer                                            |
| `z3st/utils/utils_load.py`                    | YAML loader + power-history generator                 |
| `z3st/utils/utils_plot.py`, `plot_convergence.py` | Plotting helpers                                  |
| `z3st/utils/utils_extract_{vtu,xdmf}.py`      | Field extraction from output files                    |
| `z3st/utils/utils_verification.py`            | Analytical benchmarks                                  |
| `z3st/utils/z-gui.py`, `interactive_gui.ipynb`| Interactive viewers                                   |
| `z3st/utils/geo_files/*.geo`                  | Reusable Gmsh templates                                |
| `z3st/cases/…`                                | ~40 verification / validation / demo cases            |
| `docs/source/*.rst`                           | Sphinx user & API documentation                        |

---

*Generated on 2026-04-16 for Z3ST v0.1.0.*