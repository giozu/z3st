# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
[![DOI](https://zenodo.org/badge/648784453.svg)](https://doi.org/10.5281/zenodo.17748028)
[![CI](https://github.com/giozu/z3st/actions/workflows/ci.yml/badge.svg)](https://github.com/giozu/z3st/actions/workflows/ci.yml)
![static](https://github.com/giozu/z3st/actions/workflows/static.yml/badge.svg)
<!-- ![paper build](https://github.com/giozu/z3st/actions/workflows/paper.yml/badge.svg) -->

**Z3ST** is an open-source finite-element framework for thermo-mechanical material analysis, with nuclear fuel as its driving application.
Written in Python on FEniCSx, it couples heat conduction, small- and finite-strain mechanics, plasticity, creep, hybrid phase-field fracture, gap conductance and penalty contact in multi-material domains, under stationary or transient conditions and in 1D, 2D, 3D or axisymmetry.
Materials live outside the numerical core. The card keys `k`, `E`, `nu`, `Gc`, `eigenstrain`, `radial_profile`, `axial_profile` and `stress_function` accept the dotted path of an importable Python function. A function that returns a UFL expression enters the weak form symbolically, and its consistent tangent follows by automatic differentiation.

---

## Overview

Z3ST solves coupled thermo-mechanical problems with a staggered scheme on Gmsh meshes, with one material card per region and a piecewise-linear power history.
Each physics model is a mixin of the `Spine` class (see [`docs/source/architecture.md`](docs/source/architecture.md)). External codes are coupled through the same eigenstrain and source terms the internal models use (SCIANTIX, [`z3st/coupling/sciantix/`](z3st/coupling/sciantix/README.md)).

---

## Quick installation

Clone the repository:
  ```bash
  git clone https://github.com/giozu/z3st.git
  ```

Z3ST requires FEniCSx (dolfinx), which is not installed by pip. Install it with Conda or use a Docker image.
To install miniconda:
  ```bash
  mkdir -p ~/miniconda3
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
  bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
  rm ~/miniconda3/miniconda.sh
  source ~/miniconda3/bin/activate
  conda init --all
  ```

Create and activate the Conda environment:
  ```bash
  cd z3st
  conda env create -f z3st_env.yml
  conda activate z3st
  ```

Install in editable mode:
  ```bash
  pip install -e .
  ```

Run a case from its directory, for instance:
  ```bash
  cd ~/z3st/z3st/cases/verification/mechanics/uniaxial_tension/
  gmsh mesh.geo -3 > log_mesh.md
  python3 -m z3st > log_z3st.md
  python3 non-regression.py
  python3 ../../../../utils/plot_convergence.py
  ```

`./Allrun` runs the same commands:
  ```bash
  cd ~/z3st/z3st/cases/verification/mechanics/uniaxial_tension/
  ./Allrun
  ```

Optional flags:

* `python3 -m z3st --mesh_plot` opens a PyVista window with the mesh and its facet tags before solving
* `python3 -m z3st --debug` prints every key of every material card after loading

---

# Simulation Cases

`z3st/cases` holds the simulation setups (*cases*) used for verification and demonstration.

Each case directory contains:
- YAML input configuration (`input.yaml`)
- Geometry and meshing definitions (`geometry.yaml`, `mesh.geo`)
- Boundary conditions (`boundary_conditions.yaml`)
- A non-regression test script (`non-regression.py`)
- Reference results in `output/non-regression_gold.json`

The cases are the regression tests and the worked examples.

---

## Key features

* **Coupled thermo-mechanical solver.** Heat conduction (stationary, or transient with backward Euler) and mechanics, coupled by a staggered scheme. Adaptive or Aitken Δ² relaxation, per-step form caching, optional gap-conductance damping.
* **Adaptive time-stepping.** Off by default. A step that does not converge is rolled back to the last converged state, its `dt` is halved, and it is solved as sub-steps. Output stays on the original grid. The run stops if a step fails at `dt_min`.
* **Hot-reloadable parameters.** Tolerances, relaxation factors and `max_iters` in `input.yaml` are re-read at the start of every step.
* **Regimes.** `1d`, `2d` plane strain, `3d` and `axisymmetric`, selected by one key. The axisymmetric weight `w = 2πr` and the hoop strain are built in.
* **Constitutive laws.** Small-strain isotropic Lamé, Neo-Hookean hyperelasticity (SNES Newton with line search), J2 plasticity with linear isotropic hardening, and a `custom` hook for a user UFL stress function (used by the crystal-plasticity demo).
* **Phase-field fracture.** AT1 and AT2 in the hybrid formulation of Ambati et al. (2015). The stress is fully degraded and the positive part of the energy drives the crack. Miehe spectral, Amor volumetric/deviatoric or star-convex split, irreversibility through a history field, and the hybrid constraint against damage in compression.
* **Penalty contact.** Between two concentric bodies separated by a uniform gap. One contact pressure per surface pair, from the surface-averaged gap, coupled to the gap conductance. Verified against the Lamé interference fit and, for creep relaxation, against an independent radial J2 solution (`reference_1d.py` in `z3st/cases/verification/fuel/creep_shrink_fit_2D`).
* **Porosity migration.** Temperature-gradient-driven pore advection. The porosity enters the thermal problem through a porosity-dependent conductivity and a rescaled volumetric source. Streamline-upwind (default), SUPG, or upwind DG-1 with SSP-RK3 and a vertex limiter. An SIPG diffusion block is optional and off by default.
* **Data-driven conductivity.** The Magni MA-MOX correlation, a Gaussian-process correction on its logarithmic residual, and a PyTorch neural network. Any object with the two-method `(k, dk/dT)` interface can be used, either in a lagged Picard iteration or in a Newton scheme through `dolfinx-external-operator`, evaluated at quadrature points.
* **SCIANTIX coupling.** Fission-gas behaviour computed per fuel point by SCIANTIX through the eigenstrain and source terms the internal models use. It returns gaseous swelling and fission gas release. Requires a compiled SCIANTIX shared library, see [`z3st/coupling/sciantix/`](z3st/coupling/sciantix/README.md).
* **Creep.** Implicit Norton creep with Arrhenius temperature dependence, through the incremental variational principle (radial return condensed onto the displacement space, consistent tangent by automatic differentiation). Optional flux-driven irradiation creep.
* **Fuel cracking.** Isotropic softening. The number of radial cracks follows the rod-average linear heat rate and rescales the fuel elastic constants, without recovery.
* **Fuel behaviour.** Burnup-driven solid and gaseous swelling with early-life densification as an eigenstrain. UO₂ k(T) in the modified NFI form of FRAPCON-3 at zero burnup (no burnup degradation). Rim-peaking radial and chopped-cosine or tabulated axial power profiles.
* **Multi-material domains.** Thermal, mechanical and damage properties per material, integrated with one measure per cell tag.
* **Volumetric heating.** Fissile materials receive LHR/area, optionally shaped by radial and axial profiles. Gamma heating with exponential decay in rectangular, cylindrical or spherical geometry.
* **Boundary conditions.** Thermal: Dirichlet, Neumann, Robin (convective or gap pair). Mechanical: Dirichlet (vector or one component), Neumann, Clamp, Slip. Damage: Dirichlet. Thermal Dirichlet and mechanical values can be lists with one value per time point.
* **Gap conductance.** `Fixed` or `Gas` (`k_gas = c·10⁻⁴·T_gap^0.79`) between paired facet groups. The gap width is the contact model's current gap when contact is on, otherwise the mean distance between the two surfaces from a SciPy cKDTree, recomputed at every call.
* **1D cluster dynamics.** Cluster size-distribution solver: implicit-Euler DG1 with upwind interior-facet flux and SIPG diffusion, with a rescaling that conserves the total cluster mass.
* **Python material functions.** Temperature-dependent `k(T)`, spatially graded `Gc(x)`, eigenstrains, power profiles and stress functions, imported with `importlib` at start-up.
* **Material cards.** YAML cards in `z3st/materials`: UO₂, MA-MOX, steels (austenitic, martensitic, high-carbon, T91, 15-15Ti, vessel), Zircaloy-4, ceramic, oxide, plastic, lead, H₂O. The values are representative, chosen for the demonstration and verification cases, and are not qualified design data. See [`z3st/materials/README.md`](z3st/materials/README.md).
* **Mesh input.** Gmsh `.msh` files. Gmsh templates under `z3st/utils/geo_files/`.
* **Configuration.** Three YAML files per case: `input.yaml`, `geometry.yaml`, `boundary_conditions.yaml`.
* **Parallel runs.** PETSc with MUMPS, GAMG or HYPRE BoomerAMG, MPI over `MPI.COMM_WORLD`. Neither test suite runs a case under MPI: every case in `non-regression_local.sh` and in `cases_ci.txt` runs on one rank, so a regression specific to parallel runs is not detected.
* **Output.** VTU or XDMF time series through one writer that compiles its interpolation expressions at setup. Readable by ParaView and PyVista.
* **Continuous integration.** GitHub Actions runs the cases listed in `cases_ci.txt` on every push to `main`, `develop` and `development/*` and on every pull request. Each case's `non-regression.py` compares its metrics with a version-controlled gold JSON.
* **Documentation.** Sphinx sources under `docs/source/`, built by GitHub Actions. UML class diagram in [`docs/source/architecture.md`](docs/source/architecture.md).

---

## Directory structure

```bash
z3st/                                # repository root
├── LICENSE.txt                      # Apache 2.0
├── README.md
├── CITATION.cff
├── pyproject.toml                   # installable package (PEP 621)
├── z3st_env.yml                     # Conda env recipe (FEniCSx + deps)
├── docs/                            # Sphinx documentation
│   ├── Makefile
│   └── source/
├── .github/workflows/
│   ├── ci.yml                       # non-regression CI
│   └── static.yml                   # Sphinx docs build
└── z3st/                            # Python package
    ├── __main__.py                  # CLI entry point
    ├── __init__.py                  # lazy-import facade
    ├── core/                        # FEM core
    │   ├── config.py                # YAML parser
    │   ├── spine.py                 # top-level Spine driver
    │   ├── solver.py                # staggered loop + services (physics steps
    │   │                           #   live in models/, one per mixin)
    │   ├── finite_element_setup.py  # V_t / V_m / V_d / V_c / V_pl / Q
    │   └── mesh/                    # Gmsh loader, MeshManager, PyVista preview
    │       ├── reader.py
    │       ├── manager.py
    │       └── plotter.py
    ├── models/                      # physics mixins plugged into Spine
    │   ├── thermal_model.py
    │   ├── mechanical_model.py      # lame / hyperelastic / plasticity / custom
    │   ├── damage_model.py          # AT1 / AT2, Miehe / Amor / star-convex splits, hybrid constraint
    │   ├── cohesive_model.py        # cohesive phase-field fracture (under development)
    │   ├── plasticity_model.py      # J2 + custom CP hook
    │   ├── creep_model.py           # implicit Norton + irradiation creep, AD tangent
    │   ├── cracking_model.py        # isotropic-softening fuel cracking
    │   ├── contact_model.py         # penalty pellet-clad contact (PCMI)
    │   ├── gap_model.py             # Fixed / Gas gap conductance + contact coupling
    │   ├── porosity_migration_model.py  # pore migration, CG (SU / SUPG) or DG + SSP-RK3
    │   ├── cluster_dynamic_model.py # 1D advection-diffusion (DG + SIPG + upwind)
    │   ├── nn_conductivity.py       # neural-network k(T), external operator
    │   ├── gpr_conductivity.py      # Gaussian-process correction of the Magni k
    │   └── magni_conductivity.py    # Magni MA-MOX k
    ├── coupling/                    # external-code coupling (e.g. SCIANTIX fission gas)
    ├── materials/                   # YAML cards + Python callables
    │   ├── steel.yaml, austenitic_steel.yaml, ..., vessel_steel.yaml
    │   ├── uo2.yaml, zircaloy.yaml
    │   ├── ceramic.yaml, oxide.yaml, plastic.yaml, lead.yaml, h2o.yaml
    │   └── ceramic.py, oxide.py, fuel_*.py, zircaloy_E.py  # k(T), Gc(x), swelling, E(T) callables
    ├── utils/                       # post-processing + helpers
    │   ├── writer.py                # unified VTU / XDMF OutputWriter
    │   ├── logger.py                # framework-wide logger
    │   ├── plot_convergence.py
    │   ├── utils_extract_vtu.py     # field extraction from VTU
    │   ├── utils_extract_xdmf.py    # same for XDMF
    │   ├── utils_load.py            # YAML loader + power-history generator
    │   ├── utils_plot.py            # 1D / radial plots
    │   ├── utils_verification.py    # analytical benchmarks
    │   └── geo_files/               # reusable Gmsh templates
    ├── ai/                          # agent onboarding (PROMPT.md, CONTEXT.md)
    ├── conference/                  # FEniCS 2026 materials (slides, demo, handout)
    ├── examples/                    # minimal didactic setups
    └── cases/                       # 70 cases carrying a gold, 65 in the local suite
        ├── verification/            # single-effect checks against a closed-form solution
        │   ├── thermal/             #   slabs, shells, heated box
        │   ├── mechanics/           #   Lamé, GPS, Mariotte, cylinders, cavities
        │   ├── plasticity/          #   J2 hardening, crystal-plasticity demo
        │   └── fuel/                #   swelling, burnup, creep, conductivity, law discovery
        ├── benchmarks/              # literature reproducers (Ambati SENT/SENS, McClenny
        │   └── damage/              #   pellet quench, Kamagate plate), all damage cases
        ├── regression/              # multi-physics configurations without a closed-form
        │                            # solution; some carry partial analytical references,
        │                            # all carry a gold
        ├── studies/                 # parametric sweeps (mesh sensitivity, attenuation map)
        ├── sandbox/                 # work in progress (never in the suite)
        ├── teaching/
        ├── non-regression_local.sh  # discovery-based local suite
        ├── non-regression_github.sh # CI suite (reads cases_ci.txt)
        ├── cases_ci.txt / suite_exclude.txt
        └── non-regression_summary.txt
```

The full case catalogue and per-module details are in `z3st/ai/CONTEXT.md`.

---

## Example input file

```yaml
# input.yaml
mesh_path: mesh.msh
geometry_path: geometry.yaml
boundary_conditions_path: boundary_conditions.yaml

materials:
  steel: ../../materials/steel.yaml

regime: 2d                       # 1d | 2d | 3d | axisymmetric

solver_settings:
  max_iters: 100
  relax_T: 0.9
  relax_u: 0.7
  relax_adaptive: true
  relax_growth: 1.2
  relax_shrink: 0.8
  relax_min: 0.05
  relax_max: 0.95

models:
  thermal: true
  mechanical: true
  # damage: true                  # enable phase-field fracture
  # plasticity: true              # enable J2 (or custom via plasticity.mode)
  # cluster: true                 # enable 1D cluster dynamics
  # gap_conductance:              # Robin pair-coupled gap on interfaces
  #   type: Fixed                 # or Gas
  #   value: 5000.0

mechanical:
  solver: linear                  # linear | nonlinear (SNES Newton, for hyperelastic / custom)
  linear_solver: iterative_hypre  # direct_mumps | iterative_amg | iterative_hypre
  rtol: 1.0e-5
  stag_tol: 1.0e-5
  convergence: rel_norm           # rel_norm | norm

thermal:
  analysis: stationary            # stationary | transient (backward Euler)
  solver: linear
  linear_solver: iterative_hypre
  rtol: 1.0e-6
  stag_tol: 1.0e-6
  convergence: rel_norm

# damage:                          # uncomment when damage is on
#   type: AT2                      # AT1 | AT2
#   linear_solver: iterative_hypre
#   rtol: 1.0e-6
#   stag_tol: 1.0e-3
#   convergence: rel_norm
#   lc: 1.0e-4
#   hybrid_constraint: true
#   split: amor                    # amor | miehe | star_convex (optional override; default = Amor for AT1, Miehe for AT2)
#   gamma_star: 0.0                # star-convex only: model parameter ≥ -1 controlling σ_c⁻/σ_c⁺ ratio; 0 ≡ Amor

lhr:
  - 0
time:
  - 0
n_steps: 1

# time_adaptivity:                 # optional, off by default; bisect dt on a stalled step
#   enabled: true
#   dt_min: 1.0e3                  # (s) smallest dt to attempt before aborting
#   max_cuts: 6                    # max bisection depth per original-grid step

output:
  format: vtu                     # vtu | xdmf
```

---

## Post-processing and visualisation

The output files open in ParaView and PyVista. The Python utilities below extract fields from them for scripted analysis.

| Tool                   | Description                                                                             |
| ---------------------- | --------------------------------------------------------------------------------------- |
| `writer.py`            | Unified `OutputWriter`: per-step VTU files or single-file XDMF time series              |
| `utils_extract_vtu.py` | Extracts scalar/vector fields and stress components from VTU outputs                    |
| `utils_plot.py`        | Generates 1D and radial plots (e.g. T(r), σ<sub>rr</sub>(r))                            |
| ...                    | ...                                                                                     |


---

## Python material functions

A material card can give a property as the dotted path of a Python function in place of a number (see [`z3st/materials/README.md`](z3st/materials/README.md) and the Usage page of the documentation). The function receives the current fields, for example the temperature, and returns a UFL expression that enters the weak form. Examples in the repository:

- temperature-dependent conductivity (`materials/fuel_thermal.py`, `materials/magni_mox_thermal.py`)
- a spatially graded fracture energy (`materials/oxide.py`)
- burnup-driven swelling eigenstrains (`materials/fuel_swelling.py`, `materials/sciantix_swelling.py`)
- a crystal-plasticity stress function in a case directory (`verification/plasticity/crystal_single_grain`)
- trained conductivity models (`type: neural_network` or `gpr` in the `k` block)

## Development roadmap

* Frictional and mortar contact (the current contact is penalty with a uniform pressure)
* Monolithic phase-field Newton solver for crack nucleation
* Multi-point constraints with `dolfinx_mpc`
* Material libraries through `dolfinx_materials`
* Coupling with microstructure generators (Mérope)
* Cluster dynamics with nucleation
* Coupling with rate-theory codes and Monte Carlo workflows

## Reproducing the results

The main results come from cases in this repository. Run `./Allrun` in the case directory. It meshes the case, runs Z3ST and then `non-regression.py`.

| Result | Case directory |
|---|---|
| Integral fuel rod: gap closure, PCMI, temperature | `z3st/cases/regression/pwr_rod_2D` |
| The same rod with the SCIANTIX coupling active | `z3st/cases/regression/fg_test_2D` |
| Phase-field damage after a cold-bath quench | `z3st/cases/benchmarks/damage/pellet_quench_2D_xy` |
| Radial porosity profile, CG and DG discretisations | `z3st/cases/verification/fuel/porosity_migration` (CG), `z3st/cases/verification/fuel/porosity_migration_dg` (DG) |
| Penalty contact vs. the analytical Lame interference fit | `z3st/cases/verification/fuel/shrink_fit` |
| Three-dimensional contact vs. the Lamé interference pressure | `z3st/cases/verification/fuel/shrink_fit_disk_3d` |
| Contact-pressure relaxation by creep vs. an independent radial J2 solution | `z3st/cases/verification/fuel/creep_shrink_fit_2D` |
| Three-dimensional gamma-heated spherical shell, mesh convergence | `z3st/cases/verification/thermal/spherical_shell` (`mesh_convergence.py`) |
| Two-dimensional mesh convergence of temperature and stress | `z3st/cases/studies/mesh_sensitivity_2D` |
| Gaussian-process posterior propagated to the melting margin | `z3st/cases/studies/gpr_uq_margin` (`run.py`) |

---

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Before committing, run the fifteen static checks. They take a few seconds and run no simulation:

```bash
python -m z3st.utils.audit_checks          # all; --list explains each one
python -m z3st.utils.audit_checks deps     # or one by name
```

They check the tree itself: a model no gold-carrying case reaches, a disabled or tautological
verdict, a case that writes a verdict with no gold, a NaN in a gold, an undefined name, a
dependency imported by library code but not declared, a method name colliding across
`Spine`'s 14 parent classes, a stale case path in the docs, a broken path in a driver, a shell
syntax error in any shell file, and a script the CI workflow invokes that does not exist.

To make git refuse a commit the checks reject (once per clone):

```bash
git config core.hooksPath .githooks
```

A pass means that no static defect was found. The checks run no case, compare no gold and
check no convergence. `z3st/cases/non-regression_local.sh` runs the cases.

## Contributor acknowledgements

Beyond the author, the following people have contributed to Z3ST:

* **Romain Turgis** (ENSTA Paris): the `verification/fuel/creep_shrink_fit_2D`
  case and its analysis scripts, on the relaxation of the pellet–cladding
  contact pressure by Norton creep, with the closed form of Esposito et al.,
  *Int. J. Pressure Vessels and Piping* **185** (2020) 104126.
  Internship at Politecnico di Milano, 2026-05-25 to 2026-07-31.

---

## FEniCSx project acknowledgement

Z3ST is built on the FEniCSx ecosystem. Recommended citations include:

### DOLFINX

```bibtex
@article{dolfinx2023,
  author  = {Baratta, I. A. and Dean, J. P. and Dokken, J. S. and Habera, M. and Hale, J. S. and Richardson, C. N. and Rognes, M. E. and Scroggs, M. W. and Sime, N. and Wells, G. N.},
  title   = {DOLFINx: The next generation FEniCS problem solving environment},
  year    = {2023},
  doi     = {10.5281/zenodo.10447666}
}
```

### BASIX

```bibtex
@article{basix2022a,
  author  = {Scroggs, M. W. and Dokken, J. S. and Richardson, C. N. and Wells, G. N.},
  title   = {Construction of arbitrary order finite element degree-of-freedom maps on polygonal and polyhedral cell meshes},
  journal = {ACM Transactions on Mathematical Software},
  volume  = {48},
  number  = {2},
  pages   = {18:1--18:23},
  year    = {2022},
  doi     = {10.1145/3524456}
}

@article{basix2022b,
  author  = {Scroggs, M. W. and Baratta, I. A. and Richardson, C. N. and Wells, G. N.},
  title   = {Basix: a runtime finite element basis evaluation library},
  journal = {Journal of Open Source Software},
  volume  = {7},
  number  = {73},
  pages   = {3982},
  year    = {2022},
  doi     = {10.21105/joss.03982}
}
```

### UFL

```bibtex
@article{ufl2014,
  author  = {Aln{\ae}s, M. S. and Logg, A. and {\O}lgaard, K. B. and Rognes, M. E. and Wells, G. N.},
  title   = {Unified Form Language: A domain-specific language for weak formulations of partial differential equations},
  journal = {ACM Transactions on Mathematical Software},
  volume  = {40},
  year    = {2014},
  doi     = {10.1145/2566630}
}
```

### External operators (optional, neural-network material laws)

The Newton route for data-driven conductivity builds on **dolfinx-external-operator**
(LGPL-3.0, https://github.com/a-latyshev/dolfinx-external-operator):

```bibtex
@article{latyshev2025externaloperators,
  author  = {Latyshev, Andrey and Bleyer, J{\'e}r{\'e}my and Maurini, Corrado and Hale, Jack S.},
  title   = {Expressing general constitutive models in {FEniCSx} using external operators and algorithmic automatic differentiation},
  journal = {Journal of Theoretical, Computational and Applied Mechanics},
  year    = {2025},
  doi     = {10.46298/jtcam.14449}
}
```


* **DOLFINx tutorial (J. S. Dokken):** https://jsdokken.com/dolfinx-tutorial/
* **FEniCS project documentation:** https://fenicsproject.org/documentation/

---

## License & author

If you use Z3ST in your research, please cite the archived software:

```bibtex
@software{Z3ST,
  author    = {Giovanni Zullo},
  title     = {Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.17748028},
  url       = {https://github.com/giozu/z3st}
}
```

This is the concept DOI. It covers every release and resolves to the most recent
one. To cite the exact code behind a result, use the version DOI from the Zenodo
page of the release you used, and give the version number.

A software paper describing the framework is in preparation.

* **Author:** Giovanni Zullo
* **Institution:** Politecnico di Milano
* **Version:** 0.4.0 (2026)
* **License:** Apache 2.0
