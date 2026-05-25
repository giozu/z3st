# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
[![DOI](https://zenodo.org/badge/648784453.svg)](https://doi.org/10.5281/zenodo.17748028)
[![CI](https://github.com/giozu/z3st/actions/workflows/ci.yml/badge.svg)](https://github.com/giozu/z3st/actions/workflows/ci.yml)
![static](https://github.com/giozu/z3st/actions/workflows/static.yml/badge.svg)
<!-- ![paper build](https://github.com/giozu/z3st/actions/workflows/paper.yml/badge.svg) -->

**Z3ST** is an open-source finite-element framework for thermo-mechanical material analysis.
Developed in **Python**, it leverages **FEniCSx** and provides a modular environment to couple heat conduction and linear elasticity in multi-material domains under stationary or transient conditions, with user-defined boundary conditions.

---

## Overview

Z3ST supports coupled staggered thermo-mechanical analysis, handling arbitrary geometries, materials, and time-dependent power histories.
It is designed to be **modular and extensible**, enabling rapid model development, integration with advanced solvers, and **easy coupling with external codes**.

---

## Quick installation

Clone the repository:
  ```bash
  git clone https://github.com/giozu/z3st.git
  ```

Z3ST requires a working installation of FEniCSx (dolfinx), that is not installed via pip and must be obtained via Conda (or Docker).
To install miniconda:
  ```bash
  mkdir -p ~/miniconda3
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
  bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
  rm ~/miniconda3/miniconda.sh
  source ~/miniconda3/bin/activate
  conda init --all
  ```

Create and activate the Conda environment (optional but recommended)
  ```bash
  cd z3st
  conda env create -f z3st_env.yml
  conda activate z3st
  ```

Install in editable/development mode
  ```bash
  pip install -e .
  ```

Then, it is possible to execute in each case folder, for instance:
  ```bash
  cd ~/z3st/z3st/cases/00_example/
  gmsh mesh.geo -3 > log_mesh.md
  python3 -m z3st > log_z3st.md
  python3 non-regression.py
  python3 ../../utils/plot_convergence.py
  ```

Or, with less instructions:
  ```bash
  cd ~/z3st/z3st/cases/00_example/
  ./Allrun
  ```

Optional flags:

* `python3 -m z3st --mesh_plot` → preview surface tags using PyVista

---

# Simulation Cases

This directory contains example simulation setups (*cases*) used for verification, validation,
and demonstration of the Z3ST thermo-mechanical solver.

Each subfolder is a simulation case, including:
- YAML input configuration (`input.yaml`)
- Geometry and meshing definitions (`geometry.yaml`, `mesh.geo`)
- Boundary conditions (`boundary_conditions.yaml`)
- A non-regression test script (`non-regression.py`)
- Reference results in `output/non-regression_gold.json`

These cases serve both as:
- regression tests,
- and examples for new users.

---

## Key features

* **Coupled thermo-mechanical solver** — heat conduction (stationary or transient, backward Euler) and mechanics with staggered coupling and adaptive relaxation
* **Multi-regime kinematics** — `2d` plane strain, `3d`, `axisymmetric`, and `plane_stress` available through a single configuration entry (the axisymmetric weight `w = 2πr` and cylindrical strain components are handled internally)
* **Constitutive laws** — small-strain isotropic Lamé, anisotropic Voigt (user-supplied 6×6 stiffness), Neo-Hookean hyperelasticity (SNES Newton with line search), J2 plasticity with linear isotropic hardening, and a `custom` hook for user-supplied UFL stress functions (used by the crystal-plasticity demo)
* **Phase-field fracture** — variational AT1 and AT2 models with Miehe spectral or Amor volumetric/deviatoric energy splits, irreversibility enforcement, and the Ambati-Gerasimov-De Lorenzis hybrid constraint
* **Multi-material domains** — independent thermal, mechanical, and damage properties per material; per-cell-tag integration measures handle interfaces naturally
* **Volumetric heating** — fissile (LHR/area), analytic γ-heating decay in rectangular / cylindrical / spherical geometry, or arbitrary user-defined `q'''(x)`
* **Flexible boundary conditions** — thermal (Dirichlet / Neumann / Robin with convective or gap-coupled mode), mechanical (Dirichlet vector / per-component / Neumann / Clamp / Slip), damage (Dirichlet); step-dependent value histories on mechanical BCs
* **Gap conductance model** — `Fixed` or `Gas` (`k_gas = c·10⁻⁴·T_gap^0.79`); paired facet groups, mean centroid gap-size computed once via SciPy cKDTree
* **1D cluster dynamics** — defect-cluster size-distribution solver: implicit-Euler DG1 with upwind interior-facet flux and SIPG diffusion, with mass-conservation rescaling
* **Python material modules** — temperature-dependent `k(T)`, spatially-graded `Gc(x)`, user-supplied stress functions, all loaded via `importlib` at runtime
* **Material database** — YAML-based cards (`materials/`): UO₂, multiple steel families (austenitic, martensitic, high-carbon, T91, 15-15Ti, vessel), Zircaloy, ceramics, oxides, plastic, lead, H₂O
* **Mesh input** — Gmsh `.msh` files or YAML-driven mesh builder; reusable Gmsh templates under `utils/geo_files/`
* **YAML-driven configuration** — three plain-text files per case (`input.yaml`, `geometry.yaml`, `boundary_conditions.yaml`); reproducible, diffable, version-controllable
* **Parallel performance** — PETSc with MUMPS / GAMG / HYPRE BoomerAMG, MPI via `MPI.COMM_WORLD`
* **Post-processing ecosystem** — VTU and XDMF time-series output through a unified writer that pre-compiles all interpolation expressions once at setup; ParaView- and PyVista-compatible
* **Continuous integration** — per-case `non-regression.py` vs. version-controlled gold JSON, summarised on every commit via GitHub Actions
* **Documented API** — Sphinx sources under `docs/source/`, built by GitHub Actions

---

## Directory structure

```bash
z3st/                                # repository root
├── LICENSE                          # Apache 2.0
├── README.md
├── CITATION.cff
├── CONTEXT.md                       # design + work-log document
├── pyproject.toml / setup.py
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
    │   ├── solver.py                # staggered solver, PETSc options
    │   ├── finite_element_setup.py  # V_t / V_m / V_d / V_c / V_pl / Q
    │   ├── diagnostic.py
    │   └── mesh/                    # Gmsh loader, MeshManager, PyVista preview
    │       ├── reader.py
    │       ├── manager.py
    │       └── plotter.py
    ├── models/                      # physics mixins plugged into Spine
    │   ├── thermal_model.py
    │   ├── mechanical_model.py      # lame / voigt / hyperelastic / plasticity / custom
    │   ├── damage_model.py          # AT1 / AT2, Miehe / Amor splits, hybrid constraint
    │   ├── plasticity_model.py      # J2 + custom CP hook
    │   ├── gap_model.py             # Fixed / Gas gap conductance
    │   └── cluster_dynamic_model.py # 1D advection–diffusion (DG + SIPG + upwind)
    ├── materials/                   # YAML cards + Python callables
    │   ├── steel.yaml, austenitic_steel.yaml, ..., vessel_steel.yaml
    │   ├── uo2.yaml, zircaloy.yaml
    │   ├── ceramic.yaml, oxide.yaml, plastic.yaml, lead.yaml, h2o.yaml
    │   └── ceramic.py, oxide.py     # k(T), Gc(mesh) callables
    ├── utils/                       # post-processing + helpers
    │   ├── writer.py                # unified VTU / XDMF OutputWriter
    │   ├── export_vtu.py            # legacy VTU writer (backward compat)
    │   ├── mesh_builder.py
    │   ├── plot_convergence.py
    │   ├── utils_extract_vtu.py     # field extraction from VTU
    │   ├── utils_extract_xdmf.py    # same for XDMF
    │   ├── utils_load.py            # YAML loader + power-history generator
    │   ├── utils_plot.py            # 1D / radial plots
    │   ├── utils_verification.py    # analytical benchmarks
    │   ├── output.py                # stdout / JSON helpers
    │   ├── z-gui.py                 # interactive PyVista viewer
    │   └── geo_files/               # reusable Gmsh templates
    ├── examples/                    # minimal didactic setups
    └── cases/                       # ~50 verification / validation / demo cases
        ├── 00_example/              # tutorial: uniaxial steel block 3D
        ├── 1_thin_slab_2D/          # first thermal slab
        ├── ...
        ├── 14_full_cylinder_cracking_2D_xy/   # UO2 thermal-shock + AT1 (paper flagship)
        ├── 19_single-edge_notched_*/          # SENT / SENS phase-field benchmarks
        ├── 20_plasticity_2D/
        ├── demo_CP_single_grain/             # custom crystal-plasticity demo
        ├── non-regression.sh / .py / _github.sh   # regression infrastructure
        └── non-regression_summary.txt
```

The full case catalogue and per-module details are in `CONTEXT.md`.

---

## Example input file

```yaml
# input.yaml
mesh_path: mesh.msh
geometry_path: geometry.yaml
boundary_conditions_path: boundary_conditions.yaml

materials:
  steel: ../../materials/steel.yaml

regime: 2d                       # 2d | 3d | axisymmetric | plane_stress

solver_settings:
  coupling: staggered
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
#   solver: linear
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

output:
  format: vtu                     # vtu | xdmf
```

---

## Post-processing and visualization

Z3ST provides a flexible, Python-based post-processing environment that can be used directly or extended by the user.
Displacement, temperature, and stress fields can be extracted, analyzed, and visualized through built-in utilities or custom scripts.
Full compatibility with **ParaView** and **PyVista** enables both automated and interactive visualization workflows.

| Tool                   | Description                                                                             |
| ---------------------- | --------------------------------------------------------------------------------------- |
| `export_vtu.py`        | Exports temperature, displacement, strain, and stress fields to `.vtu` format           |
| `utils_extract_vtu.py` | Extracts scalar/vector fields and stress components from VTU outputs                    |
| `utils_plot.py`        | Generates 1D and radial plots (e.g. T(r), σ<sub>rr</sub>(r)) and can be easily extended |
| `z-gui.py`             | Interactive 3D viewer built on PyVista for exploratory visualization                    |
| ...                    | ...                                                                                     |


These tools, together with other Python utilities, provide a complete post-processing ecosystem for field analysis and visualization.

---

## Capabilities beyond YAML

While Z3ST supports full configuration through YAML files for reproducibility and ease of use, the framework is deliberately designed to extend far beyond YAML-based inputs.

### External Python material models
Material behaviour does not need to be hard-coded or restricted to values inside YAML files.
Z3ST natively supports **Python-based material modules**, allowing users to implement:

- temperature-dependent constitutive laws
- empirical or semi-empirical correlations
- nonlinear stress–strain behaviour
- surrogate models and neural networks
- models depending on local fields (e.g., temperature, gradients)

These Python modules are automatically loaded and integrated into the finite-element formulation at runtime.

### Integration with FEniCSx ecosystem
Z3ST is designed to interoperate with the wider FEniCSx ecosystem, with planned integration of:

- **`dolfinx_mpc`** for multi-point constraints and advanced boundary-condition enforcement
- **`dolfinx_materials`** for standardised material libraries and automatically tabulated thermal/mechanical properties

This allows the framework to grow toward more general multiphysics and nonlinear-modelling workflows.

### Microstructure and multi-scale workflows
Z3ST is also envisaged to integrate with external tools for:

- synthetic microstructure generation (e.g., **Merope**)
- **rate-theory solvers** for irradiation-driven defect evolution
- **Monte Carlo workflows** for stochastic microstructural processes

These extensions aim to connect Z3ST to multi-scale modelling pipelines involving microstructure → mesoscale → continuum thermo-mechanics.

### Development roadmap

* Contact mechanics
* Nonlinear constitutive behavior
* Coupling with microstructure generators
* Advanced cluster dynamics (1D, nucleation)
* Coupling with rate-theory codes

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

---

## FEniCSx project acknowledgement

Z3ST is built on the **FEniCSx** ecosystem. Recommended citations include:

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


* **DOLFINx tutorial (J. S. Dokken):** https://jsdokken.com/dolfinx-tutorial/
* **FEniCS project documentation:** https://fenicsproject.org/documentation/

---

## License & author

If you use Z3ST in your research, please cite this little work.

```bibtex
@misc{Z3ST2025,
  author       = {Giovanni Zullo},
  title        = {Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis},
  year         = {2025},
  howpublished = {\url{https://github.com/giozu/z3st}},
  note         = {Version 0.1.0}
}
```

* **Author:** Giovanni Zullo
* **Institution:** Politecnico di Milano
* **Version:** 0.1.0 (2025)
* **License:** Apache 2.0
