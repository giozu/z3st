# Z3ST: A FEniCSx Framework for Thermo-Mechanical Analysis

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

Install miniconda:
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
  cd ~/z3st/cases/19_plate/
  gmsh -3 mesh.geo
  python3 -m z3st > log.z3st
  python3 non-regression.py
  python3 ../../utils/plot_convergence.py
  ```

Or, with a single instruction:
  ```bash
  cd ~/z3st/cases/19_plate/
  ./Allrun
  ```

Optional flags:

* `python3 -m z3st --mesh_plot` → preview surface tags using PyVista

---

## Key features

* **Coupled thermo-mechanical solver** — heat conduction and linear elasticity with monolithic or staggered coupling
* **Multi-material domains** — multiple materials with independent thermal and mechanical properties
* **Volumetric heating** — arbitrary, spatially dependent internal heat sources (e.g. fission, γ-heating, user-defined functions)
* **Flexible boundary conditions** — Dirichlet, Neumann, clamp, and slip BCs defined via YAML configuration files
* **Material database** — YAML-based definitions (`materials/`) including steels, zircaloys, oxides, and other structural or fuel materials
* **Geometry & mesh input** — accepts Gmsh `.msh` files or Python-generated meshes from YAML geometry descriptors; compatible with any meshing tool
* **Built-in gap conductance model** — thermal coupling between subdomains via Robin BCs or fixed conductance
* **YAML-driven configuration** — all parameters, materials, and boundary conditions defined externally for reproducibility
* **Post-processing ecosystem** — Python-based tools for field extraction and visualization, fully compatible with ParaView and PyVista
* **Fully documented API** — comprehensive Sphinx documentation under `docs/`

---

## Directory structure

```bash

z3st/
├── LICENSE
├── README.md
├── ...
│
├── z3st/                        # Python package (entry point, CLI)
│   ├── __main__.py
│   ├── __init__.py
│
├── core/                        # FEM core
│   ├── config.py
│   ├── mesh.py
│   ├── solver.py
│   ├── spine.py
│   ├── finite_element_setup.py
│   ├── diagnostic.py
│   └── __init__.py
│
├── models/                      # Physical models
│   ├── thermal_model.py
│   ├── mechanical_model.py
│   ├── gap_model.py
│   └── __init__.py
│
├── materials/                   # Material databases (YAML + helpers)
│   ├── steel.yaml
│   ├── ...
│   └── __init__.py
│
├── utils/                       # Post-processing and helpers
│   ├── export_vtu.py
│   ├── plot_convergence.py
│   ├── ...
│   └── __init__.py
│
├── cases/                       # Verification and demonstration cases
│   ├── 1_thin_thermal_slab/
│   ├── ...
│   ├── non-regression.sh
│   └── non-regression_summary.txt
│
└── docs/                        # Sphinx documentation (built by GitHub Actions)
    ├── Makefile
    ├── source/
    └── build/                   # auto-generated, not committed

```

---

## Example input file

```yaml
# input.yaml
mesh_path: mesh.msh
geometry_path: geometry.yaml
boundary_conditions_path: boundary_conditions.yaml

materials:
  steel: ../../materials/steel.yaml

solver_settings:
  coupling: staggered
  max_iters: 100
  relax_T: 0.9
  relax_u: 0.7
  relax_adaptive: True
  relax_growth: 1.2
  relax_shrink: 0.8
  relax_min: 0.05
  relax_max: 0.95

models:
  thermal: True
  mechanical: True

mechanical:
  solver: linear
  linear_solver: iterative_amg # direct_mumps | iterative_amg | iterative_hypre
  rtol: 1.0e-5
  stag_tol: 1.0e-5
  mechanical_regime: "3D" # 3D | plane_stress
  convergence: rel_norm # rel_norm | norm
  debug: False

thermal:
  solver: linear
  linear_solver: iterative_amg # direct_mumps | iterative_amg | iterative_hypre
  rtol: 1.0e-6
  stag_tol: 1.0e-6
  convergence: rel_norm # rel_norm | norm

lhr:
  - 0
time:
  - 0
n_steps: 1
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

## Development roadmap

* Monolithic thermo-mechanical solver
* Contact mechanics
* Nonlinear constitutive behavior
* Coupling with microstructure generators (e.g., Merope)
* Coupling with rate-theory codes

---

## License & author

**Author:** Giovanni Zullo
**Institution:** Politecnico di Milano
**Version:** 0.1.0 (2025)
**License:** Apache 2.0