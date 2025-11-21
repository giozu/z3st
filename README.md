# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis

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
  cd ~/z3st/cases/00_example/
  gmsh -3 mesh.geo
  python3 -m z3st > log.z3st
  python3 non-regression.py
  python3 ../../utils/plot_convergence.py
  ```

Or, with less instructions:
  ```bash
  cd ~/z3st/cases/00_example/
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

## Using z3stGPT (AI assistant)

Z3ST includes an optional AI assistant, **z3stGPT**, designed to facilitate code exploration, onboarding, and framework extension.

**Access the assistant:**  
**https://chatgpt.com/g/g-68372061d184819192c6c441e3ea4985-z3stgpt**

### What it can do

z3stGPT provides **read-only** access to the official repository and can:

- browse the complete directory structure  
- perform structured code reviews  
- explain solver logic, FEM formulations, and workflow architecture  
- analyze YAML configurations, boundary conditions, and material models  
- generate example scripts and debugging guidance  
- assist new users in understanding and extending the Z3ST framework  

### Safety and limitations

- z3stGPT **cannot modify files**, open pull requests, or push changes  
- all repository interactions are **read-only**  
- privacy policy: https://giozu.github.io/z3st/privacy.html  

This assistant supports users and contributors in learning and working with the Z3ST framework.

---

## Key features

* **Coupled thermo-mechanical solver** — heat conduction and linear elasticity with staggered coupling
* **Multi-material domains** — multiple materials with independent thermal and mechanical properties
* **Volumetric heating** — arbitrary, spatially dependent internal heat sources (e.g., γ-heating, user-defined functions)
* **Flexible boundary conditions** — Dirichlet, Neumann, clamp, and slip BCs defined via YAML configuration files
* **Material database** — YAML-based definitions (`materials/`) including steels, oxides, and other structural materials
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

* Full test and setup of the monolithic thermo-mechanical solver
* Contact mechanics
* Nonlinear constitutive behavior
* Coupling with microstructure generators (e.g., Merope)
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
  author  = {Scrogggs, M. W. and Dokken, J. S. and Richardson, C. N. and Wells, G. N.},
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