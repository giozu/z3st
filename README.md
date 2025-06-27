# ZEST: An Open-Source Framework for Thermo-Mechanical Analysis

## Overview

**ZEST** is an open-source thermo-mechanical fuel performance code built on the **FEniCSx** finite element python library.
It supports the simulation of thermal and mechanical behavior of nuclear fuel elements, with flexibility in geometry, time-dependent loading, and boundary conditions.

Developed by Giovanni Zullo, ZEST is modular, extensible, and designed to incorporate advanced physics-based models.

---

## Key Features

* **Coupled Thermo-Mechanical Solver**: Thermal conduction and linear elasticity solved in a coupled framework (staggered) with fixed-point iterations
* **Transient support**: Fully time-dependent solver with LHR history interpolated from tabulated input
* **Multi-material domain**: Arbitrary number of materials, each with independent mechanical and thermal properties
* **Volumetric heating**: Includes fission and gamma heating sources with exponential attenuation
* **Flexible boundary conditions**: Dirichlet, Neumann, clamp, and slip conditions defined via YAML
* **Geometry input via YAML**: Automatic parsing of physical tags from `.msh` mesh files
* **Mesh input via Gmsh**: Support for structured and unstructured `.msh` files (2D/3D)
* **Configurable via YAML**: All simulation parameters, materials, and boundary conditions are user-configurable
* **Relaxation schemes**: Customizable relaxation factors for temperature and displacement
* **Post-processing**: Integrated 1D/3D plotting tools (temperature, stress, displacement, residuals)
* **Output in VTU/XDMF format**: Compatible with ParaView and PyVista for full-field visualization
* **Residual tracking**: Thermal and mechanical convergence plots for every time step

---

## Directory Structure

```bash
ZEST/
├── zest.py              # Main driver
├── z-gui.py             # VTU visualizer
├── input.yaml           # Simulation input
├── README.md / LICENSE
├── Allrun / Allclean
│
├── core/                # FEM setup & solver
│   ├── config.py
│   ├── mesh.py
│   ├── solver.py
│   └── ...
│
├── models/              # Physics models
│   ├── thermal_model.py
│   └── mechanical_model.py
│
├── output/              # Plots, exports, ParaView
│   ├── output.py
│   ├── *.pvsm
│   └── ...
│
├── geometry/            # YAML + Gmsh converter
│   ├── *.yaml
│   └── msh_to_yaml.py
│
├── mesh/                # Gmsh geometries
│   ├── *.geo / *.msh
│   └── generators (*.py)
│
├── materials/           # Material properties (YAML)
│   ├── *.yaml
│   └── *.ipynb
│
├── conditions/          # Boundary conditions (YAML)
├── verification/        # Benchmarks and tests
├── zest-api/            # Render web API
└── input/               # Example runs
```

---

## Sample Input File

```yaml
models:
  thermal: true
  mechanical: true
  creep: false

time: [0.0, 5.0e6]
lhr: [25000.0, 25000.0]
n_steps: 10

geometry_path: geometry/full_cylinder.yaml
mesh_path: mesh/full_cylinder.msh
boundary_conditions_path: conditions/boundary_conditions_full_cylinder.yaml

materials:
  fuel: materials/uo2.yaml
  cladding: materials/T91.yaml

solver_settings:
  tol_u: 1e-2
  tol_T: 1e-2
  tol_mech: 1e-6
  tol_th: 1e-6
  n_iter_max: 20
```

---

## How to Run

```bash
python zest.py --debug --mesh_plot
```

---

## Postprocessing Tools

* `z-gui.py`: interactive 3D viewer with PyVista
* `plot1d_temperature(...)`: 1D temperature profiles
* `plot1d_stress(...)`: Stress components in cylindrical or Cartesian form
* `plot_residuals()`: Solver residuals over iterations
* `export_vtu()` / `export_xdmf()`: Export fields for ParaView

---

## Environment Setup

```bash
conda create -n zest python=3.12.9 -y
conda activate zest

conda install -c conda-forge fenics-dolfinx pyvista matplotlib meshio openmc jupyterlab ipywidgets pandas numpy
pip install gmsh pyyaml
```

---

## Development Roadmap

* [ ] Monolithic coupled thermo-mechanical solver
* [ ] Nonlinear effects: creep, swelling, densification
* [ ] Temperature- and stress-dependent material models
* [ ] Contact mechanics for fuel-clad gap
* [ ] Advanced gap conductance models
* [ ] Coupling with OpenMC or Serpent
* [ ] Thermochemistry modules
* [ ] ROM-based model acceleration

---

## License and Author

**Author**: Giovanni Zullo
**Institution**: Politecnico di Milano
**Year**: 2025
**License**: Apache License 2.0
