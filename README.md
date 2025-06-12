# ZEST: A FEniCSx Framework for Fuel Performance Modeling

## Overview

**ZEST** (Zero-gap Enhanced Simulation Tool) is an open-source thermomechanical fuel performance code written in Python using the **FEniCSx** finite element framework. It supports the simulation of coupled thermal and mechanical behavior of nuclear fuel elements, with special focus on multi-material configurations, flexible geometry definitions, and fully transient boundary conditions.

Developed by Giovanni Zullo, ZEST is modular, extensible, and ready to integrate advanced physics models for fuel performance simulations.

---

## Main Features

* **Coupled Thermo-Mechanical Solver**: Thermal conduction and linear elasticity solved in a coupled framework with fixed-point iterations.
* **Multi-material domain**: Support for arbitrary number of materials, each with independent mechanical and thermal properties.
* **Mesh input via Gmsh**: Import of structured/unstructured meshes in `.msh` format.
* **Geometry input via YAML**: Geometry and physical tags automatically parsed from mesh files.
* **Configurable via YAML input files**: All simulation settings, materials and boundary conditions are fully configurable.
* **Power history generator**: LHR history interpolated from tabulated data.
* **Time-dependent solver**: Full support for transient time stepping.
* **Relaxation schemes**: User-configurable relaxation factors for thermo-mechanical coupling.
* **Post-processing**: Integrated plotting utilities (temperature, displacement, stress, residuals).
* **Output in XDMF format**: Compatible with ParaView for 3D postprocessing.

---

## Directory Structure

```bash
ZEST/
 ├── zest.py                # Main driver script
 ├── core/
 │    ├── config.py         # Input parsing and configuration
 │    ├── mesh.py           # Mesh import and geometry setup
 │    ├── solver.py         # Global solver control
 │    ├── finiteElementSetup.py  # FEM spaces definition
 ├── models/
 │    ├── thermal_model.py  # Thermal model
 │    └── mechanical_model.py # Mechanical model
 ├── multi_material.py  # Main MultiMaterial class
 ├── msh_to_yaml.py     # Utility to generate YAML geometry from .msh
 ├── utils.py           # Power history, YAML I/O, log utilities
 ├── postprocessing.py  # Postprocessing and plotting
 ├── output.py          # Plotting 3D and 1D fields
 ├── materials/         # YAML material files
 ├── conditions/        # YAML boundary condition files
 ├── mesh/              # Gmsh meshes
 └── input.yaml         # Main input file for each simulation
```

---

## Input Example (input.yaml)

```yaml
models:
  thermal: true
  mechanical: true
  creep: false

time: [0.0, 5.0e6]
lhr: [25000.0, 25000.0]
n_steps: 10

geometry_path: geometry/multi_material.yaml
mesh_path: mesh/multi_material.msh
conditions_path: conditions/boundary_conditions.yaml

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

## Mesh to YAML helper

You can automatically extract physical tags from your Gmsh mesh:

```bash
python msh_to_yaml.py
```

---

## Example run

```bash
python zest.py --debug --mesh_plot
```

---

## Environment Setup

```bash
# Create environment
conda create -n zest python=3.12.9 -y
conda activate zest

# Install packages
conda install -c conda-forge fenics-dolfinx pyvista matplotlib meshio openmc
pip install gmsh pyyaml
```

---

## Postprocessing Tools

* `plot3d_temperature(problem)`: 3D temperature field (deformed mesh)
* `plot1d_temperature(problem, z_target=...)`: 1D radial temperature profiles
* `plot1d_stress(problem, z_target=...)`: 1D stress profiles
* `plot_residuals()`: Convergence plots for thermal and mechanical solvers
* `save_temperature_history()`: Exports temperature evolution in TSV format

---

## Next Development Steps

* Nonlinear creep, swelling, irradiation models
* Gap conductance updates (h\_gap(T,u))
* Nonlinear contact mechanics (fuel-clad gap)
* Material property dependencies with temperature and stress
* Support for more general geometry inputs
* Optimization of matrix assembly (reuse A)
* Integration with neutronics (OpenMC, Serpent)

---

## License and author

**License**: Apache License, Version 2.0 (see LICENSE file)

**Author**: Giovanni Zullo
**Institution**: Politecnico di Milano
**Year**: 2025
