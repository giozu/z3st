# ZEST: A FEniCSx Framework for Fuel Performance Modeling

## Overview

**ZEST** is an open-source thermomechanical fuel performance code written in Python using the **FEniCSx** framework. It allows the simulation of coupled thermal and mechanical behavior of nuclear fuel elements, supporting both steady-state and transient scenarios.

Developed by Giovanni Zullo, ZEST is modular, extensible, and supports coupling with external multiphysics codes (e.g. OpenMC, OFFBEAT).

---

## Features

* **Thermal and Mechanical Coupling**: Fully coupled thermoelastic solver for fuel and cladding.
* **FEM-based with FEniCSx**: Leverages the FEniCSx ecosystem for PDE solution and function space management.
* **Subdomain-aware**: Handles fuel and cladding regions with different material properties.
* **Time-dependent power history**: Supports LHR input interpolation from tabulated data.
* **Material YAML support**: Materials are read from easy-to-edit `.yaml` files.
* **Mesh I/O**: Structured mesh input from `.msh` files via Gmsh.
* **Logging and output**: Saves log files and results for post-processing.
* **Post-processing**: Integrated with PyVista, Matplotlib, and XDMF output for ParaView.

---

## Structure

```
ZEST/
 ├── zest.py                    # Main script
 ├── fuel_rod.py               # Core thermo-mechanical class
 ├── geometryManager.py        # Manages geometry update/deformation
 ├── FiniteElementSetups.py  # FEM setup and function spaces
 ├── utils.py                  # Helpers (power history, KSP log)
 ├── materials/
 │    ├── uo2.yaml
 │    ├── T91.yaml
 │    ├── lead.yaml
 │    └── he.yaml
 └── mesh/                    # Mesh files (.msh)
```

---

## Example Usage

```python
from fuel_rod import FuelRod
from geometry.fuel_rod import geometry
from utils import generate_power_history
from materials.loader import load_material

fuel     = load_material("materials/uo2.yaml")
cladding = load_material("materials/T91.yaml")
gas      = load_material("materials/he.yaml")

problem = FuelRod(mesh_file='mesh/fuel_cladding_3d_structured.msh')
problem.volume(geometry["fuel outer radius"], geometry["cladding inner radius"],
               geometry["cladding outer radius"], geometry["height"])

problem.parameters(fuel, cladding, gas,
                   g=9.81,
                   lhr=25000,
                   T_ref=873.0)

problem.boundary_conditions(gap_pressure=1.0e6,
                            coolant_pressure=1.5e6,
                            temperature=873.0)

problem.initialize_fields()
problem.initialize_geometry()
problem.solve()
```

---

## YAML Material Example

```yaml
name: UO2
rho: 10970           # kg/m^3
E: 1.7e11            # Pa
nu: 0.3              # -
k: 2.5               # W/m-K
cp: 300              # J/kg-K
alpha: 1.45e-5       # 1/K
```

---

## Environment Setup (WSL + Conda)

### 1. Create the Conda Environment

```bash
conda create -n zest python=3.12.9 -y
conda activate zest
```

### 2. Install Required Packages

```bash
conda install -c conda-forge openmc fenics-dolfinx jupyterlab matplotlib ipywidgets pyvista meshio
pip install gmsh pyyaml
```

### 3. Test Installation

```python
import openmc, dolfinx, pyvista, meshio, gmsh, yaml
print("Everything is working!")
```

### 4. Export/Import Environment

```bash
conda env export > zest_environment.yml
conda env create -f zest_environment.yml
```

---

## Mesh Generation

To create a Gmsh mesh:

```bash
gmsh your_geometry.geo -3 -format msh2 -o your_mesh.msh
```

---

## References and Tutorials

* [FEniCSx official demos](https://docs.fenicsproject.org/dolfinx/main/python/demos/)
* [J.S. Dokken tutorials](https://jsdokken.com/dolfinx-tutorial/)
* [OpenMC documentation](https://docs.openmc.org/)

---

## Future Work

* Time integration and transient heat transfer
* Advanced creep and irradiation models
* Multiphysics coupling (neutronics, chemistry)
* Nonlinear mechanical solvers and contact mechanics

---

## License and Contributions

ZEST is open-source. Contributions, feedback, and feature requests are welcome!

---

**Author**: Giovanni Zullo
**Institution**: Politecnico di Milano
**Project**: OperaHPC / SCIANTIX / ZEST
**Date**: 2025
