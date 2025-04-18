# Setting up a ZEST Environment with OpenMC and FEniCSx on WSL

This guide walks you through the installation of a Python environment named `zest` using Conda, specifically tailored for nuclear fuel simulations with **OpenMC** and **FEniCSx**, on **WSL (Windows Subsystem for Linux)**.

---

## ✅ Prerequisites

- You have **WSL** installed and configured (Ubuntu recommended).
- You have installed **Miniconda** or **Anaconda** in your WSL environment.
  > If not, follow these steps:
  ```bash
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ```

---

## 1. Create the Conda Environment

Open your WSL terminal and run:

```bash
conda create -n zest python=3.11 -y
conda activate zest
```

This creates and activates a new environment named `zest`.

---

## 2. Install OpenMC

Install **OpenMC** from the conda-forge channel:

```bash
conda install -c conda-forge openmc
```

This will also install all required dependencies such as `h5py`, `numpy`, etc.

---

## 3. Install FEniCSx (DOLFINx)

Install the **FEniCSx** core packages from conda-forge:

```bash
conda install -c conda-forge fenics-dolfinx
```

This includes:
- `dolfinx`
- `ufl`
- `ffcx`
- `petsc4py`
- and other necessary packages

---

## 4. Optional: Install Additional Tools

If you plan to use **JupyterLab**, **Matplotlib**, **PyVista**, **meshio**, and **Gmsh**:

```bash
conda install -c conda-forge jupyterlab matplotlib ipywidgets pyvista meshio
```

Then install `gmsh` via `pip` to ensure the Python module is available:

```bash
pip install gmsh
```

> **PyVista** is useful for advanced 3D visualization, particularly when post-processing results from FEniCSx or OpenMC.  
> **meshio** is useful for reading/writing mesh files in different formats such as `.xdmf`, `.vtk`, or `.msh`.  
> **Gmsh** is required by `dolfinx.io.gmshio.read_from_msh` to import `.msh` mesh files. The Python module must be installed via `pip`.

---

## 5. Verify the Installation

Run Python inside the `zest` environment:

```bash
python
```

Then test the libraries:

```python
import openmc
import dolfinx
import pyvista
import meshio
import gmsh
print("Everything is working!")
```

To exit Python:

```python
exit()
```

---

## 🧪 Notes

- You can export your environment for future sharing or reproducibility:

```bash
conda env export > zest_environment.yml
```

- To recreate it on another system:

```bash
conda env create -f zest_environment.yml
```

---

## 🛠️ Useful Links

- [OpenMC documentation](https://docs.openmc.org/)
- [FEniCSx documentation](https://docs.fenicsproject.org/)
- [conda-forge](https://conda-forge.org/)
- [PyVista documentation](https://docs.pyvista.org/)
- [meshio documentation](https://github.com/nschloe/meshio)
- [Gmsh documentation](https://gmsh.info/)

---

Happy Simulating! ⚛️🔥
