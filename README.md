# Finite element method

The finite element method (FEM) is a numerical method for solving engineering, physics and mathematical problems.

## FEniCSx

FEniCSx is a popular open-source computing platform, using FEM to solve partial differential equations (PDEs).
Comprehensive tutorials and useful links are:
* https://jsdokken.com/dolfinx-tutorial/
* https://github.com/Steriva/FEniCSx-tutorials
* https://fem-on-colab.github.io/packages.html
* https://github.com/RBniCS/RBniCS/tree/master/tutorials
* https://fenicsproject.org/pub/tutorial/sphinx1/
* https://docs.fenicsproject.org/dolfinx/main/python/demos/demo_gmsh.html

### How to install FEniCSx on Windows

#### WSL
On Windows, using WSL, write in Linux shell the following instruction to install gmshio:

```shell
sudo apt install python3-pip python3-gmsh
sudo apt install libgl1 xvfb
```

Then, to install FEniCSx:

```shell
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update && sudo apt full-upgrade -y
sudo apt install fenicsx -y
```

Additional modules are:

```shell
pip install pyvista numpy==1.26.4 trame
python3 -m pip install pyvista
```

If *ImportError: cannot import name '_docstring'* happens, or multiple matplotlib versions:

```shell
sudo apt remove python3-matplotlib
pip uninstall matplotlib
pip install matplotlib
```

### Documentation and references
A collection of tutorials is available at the official website https://docs.fenicsproject.org/dolfinx/main/python/demos.

# FEniCS
FEniCS is the old version of FEniCSx. FEniCSx is recommended. The latest version of legacy FEniCS (2019.1.0) was released in April 2019.