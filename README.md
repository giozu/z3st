# Finite element method

The finite element method (FEM) is a numerical method for solving engineering, physics and mathematical problems.
Several ways to use this method are available. 

## FEniCSx

FEniCSx is a popular open-source computing platform, using FEM to solve partial differential equations (PDEs).

### How to install FEniCSx on Windows

#### WSL
On Windows, using WSL, write in Linux shell the following instructions:

```shell
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update
sudo apt upgrade
sudo apt install fenicsx
```

### Documentation and references
A collection of tutorials is available at the official website https://docs.fenicsproject.org/dolfinx/main/python/demos.

# FEniCS
FEniCS is the old version of FEniCSx. FEniCSx is recommended. The latest version of legacy FEniCS (2019.1.0) was released in April 2019.