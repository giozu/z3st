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

#### Anaconda
After having installed [miniconda](https://conda.io), the latest stable release of FEniCSx can be installed through:

```shell
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista
```

Then, the new fenicsx-env environment must be selected when running a python script or notebook (e.g., in Visual Studio Code).
