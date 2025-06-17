# python modules
import dolfinx
import yaml
from mpi4py import MPI
import numpy as np
import os
from petsc4py import PETSc
import ufl
# ZEST modules
from core.config import Config
from core.mesh import Mesh
from core.solver import Solver
from core.finite_element_setup import FiniteElementSetup
from models.thermal_model import ThermalModel
from models.mechanical_model import MechanicalModel

class MultiMaterial(Config, Mesh, FiniteElementSetup, Solver, ThermalModel, MechanicalModel):
    def __init__(self, input_file, mesh_file, geometry, rank=0, gdim=3):
        Config.__init__(self, input_file)
        Mesh.__init__(self, mesh_file=mesh_file, geometry=geometry, rank=rank, gdim=gdim)
        FiniteElementSetup.__init__(self)
        Solver.__init__(self)
        ThermalModel.__init__(self)
        MechanicalModel.__init__(self)

    def parameters(self, lhr):
        self.g = 9.81 # m/s2
        self.lhr = lhr

    def materials(self, **materials):
        self.materials = {}

        for name, mat in materials.items():
            mat["lmbda"] = mat["E"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
            mat["mu"] = mat["E"] / (2 * (1 + mat["nu"]))
            mat["bulk_modulus"] = mat["E"] / (3 * (1 - 2 * mat["nu"]))
            self.materials[name] = mat

            print(f"\n{name.capitalize()} material loaded:")
            for k in sorted(mat.keys()):
                v = mat[k]
                print(f"  {k:<15} → {v} ({type(v).__name__})")

    def set_boundary_conditions(self):
        print(f"-----------------------")
        print(f"set_boundary_conditions")
        print(f"-----------------------")

        with open(self.config["boundary_conditions_path"], "r") as f:
            print(f"Loading boundary conditions from '{self.config["boundary_conditions_path"]}'")
            self.boundary_conditions = yaml.safe_load(f)

        if self.coupling == "staggered":
            self.set_thermal_boundary_conditions(self.V_t)
            self.set_mechanical_boundary_conditions(self.V_m)
        else:
            V_t_sub, V_t_map = self.W.sub(1).collapse()
            V_t_map = np.array(V_t_map, dtype=np.int32)
            V_u_sub, V_u_map = self.W.sub(0).collapse() # need to .collapse() to define BCs / define Functions from it / interpolate
            V_u_map = np.array(V_u_map, dtype=np.int32)
            self.set_thermal_boundary_conditions(V_t_sub, V_t_map)
            self.set_mechanical_boundary_conditions(V_u_sub, V_u_map)
    
    def initialize_fields(self):
        print(f"-----------------")
        print(f"initialize_fields")
        print(f"-----------------")

        # --- Temperature fields ---
        self.T = {}

        for name, mat in self.materials.items():
            T_func = dolfinx.fem.Function(self.V_t, name=f"Temperature_{name}")
            T_func.x.array[:] = mat["T_ref"]
            values = T_func.x.array
            print(f"  {name}: T min={values.min():.2f} K, T max={values.max():.2f} K, T mean={values.mean():.2f} K")
            self.T[name] = T_func

        # --- Heat source ---
        self.q_third = {}
        for name, mat in self.materials.items():
            is_fissile = bool(mat.get("fissile", False))
            self.q_third[name] = self.lhr / self.area if is_fissile else 0.0
            print(f"  q_third[{name}] = {self.q_third[name]:.2f} W/m³ (fissile: {is_fissile})")

        mu = 24 # (1/m) attenuation coefficient / macroscopic cross section of gamma rays in steels
        q_func = dolfinx.fem.Function(self.V_t)
        def f(x):
            q0 = 3.46e6 # W/m3
            return q0 * np.exp(-x[0] / mu)
        q_func.interpolate(f)
        self.q_third["steel"] = q_func

        # plot q_third(x)
        import pyvista
        V = self.V_t
        cells, types, x = dolfinx.plot.vtk_mesh(V)
        grid = pyvista.UnstructuredGrid(cells, types, x)
        grid.point_data["q"] = q_func.x.array.reshape((x.shape[0], 1))
        grid.set_active_scalars("q")
        subplotter = pyvista.Plotter(shape=(1, 1))
        subplotter.subplot(0, 0)
        subplotter.add_mesh(grid, show_edges=True, show_scalar_bar=True)
        subplotter.view_xy()
        subplotter.show()

        # --- Displacement field ---
        self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
        self.u.x.array[:] = 0.0
        values = self.u.x.array
        print(f"  u min={values.min():.2f} m, u max={values.max():.2f} m, u mean={values.mean():.2f} m")

    def solve(self, max_iters=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6, tol=1e-5):
        print("\nSolver settings:")
        print(f"  → Coupling              : {self.coupling}")

        if self.coupling == "monolithic":
            self.solve_monolithic(tol=tol, max_iter=max_iters)

        elif self.coupling == "staggered":
            self.solve_staggered(max_iter=max_iters, tol_u=tol_u, tol_T=tol_T, tol_mech=tol_mech, tol_th=tol_th)

        else:
            raise ValueError(f"Unknown coupling strategy: {self.coupling}")

    def get_results(self):
        self.strain = {}
        self.stress = {}
        
        u = self.u

        for name in self.materials:
            mat = self.materials[name]
            self.strain[name] = self.epsilon(u)
            self.stress[name] = self.sigma_mech(u, mat) - self.sigma_th(self.T[name], mat)

        return self.u, self.strain, self.stress, self.T


