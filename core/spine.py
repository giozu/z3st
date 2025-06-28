# --- Python modules ---
import dolfinx
import yaml
from mpi4py import MPI
import numpy as np
import os
from petsc4py import PETSc
import ufl
# --- ZEST modules ---
from core.config import Config
from core.mesh import Mesh
from core.solver import Solver
from core.finite_element_setup import FiniteElementSetup
from models.thermal_model import ThermalModel
from models.mechanical_model import MechanicalModel

class Spine(Config, Mesh, FiniteElementSetup, Solver, ThermalModel, MechanicalModel):
    def __init__(self, input_file, mesh_file, geometry, rank=0, gdim=3):
        Config.__init__(self, input_file)
        Mesh.__init__(self, mesh_file=mesh_file, geometry=geometry, rank=rank, gdim=gdim)
        FiniteElementSetup.__init__(self)
        Solver.__init__(self)
        ThermalModel.__init__(self)
        MechanicalModel.__init__(self)

    def parameters(self, lhr):
        print(f"----------")
        print(f"parameters")
        print(f"----------")
        self.g = 0.0 # m/s2
        self.lhr = lhr

    def materials(self, **materials):
        print(f"---------")
        print(f"materials")
        print(f"---------")
        self.materials = {}

        for name, mat in materials.items():
            mat["lmbda"] = mat["E"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
            mat["G"] = mat["E"] / (2 * (1 + mat["nu"]))
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
            mu_gamma = self.materials[name]["mu_gamma"]
            q0 = self.materials[name]["gamma_heating"]

            # --- Fissile ---
            if is_fissile:
                q_val = self.lhr / self.area if is_fissile else 0.0
                print(f"  q_third[{name}] = {q_val:.3e} W/m³ (fissile: {is_fissile})")
                print(f"  Heat flux = {self.lhr / self.perimeter:.3e} W/m2 ")

                self.q_third[name] = dolfinx.fem.Constant(self.mesh, q_val)

            # --- Gamma heating ---
            elif q0 > 0.0:
                q_func = dolfinx.fem.Function(self.V_t)
                def f(x):
                    return q0 * np.exp(-x[0] * mu_gamma)
                q_func.interpolate(f)
                q_val = q_func
                print(f"  gamma heating - q_third_0[{name}] = {q0:.3e} W/m³, attenuation: {mu_gamma} 1/m)")
            
                self.q_third[name] = q_func

            # --- No heat source ---
            else:
                q_val = 0.0
                print(f"  q_third_0[{name}] = 0.0 W/m³")

                self.q_third[name] = dolfinx.fem.Constant(self.mesh, q_val)


        # plot q_third(x)
        # import pyvista
        # V = self.V_t
        # cells, types, x = dolfinx.plot.vtk_mesh(V)
        # grid = pyvista.UnstructuredGrid(cells, types, x)
        # grid.point_data["q"] = q_func.x.array.reshape((x.shape[0], 1))
        # grid.set_active_scalars("q")
        # subplotter = pyvista.Plotter(shape=(1, 1))
        # subplotter.subplot(0, 0)
        # subplotter.add_mesh(grid, show_edges=True, show_scalar_bar=True)
        # subplotter.view_xy()
        # subplotter.show()

        # --- Displacement field ---
        self.u = {}

        for name, mat in self.materials.items():
            u_func = dolfinx.fem.Function(self.V_m, name=f"Displacement_{name}")
            u_func.x.array[:] = 0.0
            values = u_func.x.array
            print(f"  {name}: u min={values.min():.2f} m, u max={values.max():.2f} m, u mean={values.mean():.2f} m")
            self.u[name] = u_func

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
        self.stress_mech = {}
        self.stress_th = {}

        for name in self.materials:
            mat = self.materials[name]
            u_local = self.u[name]

            self.strain[name] = self.epsilon(u_local)
            self.stress[name] = self.sigma_mech(u_local, mat) - self.sigma_th(self.T[name], mat)
            self.stress_mech[name] = self.sigma_mech(u_local, mat)
            self.stress_th[name] = - self.sigma_th(self.T[name], mat)
