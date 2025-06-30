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
        # self.T = {}

        # for name, mat in self.materials.items():
        #     T_func = dolfinx.fem.Function(self.V_t, name=f"Temperature_{name}")
        #     T_func.x.array[:] = mat["T_ref"]
        #     values = T_func.x.array
        #     print(f"  {name}: T min={values.min():.2f} K, T max={values.max():.2f} K, T mean={values.mean():.2f} K")
        #     self.T[name] = T_func
    
        # --- Displacement field ---
        # self.u = {}

        # for name, mat in self.materials.items():
        #     u_func = dolfinx.fem.Function(self.V_m, name=f"Displacement_{name}")
        #     u_func.x.array[:] = 0.0
        #     values = u_func.x.array
        #     print(f"  {name}: u min={values.min():.2f} m, u max={values.max():.2f} m, u mean={values.mean():.2f} m")
        #     self.u[name] = u_func

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

        print("Initializing global temperature field...")
        self.T = dolfinx.fem.Function(self.V_t, name="Temperature")

        for name, mat in self.materials.items():
            print(f"  → Setting initial temperature for material: '{name}'")
            
            domain_tag = self.label_map[name]
            cells = self.tags.find(domain_tag)
            dofs = dolfinx.fem.locate_dofs_topological(self.V_t, self.tdim, cells)
            initial_temp = mat["T_ref"]

            self.T.x.array[dofs] = initial_temp
            
            print(f"    Set {len(dofs)} DOFs to {initial_temp:.2f} K")

        self.T.x.scatter_forward()
        T_vals = self.T.x.array
        print(f"\nGlobal initial T: min={T_vals.min():.2f} K, max={T_vals.max():.2f} K, mean={T_vals.mean():.2f} K")

        print("\nInitializing global displacement field...")
        self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
        self.u.x.array[:] = 0.0
        self.u.x.scatter_forward()
        u_vals = self.u.x.array
        print(f"\nGlobal initial u: min={u_vals.min():.2e} m, max={u_vals.max():.2e} m, mean={u_vals.mean():.2e} m")

    def solve(self, max_iters=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6, tol=1e-5):
        print("\nSolver settings:")
        print(f"  → Coupling              : {self.coupling}")

        if self.coupling == "monolithic":
            self.solve_monolithic(tol=tol, max_iter=max_iters)

        elif self.coupling == "staggered":
            self.solve_staggered(max_iter=max_iters, tol_u=tol_u, tol_T=tol_T, tol_mech=tol_mech, tol_th=tol_th)

        else:
            raise ValueError(f"Unknown coupling strategy: {self.coupling}")

    # def get_results(self):
    #     self.strain = {}
    #     self.stress = {}
    #     self.stress_mech = {}
    #     self.stress_th = {}

    #     for name in self.materials:
    #         mat = self.materials[name]
    #         u_local = self.u[name]

    #         self.strain[name] = self.epsilon(u_local)
    #         self.stress[name] = self.sigma_mech(u_local, mat) - self.sigma_th(self.T[name], mat)
    #         self.stress_mech[name] = self.sigma_mech(u_local, mat)
    #         self.stress_th[name] = - self.sigma_th(self.T[name], mat)

    def get_results(self):
        """
        Computes symbolic strain and stress expressions for each material.
        The results are stored in dictionaries keyed by material name.
        """
        print("Computing symbolic result fields (strain and stress)...")

        epsilon = lambda v: ufl.sym(ufl.grad(v))

        def sigma_mech(v, mat):
            lmbda = dolfinx.fem.Constant(self.mesh, mat["lmbda"])
            G = dolfinx.fem.Constant(self.mesh, mat["G"])
            return lmbda * ufl.tr(epsilon(v)) * ufl.Identity(self.gdim) + 2 * G * epsilon(v)

        def sigma_th(T, mat):
            lmbda = dolfinx.fem.Constant(self.mesh, mat["lmbda"])
            G = dolfinx.fem.Constant(self.mesh, mat["G"])
            alpha = dolfinx.fem.Constant(self.mesh, mat["alpha"])
            T_ref = dolfinx.fem.Constant(self.mesh, mat["T_ref"])
            return (3 * lmbda + 2 * G) * alpha * (T - T_ref) * ufl.Identity(self.gdim)

        self.strain = {}
        self.stress = {}
        self.stress_mech = {}
        self.stress_th = {}

        strain_expression = epsilon(self.u)

        for name, mat in self.materials.items():
            print(f"  → Generating expressions for material: '{name}'")
            
            self.strain[name] = strain_expression
            
            mech_stress_expr = sigma_mech(self.u, mat)
            th_stress_expr = sigma_th(self.T, mat)

            self.stress[name] = mech_stress_expr - th_stress_expr
            self.stress_mech[name] = mech_stress_expr
            self.stress_th[name] = -th_stress_expr

        print("Symbolic result fields generated successfully.")
        
        return True