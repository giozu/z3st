# --- Python modules ---
import dolfinx
import yaml
from mpi4py import MPI
import numpy as np
import os
from petsc4py import PETSc
import ufl
import importlib

# --- ZEST modules ---
from core.config import Config
from core.mesh import Mesh
from core.solver import Solver
from core.finite_element_setup import FiniteElementSetup
from models.thermal_model import ThermalModel
from models.mechanical_model import MechanicalModel
from models.gap_model import GapModel

class Spine(Config, Mesh, FiniteElementSetup, Solver, ThermalModel, MechanicalModel, GapModel):

    def __init__(self, input_file, mesh_file, geometry, rank=0, gdim=3):
        Config.__init__(self, input_file)
        Mesh.__init__(self, mesh_file=mesh_file, geometry=geometry, rank=rank, gdim=gdim)
        FiniteElementSetup.__init__(self)
        Solver.__init__(self)
        ThermalModel.__init__(self)
        MechanicalModel.__init__(self)
        GapModel.__init__(self)

    def parameters(self, lhr):
        print(f"----------")
        print(f"parameters")
        print(f"----------")
        self.g = 0.0 # m/s2
        self.lhr = lhr

    def resolve_function(self, path: str):
        module_path, func_name = path.rsplit(".", 1)
        mod = importlib.import_module(module_path)
        return getattr(mod, func_name)

    def materials(self, **materials):
        print(f"---------")
        print(f"materials")
        print(f"---------")
        self.materials = {}

        for name, mat in materials.items():
            print(f"{name.capitalize()} material loaded:")

            if "E" in mat and "nu" in mat:
                mat["E"] = float(mat["E"])
                mat["nu"] = float(mat["nu"])

                mat["lmbda"] = mat["E"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
                mat["G"] = mat["E"] / (2 * (1 + mat["nu"]))
                mat["bulk_modulus"] = mat["E"] / (3 * (1 - 2 * mat["nu"]))

                # Thermal conductivity
                if "k" in mat:
                    if isinstance(mat["k"], str):
                        print(f"  → k defined as symbolic function: {mat['k']}")
                        k_func = self.resolve_function(mat["k"])
                        bu = mat.get("burnup", 0.0)
                        gad = mat.get("gadolinia", 0.0)
                        mat["_k_func"] = (k_func, bu, gad)
                    else:
                        print(f"  → k defined as constant: {mat['k']}")
                else:
                    print(f"  → k not defined for {name}")

            else:
                print(f"  [INFO] '{name}' has no elasticity parameters — skipping mechanical properties.")

            self.materials[name] = mat

            for k in sorted(mat.keys()):
                v = mat[k]
                print(f"  {k:<15} → {v} ({type(v).__name__})")

            print("\n")

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

        self.q_third = dolfinx.fem.Function(self.V_t, name="q_third")
        self.q_third.x.array[:] = 0.0

        for name, mat in self.materials.items():
            dofs = self.locateDomainDofs(self.label_map, name, self.V_t)

            if mat.get("fissile", False): # if mat["fissile"] == True.. if "fissile" is not there, False is given, block skipped
                print(f"Fissile material")
                q_val = self.lhr / self.area
                self.q_third.x.array[dofs] = q_val
                print(f"  q_third = {q_val:.3e} W/m³ (fissile: {mat.get("fissile", False)})")
                print(f"  Heat flux = {self.lhr / self.perimeter:.3e} W/m2 ")
            
            if float(mat["gamma_heating"]) > 0.0:
                q_third_0 = float(mat["gamma_heating"])
                def f(x):
                    return q_third_0 * np.exp(-x[0] * mat["mu_gamma"])
                f_func = dolfinx.fem.Function(self.V_t)
                f_func.interpolate(f)
                self.q_third.x.array[dofs] = f_func.x.array[dofs]
                print(f"  gamma heating - q_third_0 = {q_third_0:.3e} W/m³, attenuation coefficient: {mat["mu_gamma"]} 1/m")

        print("\nInitializing global temperature field...")
        self.T = dolfinx.fem.Function(self.V_t, name="Temperature")
        for name, mat in self.materials.items():
            print(f"  → Setting initial temperature for material: '{name}'")
            dofs = self.locateDomainDofs(self.label_map, name, self.V_t)
            self.T.x.array[dofs] = mat["T_ref"]
            print(f"    Set {len(dofs)} DOFs to {mat["T_ref"]:.2f} K")

        self.T.x.scatter_forward()
        T_vals = self.T.x.array
        print(f"  Initial T: min={T_vals.min():.2f} K, max={T_vals.max():.2f} K, mean={T_vals.mean():.2f} K")

        print("\nInitializing global displacement field...")
        self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
        self.u.x.array[:] = 0.0
        self.u.x.scatter_forward()
        u_vals = self.u.x.array
        print(f"  Initial u: min={u_vals.min():.2e} m, max={u_vals.max():.2e} m, mean={u_vals.mean():.2e} m")

        # Thermal conductivity
        for name, mat in self.materials.items():
            if "_k_func" in mat:
                k_func, bu, gad = mat["_k_func"]
                mat["k"] = k_func(self.T, bu=bu, gad=gad)
                print("\nk expression for", name, "→", mat["k"])

    def solve(self, max_iters=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6, tol=1e-5):
        print("\nSolver settings:")
        print(f"  → Coupling              : {self.coupling}")

        if self.coupling == "monolithic":
            self.solve_monolithic(tol=tol, max_iter=max_iters)

        # if self.coupling == "monolithic":
        #     raise NotImplementedError("Monolithic coupling is not yet supported; use 'staggered'.")

        elif self.coupling == "staggered":
            self.solve_staggered(max_iter=max_iters, tol_u=tol_u, tol_T=tol_T, tol_mech=tol_mech, tol_th=tol_th)

        else:
            raise ValueError(f"Unknown coupling strategy: {self.coupling}")

    def get_results(self):
        """
        Computes symbolic strain and stress expressions for each material.
        The results are stored in dictionaries keyed by material name.
        """

        if self.coupling == "staggered":
            print("Computing symbolic result fields (strain and stress)...")

            self.stress = {}
            self.stress_mech = {}
            self.stress_th = {}

            self.strain = self.epsilon(self.u)

            for name, mat in self.materials.items():
                print(f"  → Generating expressions for material: '{name}'")
                
                mech_stress_expr = self.sigma_mech(self.u, mat)
                th_stress_expr = self.sigma_th(self.T, mat)

                self.stress[name] = mech_stress_expr - th_stress_expr
                self.stress_mech[name] = mech_stress_expr
                self.stress_th[name] = -th_stress_expr

            print("Symbolic result fields generated successfully.")
        
        return True