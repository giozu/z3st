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

    def parameters(self, conditions, lhr):
        self.conditions = conditions
        self.g = conditions["gravitational_acceleration"]
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
        with open("conditions/boundary_conditions.yaml", "r") as f:
            print(f"Loading boundary conditions from 'conditions/boundary_conditions.yaml'")
            self.boundary_conditions = yaml.safe_load(f)

        self.thermal_dirichlet_bcs = self.set_thermal_boundary_conditions()
        self.mechanical_dirichlet_bcs = self.set_mechanical_boundary_conditions()

        self.bcs = self.thermal_dirichlet_bcs + self.mechanical_dirichlet_bcs
    
    def initialize_fields(self):
        # --- Temperature field ---
        # v_t: Test Function
        # - Symbolic function defined as v_t = TestFunction(V_t)
        # - Represents virtual temperature variations in the weak form
        # - Used in assembling the thermal right-hand side
        # - Does not hold any actual values
        # u_t: Trial Function
        # - Symbolic unknown defined as u_t = TrialFunction(V_t)
        # - Represents the temperature field being solved for
        # - Used to assemble the thermal stiffness matrix
        # - Not used directly for numerical evaluation
        # T: Function
        # - Created with T = Function(V_t, name="Temperature")
        # - Contains the actual computed temperature field
        # - Updated by the solver during the thermal solution process
        # - Access values via T.x.array
        # - This is the field to use for output, visualization, or coupling to mechanics
        # Important: 'T' holds the real temperature solution.
        # 'u_t' and 'v_t' are symbolic and only used to define the variational forms.

        print(f"Initializing temperature fields for: {list(self.materials.keys())}")
        for name, mat in self.materials.items():
            T_func = dolfinx.fem.Function(self.V_t, name=f"Temperature_{name}")
            T_func.x.array[:] = mat["T_ref"]
            self.T[name] = T_func

        for name, T_func in self.T.items():
            values = T_func.x.array
            print(f"  [INFO] {name}: min={values.min():.2f} K, max={values.max():.2f} K, mean={values.mean():.2f} K")

        # for name, T_func in self.T.items():
        #     tag = self.label_map[name]
        #     region_cells = np.where(self.tags.values == tag)[0]
        #     self.mesh.topology.create_connectivity(self.mesh.topology.dim, 0)
        #     cell_to_vertex = self.mesh.topology.connectivity(self.mesh.topology.dim, 0)
        #     vertex_indices = np.unique(np.hstack([cell_to_vertex.links(i) for i in region_cells]))
        #     self.mesh.topology.create_connectivity(0, self.mesh.topology.dim)
        #     region_dofs = dolfinx.fem.locate_dofs_topological(self.V_t, 0, vertex_indices)
            
        #     values = T_func.x.array[region_dofs]
        #     print(f"[INFO] {name}: min T over region = {values.min():.2f} K, max = {values.max():.2f} K, mean = {values.mean():.2f} K")

        for name, mat in self.materials.items():
            is_fissile = bool(mat.get("fissile", False))
            self.q_third[name] = self.lhr / self.area if is_fissile else 0.0
            print(f"  [INFO] q_third[{name}] = {self.q_third[name]:.2f} W/m³ (fissile: {is_fissile})")

        # --- Displacement ---
        # v_m: Test Function
        # - Defined symbolically using v_m = TestFunction(V_m)
        # - Represents virtual displacements in the weak form
        # - Used when assembling the right-hand side of the variational equation
        # - Does not contain solution data
        # u_m: Trial Function
        # - Defined symbolically using u_m = TrialFunction(V_m)
        # - Represents the unknown displacement field in the variational problem
        # - Used in assembling the system matrix (stiffness matrix)
        # - Not assigned or modified during solution
        # u: Function
        # - Created with u = Function(V_m, name="Displacement")
        # - Stores the actual numerical solution of the displacement field
        # - Values are updated during solve(): solver.solve(b, u.x.petsc_vec)
        # - This is the field to use in postprocessing (e.g., stress evaluation)
        # Important: Only 'u' contains the computed displacement.
        # 'u_m' and 'v_m' are symbolic and used solely to define the variational problem

        print(f"Initializing displacement fields for: {list(self.materials.keys())}")
        for name, mat in self.materials.items():
            u_func = dolfinx.fem.Function(self.V_m, name=f"Displacement_{name}")
            u_func.x.array[:] = 0.0
            self.u[name] = u_func

        for name, u_func in self.u.items():
            values = u_func.x.array
            print(f"  [INFO] {name}: min={values.min():.2f} m, max={values.max():.2f} m, mean={values.mean():.2f} m")

        # --- Displacement / temperature ---
        w_sub_T = self.w.sub(0)
        w_sub_u = self.w.sub(1)
        w_sub_T.x.array[:] = mat["T_ref"]
        w_sub_u.x.array[:] = 0.0
        
    def solve(self, max_iters=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6, tol=1e-5):
        print("\n[INFO] Solver settings:")
        print(f"  → Max iterations        : {max_iters}")
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
        
        for name in self.u:
            u = self.u[name]
            T = self.T[name]
            mat = self.materials[name]

            self.strain[name] = self.epsilon(u)
            self.stress[name] = self.sigma_mech(u, mat) - self.sigma_th(T, mat)

        return self.u, self.strain, self.stress, self.T


