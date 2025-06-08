import dolfinx
import pyvista
import numpy as np
import matplotlib.cm as cm
from mpi4py import MPI
from petsc4py import PETSc
from core.mesh import Mesh
from core.finiteElementSetup import FiniteElementSetup
from models.thermal_model_cube import ThermalSolver
from models.mechanical_model_cube import MechanicalSolver

class MultiMaterial(Mesh, FiniteElementSetup, ThermalSolver, MechanicalSolver):
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
            
    def boundary_conditions(self):
        self.dofs_bottom1 = self.locateDofs(self.label_map["bottom1"])
        self.dofs_bottom2 = self.locateDofs(self.label_map["bottom2"])

        self.bcs_mechanical = [
            self.zero_displacement(self.mesh, self.dofs_bottom1, self.V_m),
            self.zero_displacement(self.mesh, self.dofs_bottom2, self.V_m)
        ]

        self.bcs_thermal = {
            "fuel": [
                self.fixed_temperature(self.mesh, self.dofs_bottom1, self.V_t, self.conditions["boundary_temperature"])
            ],
            "cladding": [
                self.fixed_temperature(self.mesh, self.dofs_bottom2, self.V_t, self.conditions["boundary_temperature"])
            ]
        }

        self.traction = dolfinx.fem.Constant(self.mesh, self.conditions["coolant_pressure"]) * self.normal

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
        self.T = {}  # dictionary to store temperature functions per material
        # T_fuel = self.T["fuel"]
        # T_cladding = self.T["cladding"]

        for name, mat in self.materials.items():
            T_func = dolfinx.fem.Function(self.V_t, name=f"Temperature_{name}")
            T_func.x.array[:] = mat["T_ref"]
            self.T[name] = T_func

        for name, T_func in self.T.items():
            values = T_func.x.array
            print(f"  [INFO] {name}: min={values.min():.2f} K, max={values.max():.2f} K, mean={values.mean():.2f} K")

        # --- Volumetric power ---
        self.q_third = {}

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
        self.u = {}  # dictionary to store displacement functions per material

        for name, mat in self.materials.items():
            u_func = dolfinx.fem.Function(self.V_m, name=f"Displacement_{name}")
            u_func.x.array[:] = 0.0
            self.u[name] = u_func

        for name, u_func in self.u.items():
            values = u_func.x.array
            print(f"  [INFO] {name}: min={values.min():.2f} m, max={values.max():.2f} m, mean={values.mean():.2f} m")

    def solve(self, max_iters=20, tol_u=1e-6, tol_T=1e-6):
        
        converged = False
        
        for iteration in range(max_iters):

            # Store previous states
            u_old = {name: u.x.array.copy() for name, u in self.u.items()}
            T_old = {name: T.x.array.copy() for name, T in self.T.items()}

            # Solve thermal problem
            for label in ["fuel", "cladding"]:
                if label in self.label_map:
                    tag = self.label_map[label]
                    self.assemble_thermal_problem(label, tag)
                    n_iter = self.solve_thermal_problem(label)
                    print(f"  [INFO] Thermal solver for {label}: {n_iter} iterations")
                    
            # Solve mechanical problem
            # self.assemble_mechanical_solver()
            # self.solve_mechanical_problem()

            # print(f"Max displacement: {np.max(self.u.x.array):2.2f} m")
            # print(f"Max temperature: {np.max(self.T["fuel"].x.array):2.2f} K")
            # print(f"q_third: {self.q_third:2.2f} W/m3")

            # Compute increments
            # du = self.u.x.array - u_old
            dT = {name: self.T[name].x.array - T_old[name] for name in self.T}

            for name, delta in dT.items():
                norm = np.linalg.norm(delta)
                print(f"  [INFO] ΔT[{name}] = {norm:.3e}")

            norm_dT = np.sqrt(sum(np.sum(delta**2) for delta in dT.values()))
            print(f"  [INFO] Global ||ΔT|| = {norm_dT:.3e}")

            # norm_du = np.linalg.norm(du)
            norm_du = 0.0

            print(f"Iteration {iteration}: ||Δu|| = {norm_du:.3e}, ||ΔT|| = {norm_dT:.3e}")

            if norm_dT < tol_T:
                converged = True
                print("  [INFO] Thermal problem converged.")
                break

        if not converged:
            print("  [WARNING]: Solver did not converge within the maximum number of iterations.")

    def get_results(self):
        material_name = "fuel"
        material = self.materials[material_name]

        self.strain = self.epsilon(self.u)
        self.stress = self.sigma_mech(self.u, material) - self.sigma_th(self.T["fuel"], material)

        return self.u, self.strain, self.stress, self.T["fuel"]
