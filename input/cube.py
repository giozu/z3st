import dolfinx
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from core.mesh import Mesh
from core.finite_element_setup import FiniteElementSetup
from models.thermal_model_cube import ThermalSolver
from models.mechanical_model_cube import MechanicalSolver

class Cube(Mesh, FiniteElementSetup, ThermalSolver, MechanicalSolver):
    def parameters(self, conditions, lhr):
        self.conditions = conditions
        self.g = conditions["gravitational_acceleration"]
        self.lhr = lhr
        self.q_third = lhr / self.area

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
        self.dofs_bottom = self.locateThermalDofs(self.label_map["Bottom"])
        self.dofs_top = self.locateThermalDofs(self.label_map["Top"])

        self.bcs_mechanical = [
            dolfinx.fem.dirichletbc(
                dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m
            )
        ]

        self.bcs_thermal = [
            dolfinx.fem.dirichletbc(
                dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.conditions["boundary_temperature"])), self.dofs_bottom, self.V_t
            )
        ]

        self.traction = dolfinx.fem.Constant(self.mesh, self.conditions["coolant_pressure"]) * self.normal

        # topology, cell_types, _ = dolfinx.plot.vtk_mesh(self.mesh, self.mesh.topology.dim - 1)
        # surface = pyvista.UnstructuredGrid(topology, cell_types, self.mesh.geometry.x)

        # plotter = pyvista.Plotter()
        # unique_labels = np.unique(self.ft.values)

        # print("Face labels present:", unique_labels)
        # for label in unique_labels:
        #     facets = self.ft.find(label)
        #     print(f"Label {label}: {len(facets)} facets")
        # label_map = {
        #         1: ("Bottom", "red"),
        #         2: ("Face 1", "blue"),
        #         3: ("Face 2", "green"),
        #         4: ("Face 3", "orange"),
        #         5: ("Face 4", "purple"),
        #         6: ("Top", "black")
        #     }
        # for label in np.unique(self.ft.values):
        #     facets = self.ft.find(label)
        #     surface_cells = surface.extract_cells(facets)
        #     name, color = label_map.get(label, (f"Face {label}", "white"))
        #     plotter.add_mesh(surface_cells, show_edges=True, opacity=1.0, color=color, label=name)

        # plotter.add_legend()
        # plotter.show()

    def update_state(self):
        self.q_third = self.lhr / self.area

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

        self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
        self.u.x.array[:] = 0.0

    def initialize_geometry(self, atol=1e-6):
        """
        Store original node indices for selected surfaces of the cube.
        """
        coords = self.mesh.geometry.x
        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
        edge = self.edge

        self.idx_bottom = np.where(np.isclose(z, 0.0, atol=atol))[0]
        self.idx_top = np.where(np.isclose(z, edge, atol=atol))[0]
        self.idx_left = np.where(np.isclose(x, 0.0, atol=atol))[0]
        self.idx_right = np.where(np.isclose(x, edge, atol=atol))[0]
        self.idx_front = np.where(np.isclose(y, 0.0, atol=atol))[0]
        self.idx_back = np.where(np.isclose(y, edge, atol=atol))[0]

    def solve(self, max_iters=20, tol_u=1e-6, tol_T=1e-6):
        converged = False
        for iteration in range(max_iters):
            # Store previous states
            u_old = self.u.x.array.copy()
            T_old = self.T["fuel"].x.array.copy()
            # T_old = {name: T.x.array.copy() for name, T in self.T.items()}

            # Solve thermal problem
            self.assemble_thermal_problem()
            self.solve_thermal_problem()

            # Solve mechanical problem
            self.assemble_mechanical_solver()
            self.solve_mechanical_problem()

            print(f"Max displacement: {np.max(self.u.x.array):2.2f} m")
            print(f"Max temperature: {np.max(self.T["fuel"].x.array):2.2f} K")
            print(f"q_third: {self.q_third:2.2f} W/m3")

            # Compute increments
            du = self.u.x.array - u_old
            dT = self.T["fuel"].x.array - T_old

            norm_du = np.linalg.norm(du)
            norm_dT = np.linalg.norm(dT)

            print(f"Iteration {iteration}: ||Δu|| = {norm_du:.3e}, ||ΔT|| = {norm_dT:.3e}")

            if norm_du < tol_u and norm_dT < tol_T:
                converged = True
                print("Converged.")
                break

        if not converged:
            print("Warning: Solver did not converge within the maximum number of iterations.")

    def get_results(self):
        material_name = "fuel"
        material = self.materials[material_name]

        self.strain = self.epsilon(self.u)
        self.stress = self.sigma_mech(self.u, material) - self.sigma_th(self.T["fuel"], material)

        return self.u, self.strain, self.stress, self.T["fuel"]
