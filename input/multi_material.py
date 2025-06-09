import dolfinx
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from core.mesh import Mesh
from core.finiteElementSetup import FiniteElementSetup
from models.thermal_model_cube import ThermalSolver
from models.mechanical_model_cube import MechanicalSolver

class MultiMaterial(Mesh, FiniteElementSetup, ThermalSolver, MechanicalSolver):
    def __init__(self, mesh_file, geometry, rank=0, gdim=3):
        Mesh.__init__(self, mesh_file=mesh_file, geometry=geometry, rank=rank, gdim=gdim)
        ThermalSolver.__init__(self)
        MechanicalSolver.__init__(self)
        FiniteElementSetup.__init__(self)

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

        print(f"Thermal boundary conditions")

        self.dofs_bottom1 = self.locateThermalDofs(self.label_map["bottom1"])
        self.dofs_bottom2 = self.locateThermalDofs(self.label_map["bottom2"])
        self.dofs_face2_2 = self.locateThermalDofs(self.label_map["face2_2"])
        
        facets = self.ft.find(self.label_map["bottom1"])
        print(f"[DEBUG] Found {len(facets)} facets for 'bottom1'")

        for mat_type, mat in self.materials.items():
            if mat_type == "fuel":
                dofs = self.dofs_bottom1
                print(f"  [INFO] Assigning thermal BC to '{mat_type}' using dofs_bottom1")
            elif mat_type == "cladding":
                dofs = self.dofs_face2_2
                print(f"  [INFO] Assigning thermal BC to '{mat_type}' using dofs_face2_2")
            else:
                print(f"[WARNING] No thermal boundary condition defined for material type: '{mat_type}'")
                continue

            self.bcs_thermal[mat_type] = [
                self.fixed_temperature(self.mesh, dofs, self.V_t, self.conditions["boundary_temperature"])
            ]

        print(f"Mechanical boundary conditions")

        self.dofs_bottom1 = self.locateMechanicalDofs(self.label_map["bottom1"])
        self.dofs_bottom2 = self.locateMechanicalDofs(self.label_map["bottom2"])

        facets = self.ft.find(self.label_map["bottom1"])
        print(f"[DEBUG] Found {len(facets)} facets for 'bottom1'")

        for mat_type, _ in self.materials.items():
            if mat_type == "fuel":
                dofs = self.dofs_bottom1
                print(f"  [INFO] Assigning mechanical BC to '{mat_type}' using dofs_bottom1")

                # pressure = self.conditions.get("gap_pressure", 0.0)
                pressure = 0.0
                self.traction[mat_type] = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(pressure)) * self.normal
                self.traction_labels[mat_type] = self.label_map["top2"]

                print(f"  [INFO] Assigned traction to '{mat_type}' (fuel pressure = {pressure*1e-6} MPa)")
            elif mat_type == "cladding":
                dofs = self.dofs_bottom2
                print(f"  [INFO] Assigning mechanical BC to '{mat_type}' using dofs_bottom2")
                
                # pressure = self.conditions.get("coolant_pressure", 0.0)
                pressure = 0.0
                self.traction[mat_type] = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(pressure)) * self.normal
                self.traction_labels[mat_type] = self.label_map["top2"]

                print(f"  [INFO] Assigned traction to '{mat_type}' (coolant pressure = {pressure*1e-6} MPa)")
            else:
                print(f"[WARNING] No mechanical boundary condition defined for material type: '{mat_type}'")
                continue

            self.bcs_mechanical[mat_type] = [
                self.zero_displacement(self.mesh, dofs, self.V_m)
            ]

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
            print(f"--- Iteration {iteration} ---")

            # Store previous states
            u_old = {name: u.x.array.copy() for name, u in self.u.items()}
            T_old = {name: T.x.array.copy() for name, T in self.T.items()}

            for label in self.materials:
                if label in self.label_map:
                    tag = self.label_map[label]
                    print(f"\n  [INFO] Assembling thermal problem for '{label}' (tag = {tag})")

                    self.assemble_thermal_problem(label)
                    n_iter_th = self.solve_thermal_problem(label, tol=1e-6)

                    print(f"  [INFO] Thermal solver for '{label}': {n_iter_th+1} iterations")

                    # Residual log
                    residuals = self.thermal_residuals.get(label, [])
                    for i, r in enumerate(residuals):
                        print(f"    [Residual] {label} iter {i+1}: {r:.3e}")

                    with open(f"output/residuals_thermal_{label}.csv", "w") as f:
                        f.write("iteration,residual\n")
                        for i, r in enumerate(residuals):
                            f.write(f"{i+1},{r}\n")

                    print(f"\n  [INFO] Assembling mechanical problem for '{label}' (tag = {tag})")
                    self.assemble_mechanical_solver(label)
                    n_iter_mech = self.solve_mechanical_problem(label, tol=1e-6)

                    print(f"  [INFO] Mechanical solver for '{label}': {n_iter_mech+1} iterations")

                    # Residual log
                    residuals = self.mechanical_residuals.get(label, [])
                    for i, r in enumerate(residuals):
                        print(f"    [Residual] {label} iter {i+1}: {r:.3e}")

                    with open(f"output/residuals_mechanical_{label}.csv", "w") as f:
                        f.write("iteration,residual\n")
                        for i, r in enumerate(residuals):
                            f.write(f"{i+1},{r}\n")

            # print(f"Max displacement: {np.max(self.u.x.array):2.2f} m")
            # print(f"Max temperature: {np.max(self.T["fuel"].x.array):2.2f} K")
            # print(f"q_third: {self.q_third:2.2f} W/m3")

            # --- Norm computation ---
            norm_dT_sq = 0.0
            norm_du_sq = 0.0


            for label in self.materials:
                if label in self.label_map:
                    dT = self.T[label].x.array - T_old[label]
                    du = self.u[label].x.array - u_old[label]

                    tag = self.label_map[label]
                    region_cells = np.where(self.tags.values == tag)[0]
                    self.mesh.topology.create_connectivity(0, self.tdim) # Mapping those vertices back to cells (e.g. when finding DoFs)
                    cell_to_vertex = self.mesh.topology.connectivity(self.tdim, 0)
                    vertex_indices = np.unique(np.hstack([cell_to_vertex.links(i) for i in region_cells]))
                    region_dofs = dolfinx.fem.locate_dofs_topological(self.V_t, 0, vertex_indices)

                    region_dT = dT[region_dofs]
                    region_du = du[region_dofs]
                    region_norm_sq_dT = np.sum(region_dT**2)
                    region_norm_sq_du = np.sum(region_du**2)
                    norm_dT_loc = np.sqrt(region_norm_sq_dT)
                    norm_du_loc = np.sqrt(region_norm_sq_du)

                    print(f"  [INFO] L2 norm : ΔT[{label}] = {norm_dT_loc:.3e}")
                    print(f"  [INFO] L2 norm : Δu[{label}] = {norm_du_loc:.3e}")

                    # Accumulate global squared norm
                    norm_dT_sq += region_norm_sq_dT
                    norm_du_sq += region_norm_sq_du


            norm_dT = np.sqrt(norm_dT_sq)
            norm_du = np.sqrt(norm_du_sq)

            print(f"  [INFO] Global ||ΔT|| = {norm_dT:.3e}")
            print(f"  [INFO] Global ||Δu|| = {norm_du:.3e}")

            if norm_dT < tol_T and norm_du < tol_u:
                converged = True
                print("  [INFO] Thermal/mechanical problem converged.")
                break

        if not converged:
            print("  [WARNING]: Solver did not converge within the maximum number of iterations.")

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


