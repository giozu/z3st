import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np
import os
from dolfinx import log
from mpi4py import MPI

class Solver:
    def __init__(self):
        print("__Solver initializer__")
        self.relax_T = 0.5
        self.relax_u = 0.5

        print("  [INFO] Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")

        self.bcs_info = {}    

    def solve_monolithic(self, tol=1e-8, max_iter=50):

        from dolfinx.fem.petsc import NonlinearProblem
        from dolfinx.nls.petsc import NewtonSolver
        from petsc4py import PETSc
        import ufl

        print("\n[INFO] Starting monolithic thermo-mechanical solve")

        # Extract unknowns from mixed solution
        (u_m, u_t) = ufl.split(self.sol_mixed)
        (v_m, v_t) = ufl.TestFunctions(self.W)

        F = 0  # Start assembling residual

        for label in self.materials:
            material = self.materials[label]
            tag = self.label_map[label]
            dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)

            # Material properties
            k = dolfinx.fem.Constant(self.mesh, material["k"])
            lmbda = dolfinx.fem.Constant(self.mesh, material["lmbda"])
            mu = dolfinx.fem.Constant(self.mesh, material["mu"])
            alpha = dolfinx.fem.Constant(self.mesh, material["alpha"])
            T_ref = dolfinx.fem.Constant(self.mesh, material["T_ref"])
            rho = dolfinx.default_scalar_type(material["rho"])
            g = dolfinx.default_scalar_type(self.g)
            body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))
            q_vol = self.q_third[label]

            delta_T = u_t - T_ref

            # Thermal residual
            F_thermal = (k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) - q_vol * v_t) * dx_local

            # Mechanical residual
            epsilon = ufl.sym(ufl.grad(u_m))
            sigma_m = lmbda * ufl.tr(epsilon) * ufl.Identity(self.tdim) + 2 * mu * epsilon
            sigma_th = (3 * lmbda + 2 * mu) * alpha * delta_T * ufl.Identity(self.gdim)
            sigma_tot = sigma_m - sigma_th

            F_mechanical = (ufl.inner(sigma_tot, ufl.sym(ufl.grad(v_m))) - ufl.dot(body_force, v_m)) * dx_local

            F += F_thermal # + F_mechanical

        # Neumann and traction BCs
        n_vec = ufl.FacetNormal(self.mesh)

        if 'neumann_thermal' in self.bcs_info:
            for bc_info in self.bcs_info['neumann_thermal']:
                ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                F -= bc_info['value'] * v_t * ds_local

        if 'traction' in self.bcs_info:
            for bc_info in self.bcs_info['traction']:
                ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                traction_vector = bc_info['value'] * n_vec
                F -= ufl.dot(traction_vector, v_m) * ds_local

        # Define nonlinear problem and solver
        problem = NonlinearProblem(F, self.sol_mixed, bcs=self.bcs)
        solver = NewtonSolver(self.mesh.comm, problem)

        # Set Newton solver options
        solver.atol = 1e-8
        solver.rtol = 1e-8
        solver.convergence_criterion = "incremental"
    
        # Solve
        print("[INFO] Solving monolithic nonlinear problem...")
        n, converged = solver.solve(self.sol_mixed)

        if converged:
            print(f"[INFO] Newton solver converged in {n} iterations")
        else:
            print("[WARNING] Newton solver did not converge")

        # Extract solutions
        T_sol = self.sol_mixed.sub(0).collapse()
        T_sol.name = "Temperature"
        u_sol = self.sol_mixed.sub(1).collapse()
        u_sol.name = "Displacement"

        print(f"Min/Max temperature: {T_sol.x.array.min():.2f} / {T_sol.x.array.max():.2f}")
        print(f"Min/Max displacement: {u_sol.x.array.min():.2e} / {u_sol.x.array.max():.2e}")

        # Update material fields
        for name in self.materials.keys():
            self.T[name].interpolate(T_sol)
            self.u[name].interpolate(u_sol)


    def solve_staggered(self, max_iter=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6):
            print("[INFO] Starting staggered thermo-mechanical solver")
            print(f"  → Max iterations        : {max_iter}")
            print(f"  → Global tolerance ΔT   : {tol_T}")
            print(f"  → Global tolerance Δu   : {tol_u}")
            print(f"  → Thermal solver tol    : {tol_th}")
            print(f"  → Mechanical solver tol : {tol_mech}")

            converged = False
            
            self.i_th = {name: 0 for name in self.materials}
            self.i_mech = {name: 0 for name in self.materials}
            
            for iteration in range(max_iter):
                print(f"--- Iteration {iteration} ---")

                # Store previous states
                self.u_old = {name: u.x.array.copy() for name, u in self.u.items()}
                self.T_old = {name: T.x.array.copy() for name, T in self.T.items()}

                for label in self.materials:
                    if label in self.label_map:
                        tag = self.label_map[label]

                        print(f"\n  [INFO] Assembling thermal problem for '{label}' (tag = {tag})")

                        if self.on["thermal"]:
                            self.assemble_thermal_problem(label)
                            n_iter_th = self.solve_thermal_problem(label, tol=tol_th)
                        else:
                            n_iter_th = 0

                        print(f"  [INFO] Thermal solver for '{label}': {n_iter_th+1} iterations")

                        # Relaxation
                        self.T[label].x.array[:] = self.relax_T * self.T[label].x.array[:] + (1.0 - self.relax_T) * self.T_old[label]

                        # Heat flux calculation
                        k = self.materials[label]["k"]
                        grad_T = ufl.grad(self.T[label])
                        q_vec = [-k * grad_T[i] for i in range(self.mesh.geometry.dim)]  # [qx, qy, qz]

                        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=self.label_map[label])
                        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1.0 * dx_subdomain))

                        print(f"\nAverage heat flux components in '{label}':")
                        for i, q_comp in enumerate(q_vec):
                            integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_comp * dx_subdomain))
                            q_avg = integral / volume
                            print(f"  q_{'xyz'[i]} = {q_avg:.4e} W/m²")

                        # Residual log
                        residuals = self.thermal_residuals.get(label, [])
                        if self.debug_mode is True:
                            for i, r in enumerate(residuals):
                                print(f"    [Residual] {label} iter {i+1}: {r:.3e}")

                        filename = f"output/residuals_thermal_{label}.csv"
                        file_exists = os.path.isfile(filename)

                        with open(filename, "a") as f:
                            if not file_exists:
                                f.write("iteration,residual\n")
                            for r in residuals:
                                self.i_th[label] += 1
                                f.write(f"{self.i_th[label]},{r}\n")

                        print(f"\n  [INFO] Assembling mechanical problem for '{label}' (tag = {tag})")

                        if self.on["mechanical"]:
                            self.assemble_mechanical_solver(label)
                            n_iter_mech = self.solve_mechanical_problem(label, tol=tol_mech)
                        else:
                            n_iter_mech = 0

                        print(f"  [INFO] Mechanical solver for '{label}': {n_iter_mech+1} iterations")

                        # Relaxation
                        self.u[label].x.array[:] = self.relax_u * self.u[label].x.array[:] + (1.0 - self.relax_u) * self.u_old[label]

                        # Residual log
                        residuals = self.mechanical_residuals.get(label, [])
                        if self.debug_mode is True:
                            for i, r in enumerate(residuals):
                                print(f"    [Residual] {label} iter {i+1}: {r:.3e}")

                        filename = f"output/residuals_mechanical_{label}.csv"
                        file_exists = os.path.isfile(filename)

                        with open(filename, "a") as f:
                            if not file_exists:
                                f.write("iteration,residual\n")
                            for r in residuals:
                                self.i_mech[label] += 1
                                f.write(f"{self.i_mech[label]},{r}\n")
        
                # print(f"Max displacement: {np.max(self.u.x.array):2.2f} m")
                # print(f"Max temperature: {np.max(self.T["fuel"].x.array):2.2f} K")
                # print(f"q_third: {self.q_third:2.2f} W/m3")

                # --- Norm computation ---
                norm_dT_sq = 0.0
                norm_du_sq = 0.0

                for label in self.materials:
                    if label in self.label_map:
                        dT = self.T[label].x.array - self.T_old[label]
                        du = self.u[label].x.array - self.u_old[label]

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
                    print(f"  [INFO] Thermal/mechanical problem converged in {iteration+1} iterations.")
                    break

            if not converged:
                print("  [WARNING]: Solver did not converge within the maximum number of iterations.")
