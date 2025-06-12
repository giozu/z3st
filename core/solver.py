import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np
import os

class Solver:
    def __init__(self):
        print("__Solver initializer__")
        self.relax_T = 0.5
        self.relax_u = 0.5

        print("  [INFO] Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")

    def solve(self, max_iters=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6):
        print("\n[INFO] Solver settings:")
        print(f"  → Max iterations        : {max_iters}")
        print(f"  → Coupling              : {self.coupling}")

        if self.coupling == "monolithic":
            self.solve_monolithic(tol=tol_T, max_iter=max_iters)
        elif self.coupling == "staggered":
            self.solve_staggered(max_iters=max_iters, tol_u=tol_u, tol_T=tol_T, tol_mech=tol_mech, tol_th=tol_th)
        else:
            raise ValueError(f"Unknown coupling strategy: {self.coupling}")

    def solve_monolithic(self, tol=1e-8, max_iter=20):
        print("\n[INFO] Starting monolithic thermo-mechanical solve")
        print(f"  → Max iterations        : {max_iter}")
        print(f"  → Tolerance             : {tol}")

        for label in self.materials:
            material = self.materials[label]
            tag = self.label_map[label]

            from basix.ufl import mixed_element

            element_T = self.V_t.ufl_element()
            element_u = self.V_m.ufl_element()

            mixed = mixed_element([element_T, element_u])
            w = dolfinx.fem.functionspace(self.mesh, mixed)

            (T, u) = ufl.split(w)
            (v_T, v_u) = ufl.TestFunctions(w)

            # Measures
            dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)

            # Thermal form
            k = material["k"]
            q_vol = self.q_third[label]
            F_T = (k * ufl.inner(ufl.grad(T), ufl.grad(v_T)) - q_vol * v_T) * dx_local

            # Mechanical form
            epsilon = ufl.sym(ufl.grad(u))
            sigma_m = material["lmbda"] * ufl.tr(epsilon) * ufl.Identity(self.gdim) + 2 * material["mu"] * epsilon
            sigma_th = (3*material["lmbda"] + 2*material["mu"]) * material["alpha"] * (T - material["T_ref"]) * ufl.Identity(self.gdim)
            sigma_tot = sigma_m - sigma_th

            body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -material["rho"] * self.g))
            F_u = (ufl.inner(sigma_tot, ufl.sym(ufl.grad(v_u))) - ufl.dot(body_force, v_u)) * dx_local

            # Total residual
            F = F_T + F_u

            # Assemble nonlinear problem
            problem = dolfinx.fem.petsc.NonlinearProblem(F, w)
            solver = dolfinx.nls.petsc.NewtonSolver(self.mesh.comm, problem)

            solver.convergence_criterion = "incremental"
            solver.rtol = tol
            solver.report = True

            n, converged = solver.solve(w)

            if not converged:
                print(f"[WARNING] Newton solver did not converge for {label}")
            else:
                print(f"[INFO] Monolithic Newton solver converged in {n} iterations")

            # Update fields after solution
            T_sol, u_sol = w.split()

            self.T[label].x.array[:] = T_sol.x.array
            self.u[label].x.array[:] = u_sol.x.array

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
