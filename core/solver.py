import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np
import os
from dolfinx import log
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem

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
        """Solves the coupled thermo-mechanical problem using a monolithic approach.

        This method formulates and solves the fully coupled system of equations
        for thermal conduction and linear elasticity. It assembles the residual
        (the weak form of the PDE system) in a single matrix problem.

        The problem is defined on a mixed function space `self.W`, which, based on
        the `_setup_function_space` method, is ordered as:
        - Sub-space 0: Vector space for Mechanical Displacement.
        - Sub-space 1: Scalar space for Temperature.

        The method correctly assembles the weak form `F` using this ordering.
        The mechanical residual (`F_mechanical`) is associated with the displacement
        test function `v_m` (from sub-space 0), and the thermal residual
        (`F_thermal`) is associated with the temperature test function `v_t`
        (from sub-space 1).

        The resulting nonlinear system is solved using a Newton-Raphson method
        provided by PETSc. After convergence, the individual solution fields
        (Displacement and Temperature) are extracted from the mixed solution vector.

        Note:
            The original implementation had a critical bug in the solution extraction
            step, where the sub-spaces were incorrectly indexed, leading to a swap
            of the displacement and temperature fields. This version corrects that bug.
            Also, the `tol` and `max_iter` arguments are not passed to the solver,
            which uses hardcoded values instead.

        Args:
            tol (float, optional): Desired tolerance for the solver. Defaults to 1e-8.
                                (Currently unused in implementation).
            max_iter (int, optional): Maximum number of iterations for the Newton
                                    solver. Defaults to 50. (Currently unused).

        Side Effects:
            - self.sol_mixed (dolfinx.fem.Function): Updated in-place with the
            computed solution vector for the mixed problem.
            - self.T (dict): The dictionary of temperature Functions is updated
            by interpolating the final temperature solution.
            - self.u (dict): The dictionary of displacement Functions is updated
            by interpolating the final displacement solution.
        """
        from dolfinx.fem.petsc import NonlinearProblem
        from dolfinx.nls.petsc import NewtonSolver
        from petsc4py import PETSc
        import ufl

        print("\n[INFO] Starting monolithic thermo-mechanical solve")
        print(f"  → Tolerance              : {tol}")

        # Extract unknowns from mixed solution. Based on the definition of self.W:
        # u_m is the function for Displacement (from sub-space 0)
        # u_t is the function for Temperature (from sub-space 1)
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

            # Thermal residual (associated with v_t, from sub-space 1)
            F_thermal = (k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) - q_vol * v_t) * dx_local

            if self.on["thermal"] is True:
                print("Thermal model ON")
                F += F_thermal

            # Mechanical residual
            epsilon = ufl.sym(ufl.grad(u_m))
            sigma_m = lmbda * ufl.tr(epsilon) * ufl.Identity(self.tdim) + 2 * mu * epsilon
            sigma_th = (3 * lmbda + 2 * mu) * alpha * delta_T * ufl.Identity(self.gdim)
            sigma_tot = sigma_m - sigma_th

            if self.on["mechanical"] is True:
                print("Mechanical model ON")
                F += (ufl.inner(sigma_tot, ufl.sym(ufl.grad(v_m))) - ufl.dot(body_force, v_m)) * dx_local

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

        ksp = solver.krylov_solver
        opts = PETSc.Options()
        option_prefix = ksp.getOptionsPrefix()

        # Direct solver (LU)
        print("ksp_type, preonly")
        opts[f"{option_prefix}ksp_type"] = "preonly"
        opts[f"{option_prefix}pc_type"] = "lu"
        opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
        opts[f"{option_prefix}ksp_monitor"] = None

        # snes_prefix = solver.getOptionsPrefix()
        # opts[f"{snes_prefix}snes_linesearch_type"] = "basic"
        # opts[f"{snes_prefix}snes_monitor"] = None

        # GAMG
        # print("pc_gamg_type, agg")
        # opts[f"{option_prefix}pc_gamg_type"] = "agg"
        # opts[f"{option_prefix}pc_gamg_coarse_eq_limit"] = 1000
        # opts[f"{option_prefix}mg_levels_ksp_type"] = "chebyshev"
        # opts[f"{option_prefix}mg_levels_pc_type"] = "sor"
        # opts[f"{option_prefix}mg_levels_esteig_ksp_type"] = "cg"
        # opts[f"{option_prefix}mg_levels_ksp_chebyshev_esteig_steps"] = 20
        # opts[f"{option_prefix}ksp_monitor"] = None

        ksp.setFromOptions()

        # Set Newton solver options
        solver.atol = tol
        solver.rtol = tol
        solver.convergence_criterion = "incremental"
        solver.max_it = max_iter

        # Solve
        print("[INFO] Solving monolithic nonlinear problem...")
        n, converged = solver.solve(self.sol_mixed)

        if converged:
            print(f"[INFO] Newton solver converged in {n} iterations")
        else:
            print(f"[WARNING] Newton solver did not converge after {n} iterations")

        # --- CORRECTED SOLUTION EXTRACTION ---
        # Sub-space 0 corresponds to Displacement (u_m)
        # Sub-space 1 corresponds to Temperature (u_t)
        u_sol = self.sol_mixed.sub(0).collapse()
        u_sol.name = "Displacement"
        T_sol = self.sol_mixed.sub(1).collapse()
        T_sol.name = "Temperature"

        print(f"Min/Max temperature: {T_sol.x.array.min():.2f} / {T_sol.x.array.max():.2f}")
        print(f"Min/Max displacement magnitude: {np.linalg.norm(u_sol.x.array.reshape(-1, self.gdim), axis=1).min():.2e} / {np.linalg.norm(u_sol.x.array.reshape(-1, self.gdim), axis=1).max():.2e}")

        # Update material fields for post-processing or further steps
        for name in self.materials.keys():
            self.T[name].interpolate(T_sol)
            self.u[name].interpolate(u_sol)
    def solve_staggered(self, max_iter=10, tol_u=1e-6, tol_T=1e-6, tol_mech=1e-8, tol_th=1e-8):
        """
        Solves the coupled thermo-mechanical problem using a staggered (partitioned)
        approach with individual tolerances.

        Args:
            max_iter (int): Maximum number of staggering iterations between physics.
            tol_u (float): Convergence tolerance for the displacement field norm
                           between staggering iterations.
            tol_T (float): Convergence tolerance for the temperature field norm
                           between staggering iterations.
            tol_mech (float): Relative tolerance for the linear mechanical solver (KSP).
            tol_th (float): Relative tolerance for the linear thermal solver (KSP).
        """
        print("\n[INFO] Starting staggered thermo-mechanical solve")

        # --- 1. Preparazione dei campi e delle BC ---
        T_sol = self.T['fuel']
        u_sol = self.u['fuel']
        
        T_old = dolfinx.fem.Function(T_sol.function_space)
        u_old = dolfinx.fem.Function(u_sol.function_space)
        T_old.x.array[:] = T_sol.x.array
        u_old.x.array[:] = u_sol.x.array
        
        bcs_t = [bc for bc in self.bcs if bc.function_space == self.V_t]
        bcs_m = [bc for bc in self.bcs if bc.function_space.contains(self.V_m)]

        print(f"[INFO] Found {len(bcs_t)} thermal BCs and {len(bcs_m)} mechanical BCs.")
        print(f"[INFO] Staggering tolerances: dT={tol_T:.1e}, du={tol_u:.1e}")
        print(f"[INFO] Linear solver tolerances: thermal={tol_th:.1e}, mechanical={tol_mech:.1e}")

        # --- 2. Ciclo di Staggering ---

        for i in range(max_iter):
            print(f"--- Staggered Iteration {i+1}/{max_iter} ---")

            # === A) SOTTO-PROBLEMA TERMICO ===
            print("[INFO] Solving thermal sub-problem...")
            u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
            a_t, L_t = 0, 0
            # ... (assemblaggio di a_t e L_t come prima) ...
            for label in self.materials:
                # ... (codice di assemblaggio per materiale) ...
                material = self.materials[label]
                tag = self.label_map[label]
                dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)
                k = dolfinx.fem.Constant(self.mesh, material["k"])
                q_vol = self.q_third[label]
                a_t += k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx_local
                L_t += q_vol * v_t * dx_local
            if 'neumann_thermal' in self.bcs_info:
                # ... (codice di assemblaggio per Neumann BC) ...
                for bc_info in self.bcs_info['neumann_thermal']:
                    ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                    L_t += bc_info['value'] * v_t * ds_local

            # Risolvi il problema lineare termico con la tolleranza specificata
            petsc_options_th = {"ksp_type": "cg", "pc_type": "gamg", "ksp_rtol": tol_th}
            problem_t = LinearProblem(a_t, L_t, bcs=bcs_t, u=T_sol, petsc_options=petsc_options_th)
            problem_t.solve()
            
            # === B) SOTTO-PROBLEMA MECCANICO ===
            print("[INFO] Solving mechanical sub-problem...")
            u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
            a_m, L_m = 0, 0
            # ... (assemblaggio di a_m e L_m come prima, usando T_sol) ...
            for label in self.materials:
                # ... (codice di assemblaggio per materiale) ...
                material = self.materials[label]
                tag = self.label_map[label]
                dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)
                lmbda, mu = dolfinx.fem.Constant(self.mesh, material["lmbda"]), dolfinx.fem.Constant(self.mesh, material["mu"])
                alpha, T_ref = dolfinx.fem.Constant(self.mesh, material["alpha"]), dolfinx.fem.Constant(self.mesh, material["T_ref"])
                rho, g = dolfinx.default_scalar_type(material["rho"]), dolfinx.default_scalar_type(self.g)
                body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))
                epsilon = lambda u: ufl.sym(ufl.grad(u))
                sigma_m = lambda u: lmbda * ufl.tr(epsilon(u)) * ufl.Identity(self.tdim) + 2 * mu * epsilon(u)
                delta_T = T_sol - T_ref
                sigma_th = (3 * lmbda + 2 * mu) * alpha * delta_T * ufl.Identity(self.gdim)
                a_m += ufl.inner(sigma_m(u_m), epsilon(v_m)) * dx_local
                L_m += ufl.dot(body_force, v_m) * dx_local + ufl.inner(sigma_th, epsilon(v_m)) * dx_local
            if 'traction' in self.bcs_info:
                 # ... (codice di assemblaggio per Traction BC) ...
                for bc_info in self.bcs_info['traction']:
                    ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                    traction_vector = bc_info['value'] * ufl.FacetNormal(self.mesh)
                    L_m += ufl.dot(traction_vector, v_m) * ds_local

            # Risolvi il problema lineare meccanico con la tolleranza specificata
            petsc_options_mech = {"ksp_type": "cg", "pc_type": "gamg", "ksp_rtol": tol_mech, "pc_gamg_coarse_eq_limit": 1000}
            problem_m = LinearProblem(a_m, L_m, bcs=bcs_m, u=u_sol, petsc_options=petsc_options_mech)
            problem_m.solve()
            
            # === C) CONTROLLO CONVERGENZA ===
            # Calcola l'errore relativo tra le iterazioni
            # Usiamo dolfinx.la.Norm per calcoli in parallelo corretti
            comm = self.mesh.comm
            T_error = Norm.L2(comm, T_sol.x.array - T_old.x.array)
            T_norm = Norm.L2(comm, T_sol.x.array)
            u_error = Norm.L2(comm, u_sol.x.array - u_old.x.array)
            u_norm = Norm.L2(comm, u_sol.x.array)
            
            # Evita divisione per zero se la soluzione è nulla
            T_rel_error = T_error / T_norm if T_norm > 1e-12 else T_error
            u_rel_error = u_error / u_norm if u_norm > 1e-12 else u_error
            
            print(f"  Convergence check: rel_dT = {T_rel_error:.2e}, rel_du = {u_rel_error:.2e}")

            # Aggiorna i campi vecchi
            T_old.x.array[:] = T_sol.x.array
            u_old.x.array[:] = u_sol.x.array
            
            if T_rel_error < tol_T and u_rel_error < tol_u:
                print(f"[INFO] Staggered scheme converged in {i+1} iterations.")
                break
            
            # Per il tuo attuale problema a senso unico, questo è comunque vero
            if i == 0:
                 print("[INFO] One-way coupling problem: convergence is guaranteed after 1 iteration.")
                 break
        else: # Questo blocco viene eseguito se il loop 'for' termina senza 'break'
            print(f"[WARNING] Staggered scheme did not converge after {max_iter} iterations.")
            
        # --- 3. Post-processing ---
        print(f"\n[INFO] Staggered solve finished.")
        print(f"Min/Max temperature: {T_sol.x.array.min():.2f} / {T_sol.x.array.max():.2f}")
        
        u_magnitude = np.sqrt(np.sum(u_sol.x.array.reshape(-1, self.gdim)**2, axis=1))
        print(f"Min/Max displacement magnitude: {u_magnitude.min():.2e} / {u_magnitude.max():.2e}")
        
        # Aggiorna tutti i campi nei dizionari self.T e self.u
        for name in self.materials.keys():
            self.T[name].interpolate(T_sol)
            self.u[name].interpolate(u_sol)

    def solve_staggered_old(self, max_iter=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6):
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
