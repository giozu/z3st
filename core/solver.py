import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np
import os
from dolfinx import log
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.la import Norm 

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

    def solve_staggered(self, max_iter=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6):
        print("[INFO] Starting staggered thermo-mechanical solver")
        print(f"  → Max iterations        : {max_iter}")
        print(f"  → Staggering tolerance |ΔT| : {tol_T:.1e}")
        print(f"  → Staggering tolerance |Δu| : {tol_u:.1e}")
        print(f"  → Thermal KSP tolerance   : {tol_th:.1e}")
        print(f"  → Mechanical KSP tolerance: {tol_mech:.1e}")

        # --- Setup Iniziale ---
        comm = self.mesh.comm
        T_sol = self.T['fuel']  # Funzione di riferimento per la soluzione
        u_sol = self.u['fuel']  # Funzione di riferimento per la soluzione

        # Funzioni per salvare lo stato precedente
        T_old = dolfinx.fem.Function(T_sol.function_space)
        u_old = dolfinx.fem.Function(u_sol.function_space)

        # Contatori per il logging dei residui
        self.i_th = {name: 0 for name in self.materials}
        self.i_mech = {name: 0 for name in self.materials}

        W_m_sub, W_t_sub = self.W.sub(0), self.W.sub(1) 

        # Estrai le BC per ogni fisica (una sola volta)
        bcs_t = []
        bcs_m = []
        for bc in self.bcs:
            # bc.function_space è lo spazio su cui la BC è definita (probabilmente self.W)
            # W_t_sub.collapse() restituisce lo spazio funzionale "collassato" che possiamo confrontare
            # per vedere se la BC agisce su quella parte del problema misto.
            if bc.function_space == W_t_sub.collapse()[0]:
                bcs_t.append(bc)
            elif bc.function_space == W_m_sub.collapse()[0]:
                bcs_m.append(bc)

        # --- Ciclo di Staggering ---
        for iteration in range(max_iter):
            print(f"\n--- Staggering Iteration {iteration+1}/{max_iter} ---")

            # Salva lo stato precedente per il controllo di convergenza e il rilassamento
            T_old.x.array[:] = T_sol.x.array
            u_old.x.array[:] = u_sol.x.array

            # =========================================================
            #  1. RISOLUZIONE DEL PROBLEMA TERMICO (globale)
            # =========================================================
            if self.on["thermal"]:
                print("[INFO] Assembling and solving GLOBAL thermal problem...")
                
                # --- INIZIO MODIFICA ---
                u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
                # Inizializza a 0. Le forme verranno popolate nel ciclo.
                a_t = 0
                L_t = 0
                # --- FINE MODIFICA ---

                # Assembla le forme variazionali su tutti i materiali
                for label, material in self.materials.items():
                    tag = self.label_map[label]
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)
                    k = dolfinx.fem.Constant(self.mesh, material["k"])
                    q_vol = self.q_third[label]
                    
                    # Ora usiamo ufl.inner per la sorgente volumetrica
                    a_t += k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
                    L_t += ufl.inner(dolfinx.fem.Constant(self.mesh, q_vol), v_t) * dx
                
                # Applica Neumann BC
                if 'neumann_thermal' in self.bcs_info:
                    for bc_info in self.bcs_info['neumann_thermal']:
                        ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        # E anche qui, per coerenza
                        L_t += ufl.inner(dolfinx.fem.Constant(self.mesh, bc_info['value']), v_t) * ds
                
                # Definisci e risolvi il problema lineare
                petsc_opts = {"ksp_type": "cg", "pc_type": "gamg", "ksp_rtol": tol_th}
                problem_t = LinearProblem(a_t, L_t, bcs=bcs_t, u=T_sol, petsc_options=petsc_opts)
                problem_t.solve()
            
            # Applica rilassamento
            T_sol.x.array[:] = self.relax_T * T_sol.x.array + (1 - self.relax_T) * T_old.x.array

            # =========================================================
            #  2. RISOLUZIONE DEL PROBLEMA MECCANICO (globale)
            # =========================================================
            if self.on["mechanical"]:
                print("[INFO] Assembling and solving GLOBAL mechanical problem...")
                # --- INIZIO MODIFICA ---
                u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
                # Inizializza a 0
                a_m = 0
                L_m = 0
                # --- FINE MODIFICA ---
                # Assembla le forme variazionali su tutti i materiali
                for label, material in self.materials.items():
                    tag = self.label_map[label]
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)
                    lmbda = dolfinx.fem.Constant(self.mesh, material["lmbda"])
                    mu = dolfinx.fem.Constant(self.mesh, material["mu"])
                    alpha = dolfinx.fem.Constant(self.mesh, material["alpha"])
                    T_ref = dolfinx.fem.Constant(self.mesh, material["T_ref"])
                    rho = dolfinx.default_scalar_type(material["rho"])
                    g = dolfinx.default_scalar_type(self.g)
                    body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))
                    
                    epsilon = lambda u: ufl.sym(ufl.grad(u))
                    sigma_m = lambda u: lmbda * ufl.tr(epsilon(u)) * ufl.Identity(self.tdim) + 2 * mu * epsilon(u)
                    
                    delta_T = T_sol - T_ref # Usa la T appena calcolata
                    sigma_th = (3 * lmbda + 2 * mu) * alpha * delta_T * ufl.Identity(self.gdim)
                    
                    a_m += ufl.inner(sigma_m(u_m), epsilon(v_m)) * dx
                    L_m += ufl.dot(body_force, v_m) * dx + ufl.inner(sigma_th, epsilon(v_m)) * dx

                # Applica Traction BC
                if 'traction' in self.bcs_info:
                    for bc_info in self.bcs_info['traction']:
                        ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_m += ufl.dot(bc_info['value'], v_m) * ds # Assumendo che 'value' sia già un vettore

                # Definisci e risolvi il problema lineare
                petsc_opts = {"ksp_type": "cg", "pc_type": "gamg", "ksp_rtol": tol_mech, "pc_gamg_coarse_eq_limit": 1000}
                problem_m = LinearProblem(a_m, L_m, bcs=bcs_m, u=u_sol, petsc_options=petsc_opts)
                problem_m.solve()

                # Applica rilassamento
                u_sol.x.array[:] = self.relax_u * u_sol.x.array + (1 - self.relax_u) * u_old.x.array
            
            # =========================================================
            #  3. CONTROLLO CONVERGENZA e LOGGING
            # =========================================================
            

            # Sincronizza i dati dai processi ghost prima delle operazioni vettoriali
            T_sol.x.scatter_forward()
            T_old.x.scatter_forward()
            u_sol.x.scatter_forward()
            u_old.x.scatter_forward()

            # Accedi ai vettori PETSc sottostanti
            vec_T_now = T_sol.x.petsc_vec
            vec_T_old = T_old.x.petsc_vec
            vec_u_now = u_sol.x.petsc_vec
            vec_u_old = u_old.x.petsc_vec

            # Calcola la differenza: diff = vec_now - vec_old
            # 1. Crea una copia del vettore attuale
            diff_T = vec_T_now.copy()
            diff_u = vec_u_now.copy()
            
            # 2. Sottrai il vecchio vettore (diff = 1*diff - 1*vec_old)
            diff_T.axpy(-1.0, vec_T_old)
            diff_u.axpy(-1.0, vec_u_old)

            # Calcola la norma L2 del vettore differenza.
            # Questo è gestito correttamente in parallelo da PETSc.
            norm_dT = diff_T.norm(PETSc.NormType.NORM_2)
            norm_du = diff_u.norm(PETSc.NormType.NORM_2)
            
            # --- FINE BLOCCO DEFINITIVO ---
            
            print(f"\n[Convergence] Global ||ΔT|| = {norm_dT:.3e} (target: {tol_T:.1e})")
            print(f"[Convergence] Global ||Δu|| = {norm_du:.3e} (target: {tol_u:.1e})")


            if norm_dT < tol_T and norm_du < tol_u:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")
                # Aggiorna tutti i campi nei dizionari self.T e self.u
                for name in self.materials.keys():
                    self.T[name].interpolate(T_sol)
                    self.u[name].interpolate(u_sol)
                return True # Ritorna successo

            # Il tuo problema a senso unico dovrebbe convergere qui
            if iteration == 0:
                print("[INFO] One-way coupling: convergence should be achieved.")
                # Puoi decidere se forzare un break qui o lasciare che il ciclo verifichi la tolleranza.

        # Se il ciclo finisce senza convergenza
        print("\n[WARNING] Staggered solver did not converge within the maximum number of iterations.")
        # Aggiorna comunque i campi con l'ultima soluzione calcolata
        for name in self.materials.keys():
            self.T[name].interpolate(T_sol)
            self.u[name].interpolate(u_sol)
        return False # Ritorna fallimento

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
