import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from core.diagnostic import *

class Solver:
    def __init__(self):
        print("__Solver initializer__")
        solver_settings = self.config.get("solver_settings", {})

        self.relax_T = float(solver_settings.get("relax_T", 0.2))
        self.relax_u = float(solver_settings.get("relax_u", 0.5))

        self.thermal_solver = solver_settings.get("thermal_solver", "default")
        self.mechanical_solver = solver_settings.get("mechanical_solver", "default")

        print("  Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")

    def get_solver_options(self, physics, solver_type="iterative_amg", rtol=1e-10):
        """
        Returns PETSc options for the linear solver based on the physics.

        Parameters:
            physics (str): "thermal" or "mechanical". This primarily determines
                        if a symmetric (CG) or non-symmetric (GMRES) solver is used.
            solver_type (str): Type of solver/preconditioner strategy. Supported:
                - "direct_mumps": Robust but memory-intensive direct solver.
                - "iterative_amg": Default iterative solver with Algebraic Multigrid (PETSc GAMG).
                - "iterative_hypre": Iterative solver with HYPRE's BoomerAMG preconditioner.
            rtol (float): Relative tolerance for KSP.

        Returns:
            dict: PETSc options dictionary.
        """
        if physics not in ["thermal", "mechanical"]:
            raise ValueError(f"Unknown physics '{physics}'. Must be 'thermal' or 'mechanical'.")

        # 1. Choose KSP type based on matrix symmetry
        if physics == "thermal":
            ksp_type = "cg"  # Symmetric Positive Definite matrix
        elif physics == "mechanical":
            ksp_type = "gmres" # Non-symmetric matrix
        
        # 2. Choose preconditioner based on solver_type
        if solver_type == "direct_mumps":
            return {
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps" # MUMPS, direct solver
            }

        elif solver_type == "iterative_amg":
            coarse_eq_limit = 500 if physics == "thermal" else 1000
            return {
                "ksp_type": ksp_type,
                "pc_type": "gamg", # PETSc's built-in Algebraic Multigrid (GAMG)
                "ksp_rtol": rtol,
                "pc_gamg_coarse_eq_limit": coarse_eq_limit,
                # "ksp_gmres_restart": 30 # Restart for GMRES can be useful for memory
            }

        elif solver_type == "iterative_hypre":
            return {
                "ksp_type": ksp_type,
                "pc_type": "hypre",
                "pc_hypre_type": "boomeramg", # HYPRE's BoomerAMG, often more robust than GAMG
                "ksp_rtol": rtol
            }
        
        else:
            raise ValueError(f"Unknown solver_type '{solver_type}'.")
        
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
            dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.volume_tags, subdomain_id=tag)

            # Material properties
            k = dolfinx.fem.Constant(self.mesh, material["k"])
            lmbda = dolfinx.fem.Constant(self.mesh, material["lmbda"])
            G = dolfinx.fem.Constant(self.mesh, material["G"])
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
            sigma_m = lmbda * ufl.tr(self.epsilon(u_m)) * ufl.Identity(self.tdim) + 2 * G * self.epsilon(u_m)
            sigma_th = (3 * lmbda + 2 * G) * alpha * delta_T * ufl.Identity(self.gdim)
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
        # print("ksp_type, preonly")
        # opts[f"{option_prefix}ksp_type"] = "preonly"
        # opts[f"{option_prefix}pc_type"] = "lu"
        # opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
        # opts[f"{option_prefix}ksp_monitor"] = None

        # GAMG
        print("pc_gamg_type, agg")
        opts[f"{option_prefix}pc_gamg_type"] = "agg"
        opts[f"{option_prefix}pc_gamg_coarse_eq_limit"] = 1000
        opts[f"{option_prefix}mg_levels_ksp_type"] = "chebyshev"
        opts[f"{option_prefix}mg_levels_pc_type"] = "sor"
        opts[f"{option_prefix}mg_levels_esteig_ksp_type"] = "cg"
        opts[f"{option_prefix}mg_levels_ksp_chebyshev_esteig_steps"] = 20
        opts[f"{option_prefix}ksp_monitor"] = None

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
        print(f"  → Max iterations        : {max_iter}")
        print(f"  → Staggering tolerance |ΔT| : {tol_T:.1e}")
        print(f"  → Staggering tolerance |Δu| : {tol_u:.1e}")
        print(f"  → Thermal KSP tolerance   : {tol_th:.1e}")
        print(f"  → Mechanical KSP tolerance: {tol_mech:.1e}")

        # --- Initial setup ---
        # Use temporary functions for the current iteration's solution
        T_n = dolfinx.fem.Function(self.V_t)
        T_n.x.array[:] = self.T.x.array # Start with the initial/previous converged state

        u_n = dolfinx.fem.Function(self.V_m)
        u_n.x.array[:] = self.u.x.array

        # Functions to store the previous state for convergence check
        T_old = dolfinx.fem.Function(self.V_t)
        u_old = dolfinx.fem.Function(self.V_m)
        
        # Assemble all Dirichlet BCs
        bcs_t = []
        for bc_list in self.dirichlet_thermal.values():
            bcs_t.extend(bc_list)

        bcs_m = []
        for bc_list in self.dirichlet_mechanical.values():
            bcs_m.extend(bc_list)

        # --- Staggering loop ---
        for iteration in range(max_iter):
            print(f"\n--- Staggering iteration {iteration+1}/{max_iter} ---")

            # Store previous state for convergence check
            T_old.x.array[:] = T_n.x.array
            u_old.x.array[:] = u_n.x.array

            # =========================================================
            #  1. SOLVE THERMAL SUB-PROBLEM
            # =========================================================

            if self.on["thermal"]:
                print("\n[INFO] Assembling thermal problem...")
                u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
                a_t = 0
                L_t = 0

                # --- Volume integrals ---
                for label, material in self.materials.items():
                    tag = self.label_map[label] # in this case, volume tag
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.volume_tags, subdomain_id=tag) # volume integral

                    print(f"  Building weak form, volume integrals (dx) for {label}, tag = {tag}")

                    k = dolfinx.fem.Constant(self.mesh, material["k"])
                    a_t += k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
                    L_t += self.q_third * v_t * dx

                    # q_third value in the volume "label" / "tag"
                    dofs = self.locateDomainDofs(self.label_map, label, self.V_t)
                    q_vals = self.q_third.x.array[dofs]
                    print(f"  → q_third[{label}] min = {q_vals.min():.2e}, max = {q_vals.max():.2e}, mean = {q_vals.mean():.2e}")

                # --- Neumann BCs ---
                for label in self.materials:
                    for bc_info in self.neumann_thermal[label]:
                        print(f"  Applying thermal Neumann BC on subdomain id = {bc_info['id']}")
                        ds_neumann = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_t += bc_info['value'] * v_t * ds_neumann

                for label in self.materials:
                    for bc_info in self.robin_thermal[label]:
                        region_id = bc_info["id"]
                        pair_region = bc_info["pair"]

                        dofs_here = self.locateFacetDofs(region_id, self.V_t)
                        dofs_other = self.locateFacetDofs(self.label_map[pair_region], self.V_t)

                        T_here = T_n.x.array[dofs_here]
                        T_other = T_n.x.array[dofs_other]
                        self.gap_temperature = 0.5 * (T_here.mean() + T_other.mean())

                        # print(f"Average gap temperature between {label} and {pair_region}: {self.gap_temperature:.2f} K")

                        break
                    break

                h_gap = self.set_gap_conductance()     

                # --- Robin BCs ---
                for label in self.materials:
                    for bc_info in self.robin_thermal[label]:
                        region_id = bc_info["id"] # e.g., fuel_interface
                        pair_region = bc_info["pair"] # e.g., clad_interface

                        ds_interface = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=region_id) # e.g., fuel_interface surface

                        # 1) project
                        T_other = dolfinx.fem.Function(self.V_t) # e.g., T_clad
                        dofs_here = self.locateFacetDofs(region_id, self.V_t) # e.g., dofs of fuel_interface
                        dofs_other = self.locateFacetDofs(self.label_map[pair_region], self.V_t) # e.g., dofs of clad_interface

                        T_other.x.array[dofs_here] = T_n.x.array[dofs_other] # e.g., T_clad projected on fuel_interface dofs

                        # 2) mean
                        # T_other = T_n.x.array[dofs_other].mean()

                        # -) bias
                        # if label == 'fuel':
                        #     print(f"Biasing cladding temperature")
                        #     T_other.x.array[dofs_here[:200]] += 500

                        a_t += h_gap * u_t * v_t * ds_interface
                        L_t += h_gap * T_other * v_t * ds_interface

                        print(f"  [INFO] Gap Robin BC between '{label}' (region={region_id}) and '{pair_region}' (region={self.label_map[pair_region]})")

                # --- Solve the system ---
                petsc_opts_thermal = self.get_solver_options(
                    solver_type=self.thermal_solver, physics="thermal", rtol=tol_th
                )
                problem_t = LinearProblem(a_t, L_t, bcs=bcs_t, u=T_n, petsc_options=petsc_opts_thermal)
                problem_t.solve()

                # Apply relaxation
                T_n.x.array[:] = self.relax_T * T_n.x.array + (1 - self.relax_T) * T_old.x.array

                # Print heat flux
                # self.heat_flux(T_n)
            
            # =========================================================
            #  2. SOLVE MECHANICAL SUB-PROBLEM
            # =========================================================

            if self.on["mechanical"]:
                print("\n[INFO] Assembling mechanical problem...")
                u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
                a_m = 0
                L_m = 0
                
                # --- Volume integrals ---
                for label, material in self.materials.items():
                    tag = self.label_map[label]
                    print(f"  Building weak form, volume integrals (dx) for {label}, tag = {tag}")
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.volume_tags, subdomain_id=tag)

                    # Body force
                    rho = dolfinx.default_scalar_type(material["rho"])
                    g = dolfinx.default_scalar_type(self.g)
                    body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))

                    # Add contributions to global forms
                    a_m += ufl.inner(self.sigma_mech(u_m, material), self.epsilon(v_m)) * dx
                    L_m += ufl.dot(body_force, v_m) * dx + ufl.inner(self.sigma_th(T_n, material), self.epsilon(v_m)) * dx # dot (vectors, scalar product) inner (tensors, double contraction)
                    
                # --- Neumann BCs (Traction) ---
                for label in self.materials:
                    for bc_info in self.traction[label]:
                        print(f"  Applying mechanical traction on subdomain id = {bc_info['id']}")
                        ds_traction = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_m += ufl.dot(bc_info['value'], v_m) * ds_traction

                # --- Solve the system ---
                petsc_opts_mech = self.get_solver_options(
                    solver_type=self.mechanical_solver, physics="mechanical", rtol=tol_mech
                )
                problem_m = LinearProblem(a_m, L_m, bcs=bcs_m, u=u_n, petsc_options=petsc_opts_mech)
                problem_m.solve()

                # Apply relaxation
                u_n.x.array[:] = self.relax_u * u_n.x.array + (1 - self.relax_u) * u_old.x.array

            # =========================================================
            #  3. CONVERGENCE CHECK and LOGGING
            # =========================================================

            print(f"\nConvergence check")
            
            # Synchronize ghost values before vector operations
            T_n.x.scatter_forward()
            T_old.x.scatter_forward()
            u_n.x.scatter_forward()
            u_old.x.scatter_forward()

            # --- Temperature convergence ---
            # Get the PETSc vectors from the function's data layout '.x'
            vec_T_new = T_n.x.petsc_vec
            vec_T_old = T_old.x.petsc_vec

            # Compute difference vector: diff_T = T_now - T_old
            diff_T = vec_T_new.copy()
            diff_T.axpy(-1.0, vec_T_old) 

            # Compute norms
            norm_dT = diff_T.norm(PETSc.NormType.NORM_2)
            norm_T = vec_T_new.norm(PETSc.NormType.NORM_2)
            rel_norm_dT = norm_dT / norm_T if norm_T > 1e-12 else norm_dT

            # --- Displacement convergence ---
            # Get the PETSc vectors
            vec_u_new = u_n.x.petsc_vec
            vec_u_old = u_old.x.petsc_vec

            # Compute difference vector: diff_u = u_now - u_old
            diff_u = vec_u_new.copy()
            diff_u.axpy(-1.0, vec_u_old)

            # Compute norms
            norm_du = diff_u.norm(PETSc.NormType.NORM_2)
            norm_u = vec_u_new.norm(PETSc.NormType.NORM_2)
            rel_norm_du = norm_du / norm_u if norm_u > 1e-12 else norm_du

            print(f"  ||ΔT|| = {norm_dT:.3e}, ||Δu|| = {norm_du:.3e}")
            print(f"  ||ΔT||/||T|| = {rel_norm_dT:.3e}, ||Δu||/||u|| = {rel_norm_du:.3e}")

            if rel_norm_dT < tol_T and rel_norm_du < tol_u:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")

                # --- Update the main class state ONLY on successful convergence ---
                self.T.x.array[:] = T_n.x.array
                self.u.x.array[:] = u_n.x.array
                
                T_min, T_max = self.T.x.array.min(), self.T.x.array.max()
                print(f"Global min/max temperature: {T_min:.2f} / {T_max:.2f} K")
                
                # Get the raw data from the displacement solution vector
                u_array = u_n.x.array
                u_vectors = u_array.reshape(-1, self.gdim)
                magnitude_values = np.sqrt(np.sum(u_vectors**2, axis=1))
                u_min = magnitude_values.min()
                u_max = magnitude_values.max()
                print(f"Global min/max displacement magnitude: {u_min:.2e} / {u_max:.2e} m")

                return True

        print("\n[WARNING] Staggered solver did not converge. The state from the last iteration is used.")

        # Update the main state with the last (non-converged) result
        self.T.x.array[:] = T_n.x.array
        self.u.x.array[:] = u_n.x.array
        return False