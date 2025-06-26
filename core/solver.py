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
            dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)

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

        # --- Initial Setup ---
        T_sol = self.T  # Reference function for the temperature solution
        u_sol = self.u  # Reference function for the displacement solution

        # Functions to store the previous state
        T_old = {}
        for label, T_func in self.T.items():
            T_old[label] = dolfinx.fem.Function(T_func.function_space)
        u_old = dolfinx.fem.Function(u_sol.function_space)

        # Assemble all Dirichlet BCs from per-material dictionaries
        bcs_t = []
        for bc_list in self.dirichlet_thermal.values():
            bcs_t.extend(bc_list)

        bcs_m = []
        for bc_list in self.dirichlet_mechanical.values():
            bcs_m.extend(bc_list)
        # bcs_m.append(self.constrain_point(self.mesh, self.V_m, np.array([0.1, 0.1, 0.5])))

        print(f"[INFO] Loaded {len(bcs_t)} thermal and {len(bcs_m)} mechanical Dirichlet BCs for the staggered solver.")

        # --- Staggering Loop ---
        for iteration in range(max_iter):
            print(f"\n--- Staggering iteration {iteration+1}/{max_iter} ---")

            # Store previous state for convergence check and relaxation
            for label in self.T:
                T_old[label].x.array[:] = T_sol[label].x.array

            u_old.x.array[:] = u_sol.x.array

            # =========================================================
            #  1. SOLVE THERMAL SUB-PROBLEM
            # =========================================================
            if self.on["thermal"]:

                u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
                a_t = {}
                L_t = {}

                for label, material in self.materials.items():
                    print(f"[INFO] Assembling thermal problem for material: {label}")

                    if label not in self.label_map:
                        print(f"[WARNING] Material label '{label}' not found in label_map → skipping.")
                        continue

                    a_t[label], L_t[label] = 0, 0

                    tag = self.label_map[label]
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)

                    k = dolfinx.fem.Constant(self.mesh, material["k"])
                    q_vol = dolfinx.fem.Constant(self.mesh, self.q_third[label])

                    a_t[label] += k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
                    L_t[label] += q_vol * v_t * dx

                    # Apply Neumann BCs
                    for bc_info in self.neumann_thermal[label]:
                        ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_t[label] += bc_info['value'] * v_t * ds

                    # Thermal coupling fuel-cladding interface
                    h_gap_value = 10.0
                    h_gap = dolfinx.fem.Constant(self.mesh, h_gap_value)

                    # Fuel contribution
                    if "fuel" in self.materials and label == "fuel":
                        ds_fuel_interface = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=self.label_map["face_xx_1_1"])
                        a_t[label] += h_gap * ufl.inner(u_t, v_t) * ds_fuel_interface
                        L_t[label] += h_gap * self.T["cladding"] * v_t * ds_fuel_interface

                    # Cladding contribution
                    if "cladding" in self.materials and label == "cladding":
                        ds_clad_interface = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=self.label_map["face_xx_0_2"])
                        a_t[label] -= h_gap * ufl.inner(u_t, v_t) * ds_clad_interface
                        L_t[label] -= h_gap * self.T["fuel"] * v_t * ds_clad_interface

                    # Sparsity pattern
                    # assemble_and_analyze_matrix(a_t[label], L_t[label], title=f"Thermal stiffness matrix - {label}")

                for label, material in self.materials.items():

                    petsc_opts_thermal = self.get_solver_options(
                        solver_type=self.thermal_solver,
                        physics="thermal",
                        rtol=tol_th
                    )

                    problem_t = LinearProblem(a_t[label], L_t[label],
                        bcs=self.dirichlet_thermal[label],
                        u=self.T[label],
                                            petsc_options=petsc_opts_thermal)
                    problem_t.solve()

                # Apply relaxation
                for label in self.T:
                    T_sol[label].x.array[:] = self.relax_T * T_sol[label].x.array + (1 - self.relax_T) * T_old[label].x.array
                
            # =========================================================
            #  2. SOLVE MECHANICAL SUB-PROBLEM (global)
            # =========================================================
            if self.on["mechanical"]:
                print("[INFO] Assembling and solving GLOBAL mechanical problem...")
                u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
                a_m, L_m = 0, 0

                for label, material in self.materials.items():
                    tag = self.label_map[label]
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=tag)

                    lmbda = dolfinx.fem.Constant(self.mesh, material["lmbda"])
                    G = dolfinx.fem.Constant(self.mesh, material["G"])
                    alpha = dolfinx.fem.Constant(self.mesh, material["alpha"])
                    T_ref = dolfinx.fem.Constant(self.mesh, material["T_ref"])
                    rho = dolfinx.default_scalar_type(material["rho"])
                    g = dolfinx.default_scalar_type(self.g)
                    body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))

                    epsilon = lambda u: ufl.sym(ufl.grad(u))
                    sigma_m = lambda u: lmbda * ufl.tr(epsilon(u)) * ufl.Identity(self.tdim) + 2 * G * epsilon(u)

                    sigma_th = (3 * lmbda + 2 * G) * alpha * (T_sol[label] - T_ref) * ufl.Identity(self.gdim)

                    a_m += ufl.inner(sigma_m(u_m), epsilon(v_m)) * dx
                    L_m += ufl.dot(body_force, v_m) * dx + ufl.inner(sigma_th, epsilon(v_m)) * dx

                    for bc_info in self.traction[label]:
                        ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_m += ufl.dot(bc_info['value'], v_m) * ds

                # Sparsity pattern
                # assemble_and_analyze_matrix(a_m, L_m, title=f"Mechanical stiffness matrix")

                petsc_opts_mech = self.get_solver_options(
                    solver_type=self.mechanical_solver,
                    physics="mechanical",
                    rtol=tol_mech
                )

                problem_m = LinearProblem(a_m, L_m, bcs=bcs_m, u=u_sol, petsc_options=petsc_opts_mech)
                problem_m.solve()

                u_sol.x.array[:] = self.relax_u * u_sol.x.array + (1 - self.relax_u) * u_old.x.array
                
            # =========================================================
            #  3. CONVERGENCE CHECK and LOGGING
            # =========================================================
            # Synchronize ghost values before vector operations
            for label in self.T:
                T_sol[label].x.scatter_forward()
                T_old[label].x.scatter_forward()

            u_sol.x.scatter_forward()
            u_old.x.scatter_forward()

            norm_dT_total = 0.0

            for label in self.T:
                vec_T_now = self.T[label].x.petsc_vec
                vec_T_old = T_old[label].x.petsc_vec

                diff_T = vec_T_now.copy()
                diff_T.axpy(-1.0, vec_T_old)

                norm_dT_local = diff_T.norm(PETSc.NormType.NORM_2)
                norm_dT_total += norm_dT_local**2

            norm_dT = norm_dT_total**0.5

            vec_u_now, vec_u_old = u_sol.x.petsc_vec, u_old.x.petsc_vec
            diff_u = vec_u_now.copy(); diff_u.axpy(-1.0, vec_u_old)
            norm_du = diff_u.norm(PETSc.NormType.NORM_2)

            print(f"||ΔT|| = {norm_dT:.3e} | ||Δu|| = {norm_du:.3e}")

            if norm_dT < tol_T and norm_du < tol_u:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")

                for label in self.T:
                    T_min = self.T[label].x.array.min()
                    T_max = self.T[label].x.array.max()
                    print(f"Material {label}: Min/Max temperature: {T_min:.2f} / {T_max:.2f} K")

                u_magnitude = np.sqrt(np.sum(u_sol.x.array.reshape(-1, 3)**2, axis=1))
                print(f"Min/Max displacement magnitude: {u_magnitude.min():.3e} / {u_magnitude.max():.3e} m")
                print("----------------------------\n")

                # Update all fields in the self.T and self.u dictionaries
                self.u.interpolate(u_sol)

                return True

        print("\n[WARNING] Staggered solver did not converge within the maximum number of iterations.")
        self.u.interpolate(u_sol)

        return False
