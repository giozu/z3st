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
        """
        Monolithic thermo-mechanical solve (MVP, Phase 1):

        Supported features:
        - Dirichlet/Neumann boundary conditions for both thermal and mechanical problems
          (including Clamp_x/y/z constraints).
        - Material thermal conductivity `k` can be constant or UFL symbolic expression
          (already stored in self.materials[*]["k"]).

        Not included in this phase:
        - Robin-type gap conduction
        - Interface projection operators
        - "Gas" gap conductance h_gap

        NOTE:
        - Slip_* boundary conditions and Robin conditions are not yet implemented here.
        - Dirichlet BCs are reconstructed on the mixed space W (W.sub(1) = temperature, W.sub(0) = displacement).
        """

        print("\n[INFO] Starting monolithic thermo-mechanical solve (Phase 1)")
        print(f"  → Tolerance (Newton)     : {tol}")
        print(f"  → Max iterations (Newton): {max_iter}")

        # Split (initialized) mixed unknowns into mechanical (u_m) and thermal (u_t) parts
        (u_m, u_t) = ufl.split(self.sol_mixed)
        (v_m, v_t) = ufl.TestFunctions(self.W)

        # Initialize total residual
        F = 0

        # ===== VOLUME CONTRIBUTIONS: loop over materials =====
        for label, material in self.materials.items():
            tag = self.label_map[label]
            dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.volume_tags, subdomain_id=tag)

            # Material properties (may be constants or UFL expressions consistent with self.T)
            k = material["k"]

            rho = dolfinx.default_scalar_type(material["rho"])
            g = dolfinx.default_scalar_type(self.g)
            body_force = dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -rho * g))

            # Volumetric heat source: global Function self.q_third
            q_vol = self.q_third

            # Thermal contribution
            F_thermal = (k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) - q_vol * v_t) * dx_local
            if self.on["thermal"]:
                F += F_thermal

            # Mechanical stress contributions (mechanical + thermal)
            sigma_m = self.sigma_mech(u_m, material)
            if self.on["thermal"]:
                sigma_th = self.sigma_th(u_t, material)
                sigma_tot = sigma_m - sigma_th
            else:
                sigma_tot = sigma_m

            if self.on["mechanical"]:
                F += (ufl.inner(sigma_tot, ufl.sym(ufl.grad(v_m))) - ufl.dot(body_force, v_m)) * dx_local

        # ===== BOUNDARY CONTRIBUTIONS =====
        n_vec = ufl.FacetNormal(self.mesh)

        # Thermal Neumann BCs (heat flux) — reuse lists from set_thermal_boundary_conditions
        for label in self.materials:
            for bc_info in self.neumann_thermal.get(label, []):
                ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info["id"])
                F -= bc_info["value"] * v_t * ds_local  # outward heat flux term

        # Mechanical traction BCs — reuse lists from set_mechanical_boundary_conditions
        for label in self.materials:
            for bc_info in self.traction.get(label, []):
                ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info["id"])
                # bc_info["value"] is already a traction vector (Constant * n), use directly
                F -= ufl.dot(bc_info["value"], v_m) * ds_local

        # ===== DIRICHLET CONDITIONS on mixed space W.sub(1) (T) and W.sub(0) (u) =====
        bcs_mixed = []

        # Helper: locate DOFs on a subspace for a given facet region
        def locate_dofs_on_sub(subspace, region_id):
            facets = self.ft.find(region_id)
            return dolfinx.fem.locate_dofs_topological(subspace, self.fdim, facets)

        # Thermal Dirichlet BCs reconstructed on W.sub(1)
        thermal_bcs_defs = self.boundary_conditions.get("thermal_bcs", {})
        for mat, bc_list in thermal_bcs_defs.items():
            for bc in bc_list:
                if bc.get("type") != "Dirichlet":
                    continue
                region = bc.get("region")
                Tval = bc.get("temperature", None)
                if region is None or Tval is None:
                    continue
                region_id = self.label_map.get(region)
                if region_id is None:
                    continue
                dofs_T = locate_dofs_on_sub(self.W.sub(1), region_id)
                T_d = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(Tval))
                bcs_mixed.append(dolfinx.fem.dirichletbc(T_d, dofs_T, self.W.sub(1)))

        # Mechanical Dirichlet BCs: full vector displacement or single-component clamps
        mech_bcs_defs = self.boundary_conditions.get("mechanical_bcs", {})
        for mat, bc_list in mech_bcs_defs.items():
            for bc in bc_list:
                bc_type = bc.get("type")
                region = bc.get("region")
                if region is None or bc_type is None:
                    continue
                region_id = self.label_map.get(region)
                if region_id is None:
                    continue

                if bc_type == "Dirichlet":
                    disp = bc.get("displacement", None)
                    if disp is None:
                        continue
                    # Collapse vector subspace (u) to a full FunctionSpace
                    V_u_sub, _ = self.W.sub(0).collapse()
                    # Create constant displacement Function
                    vec = np.asarray(disp, dtype=dolfinx.default_scalar_type).reshape(self.gdim)
                    u_d_fun = dolfinx.fem.Function(V_u_sub, name="u_dirichlet_const")
                    u_d_fun.interpolate(lambda x, v=vec: np.tile(v.reshape(self.gdim, 1), (1, x.shape[1])))
                    # Locate DOFs mapping mixed space → collapsed space
                    facets = self.ft.find(region_id)
                    dofs_mixed = dolfinx.fem.locate_dofs_topological((self.W.sub(0), V_u_sub), self.fdim, facets)
                    bcs_mixed.append(dolfinx.fem.dirichletbc(u_d_fun, dofs_mixed, self.W.sub(0)))

                elif bc_type == "Clamp_x":
                    dofs_x = locate_dofs_on_sub(self.W.sub(0).sub(0), region_id)
                    bcs_mixed.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0.0), dofs_x, self.W.sub(0).sub(0)))

                elif bc_type == "Clamp_y":
                    dofs_y = locate_dofs_on_sub(self.W.sub(0).sub(1), region_id)
                    bcs_mixed.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0.0), dofs_y, self.W.sub(0).sub(1)))

                elif bc_type == "Clamp_z":
                    val = float(bc.get("value", 0.0))
                    dofs_z = locate_dofs_on_sub(self.W.sub(0).sub(2), region_id)
                    bcs_mixed.append(dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(val), dofs_z, self.W.sub(0).sub(2)))

        # ===== NONLINEAR PROBLEM & NEWTON SOLVER =====
        problem = NonlinearProblem(F, self.sol_mixed, bcs=bcs_mixed)
        solver = NewtonSolver(self.mesh.comm, problem)

        # PETSc KSP options for GAMG, similar to staggered solver settings
        ksp = solver.krylov_solver
        opts = PETSc.Options()
        pref = ksp.getOptionsPrefix()
        opts[f"{pref}pc_gamg_type"] = "agg"
        opts[f"{pref}pc_gamg_coarse_eq_limit"] = 1000
        opts[f"{pref}mg_levels_ksp_type"] = "chebyshev"
        opts[f"{pref}mg_levels_pc_type"] = "sor"
        opts[f"{pref}mg_levels_esteig_ksp_type"] = "cg"
        opts[f"{pref}mg_levels_ksp_chebyshev_esteig_steps"] = 20
        opts[f"{pref}ksp_monitor"] = None
        opts["options_left"] = "0"
        ksp.setFromOptions()

        solver.atol = tol
        solver.rtol = tol
        solver.convergence_criterion = "incremental"
        solver.max_it = max_iter

        print("[INFO] Solving monolithic nonlinear problem...")
        iters, converged = solver.solve(self.sol_mixed)

        if converged:
            print(f"[INFO] Newton solver converged in {iters} iterations")
        else:
            print(f"[WARNING] Newton solver did not converge after {iters} iterations")

        # ===== UPDATE GLOBAL STATE FIELDS =====
        u_sol = self.sol_mixed.sub(0).collapse()
        u_sol.name = "Displacement"
        T_sol = self.sol_mixed.sub(1).collapse()
        T_sol.name = "Temperature"

        # Copy results back into global Functions self.u / self.T
        try:
            self.u.interpolate(u_sol)
        except Exception:
            self.u.x.array[:] = u_sol.x.array
        try:
            self.T.interpolate(T_sol)
        except Exception:
            self.T.x.array[:] = T_sol.x.array

        print(f"Min/Max temperature: {self.T.x.array.min():.2f} / {self.T.x.array.max():.2f} K")
        u_vec = self.u.x.array.reshape(-1, self.gdim)
        umag = np.linalg.norm(u_vec, axis=1)
        print(f"Min/Max displacement magnitude: {umag.min():.2e} / {umag.max():.2e}")

    def solve_staggered(self, max_iter=20, tol_u=1e-2, tol_T=1e-2, tol_mech=1e-6, tol_th=1e-6, tol_d=1e-3):
        print(f"  → Max iterations              : {max_iter}")
        print(f"  → Staggering tolerance |ΔT|   : {tol_T:.1e}")
        print(f"  → Staggering tolerance |Δu|   : {tol_u:.1e}")
        print(f"  → Staggering tolerance |Δd|   : {tol_d:.1e}")
        print(f"  → Thermal KSP tolerance       : {tol_th:.1e}")
        print(f"  → Mechanical KSP tolerance    : {tol_mech:.1e}")

        # --- Initial setup ---
        # Use temporary functions for the current iteration's solution
        T_n = dolfinx.fem.Function(self.V_t)
        T_n.x.array[:] = self.T.x.array # Start with the initial/previous converged state

        u_n = dolfinx.fem.Function(self.V_m)
        u_n.x.array[:] = self.u.x.array

        # Functions to store the previous state for convergence check
        T_old = dolfinx.fem.Function(self.V_t)
        u_old = dolfinx.fem.Function(self.V_m)

        # d (damage)
        if self.on.get("fracture", False):
            d_old = dolfinx.fem.Function(self.V_d)
            d_old.x.array[:] = self.d.x.array
        else:
            d_old = None

        # Assemble all Dirichlet BCs
        bcs_t = [];  [bcs_t.extend(bc_list) for bc_list in self.dirichlet_thermal.values()]
        bcs_m = [];  [bcs_m.extend(bc_list) for bc_list in self.dirichlet_mechanical.values()]

        # --- Staggering loop ---
        for iteration in range(max_iter):
            print(f"\n--- Staggering iteration {iteration+1}/{max_iter} ---")

            # Store previous state for convergence check
            T_old.x.array[:] = T_n.x.array
            u_old.x.array[:] = u_n.x.array
            if self.on.get("fracture", False):
                d_old.x.array[:] = self.d.x.array
                self.d.x.array[:] = np.maximum(self.d.x.array, self.d_n.x.array)

            # 1) THERMAL
            if self.on.get("thermal", False):
                print("\n[INFO] Assembling thermal problem...")
                u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
                a_t = 0
                L_t = 0

                # --- Volume integrals ---
                for label, material in self.materials.items():
                    tag = self.label_map[label] # in this case, volume tag
                    dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.volume_tags, subdomain_id=tag) # volume integral

                    print(f"\n  Building weak form, volume integrals (dx) for {label}, tag = {tag}")

                    # Thermal conductivity
                    k = material["k"]

                    # Output the value of the average thermal conductivity
                    if isinstance(k, (int, float)):
                        k_avg = k
                    else:
                        V0 = dolfinx.fem.functionspace(self.mesh, ("CG", 1))
                        k_func = dolfinx.fem.Function(V0, name=f"k_{label}")
                        expr = dolfinx.fem.Expression(k, V0.element.interpolation_points())
                        k_func.interpolate(expr)

                        num = dolfinx.fem.assemble_scalar(dolfinx.fem.form(k_func * dx))
                        den = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx))
                        k_avg = num / den

                    print(f"  Average conductivity for {label} = {k_avg:.3f} W/m·K")

                    a_t += k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
                    L_t += self.q_third * v_t * dx

                    # q_third value in the volume "label" / "tag"
                    dofs = self.locateDomainDofs(self.label_map, label, self.V_t)
                    q_vals = self.q_third.x.array[dofs]
                    print(f"  → q_third[{label}] min = {q_vals.min():.2e}, max = {q_vals.max():.2e}, mean = {q_vals.mean():.2e}")

                # Neumann
                for label in self.materials:
                    for bc_info in self.neumann_thermal[label]:
                        print(f"  Applying thermal Neumann BC on subdomain id = {bc_info['id']}")
                        ds_neumann = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_t += bc_info['value'] * v_t * ds_neumann

                # Gap (Robin)
                h_gap = self.set_gap_conductance(T_n)

                for label in self.materials:
                    for bc_info in self.robin_thermal[label]:
                        print(f"  Applying thermal Robin BC on subdomain id = {bc_info['id']}")

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

                # Check T_fo
                # dofs = self.locateFacetDofs(self.label_map["lateral"], self.V_t)
                # T = T_n.x.array[dofs].mean()
                # print(f"  T_lateral = {T:.2f} K")

                # Apply relaxation
                T_n.x.array[:] = self.relax_T * T_n.x.array + (1 - self.relax_T) * T_old.x.array

                # Print heat flux
                # self.heat_flux(T_n)

            # 2) PHASE FIELD (damage)
            if self.on.get("fracture", False):
                print("\n[INFO] Assembling phase-field (damage) problem...")
                self.update_history(u_n)      # H(u_n) irreversibilità
                self.build_damage_forms()     # forme AT-2 con H aggiornato
                self.assemble_solve_damage()  # solve per d (con clamp [0,1])
                np.maximum(self.d.x.array, self.d_n.x.array, out=self.d.x.array)
                self.d_n.x.array[:] = self.d.x.array

            # 3) MECHANICAL SUB-PROBLEM
            if self.on.get("mechanical", False):
                print("\n[INFO] Assembling mechanical problem...")

                if self.on.get("fracture", False):
                    self.update_split_indicator(u_n)

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

                    # Mechanical stress, linearized, degradated
                    if self.on.get("fracture", False):
                        sigma_el = self.sigma_damage_linear(u_m, self.d, material)
                    else:
                        sigma_el = self.sigma_mech(u_m, material)

                    a_m += ufl.inner(sigma_el, self.epsilon(v_m)) * dx
                    L_m += ufl.dot(body_force, v_m) * dx
                    L_m += ufl.inner(self.sigma_th(T_n, material), self.epsilon(v_m)) * dx
                    
                # --- Neumann BCs (Traction) ---
                for label in self.materials:
                    for bc_info in self.traction[label]:
                        print(f"  Applying mechanical traction on subdomain id = {bc_info['id']}")
                        ds_traction = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=bc_info['id'])
                        L_m += ufl.dot(bc_info['value'], v_m) * ds_traction

                if self.on.get("fracture", False):
                    self.update_split_indicator(u_n)

                # --- Solve the system ---
                petsc_opts_mech = self.get_solver_options(
                    solver_type=self.mechanical_solver, physics="mechanical", rtol=tol_mech
                )
                problem_m = LinearProblem(a_m, L_m, bcs=bcs_m, u=u_n, petsc_options=petsc_opts_mech)
                problem_m.solve()

                # Apply relaxation
                u_n.x.array[:] = self.relax_u * u_n.x.array + (1 - self.relax_u) * u_old.x.array

            # 4) CONVERGENCE CHECK and LOGGING
            print(f"\nConvergence check")
            
            # Synchronize ghost values before vector operations
            T_n.x.scatter_forward()
            T_old.x.scatter_forward()
            u_n.x.scatter_forward()
            u_old.x.scatter_forward()
            if self.on.get("fracture", False):
                self.d.x.scatter_forward(); d_old.x.scatter_forward()

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

            conv_d = True
            d_diff_inf = 0.0
            if self.on.get("fracture", False):
                d_diff_inf = np.linalg.norm(self.d.x.array - d_old.x.array, ord=np.inf)
                conv_d = d_diff_inf < tol_d

            print(f"  ||ΔT|| = {norm_dT:.3e}, ||Δu|| = {norm_du:.3e}, ||Δd||_inf = {d_diff_inf:.3e}")
            print(f"  ||ΔT||/||T|| = {rel_norm_dT:.3e}, ||Δu||/||u|| = {rel_norm_du:.3e}  (tol_T={tol_T:.1e}, tol_u={tol_u:.1e}, tol_d={tol_d:.1e})")
            if self.on.get("fracture", False):
                print(f"Damage d: min={self.d.x.array.min():.3e}, max={self.d.x.array.max():.3e}")
                print(f"History H: min={self.H.x.array.min():.3e}, max={self.H.x.array.max():.3e}")

            if (rel_norm_dT < tol_T or not self.on.get("thermal", False)) \
            and (rel_norm_du < tol_u or not self.on.get("mechanical", False)) \
            and conv_d:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")
                self.T.x.array[:] = T_n.x.array
                self.u.x.array[:] = u_n.x.array
                if self.on.get("fracture", False):
                    self.d.x.array[:] = self.d.x.array
                    self.d_n.x.array[:] = self.d.x.array

                T_min, T_max = self.T.x.array.min(), self.T.x.array.max()
                print(f"Global min/max temperature: {T_min:.2f} / {T_max:.2f} K")

                # # Get the raw data from the displacement solution vector
                # u_array = u_n.x.array
                # u_vectors = u_array.reshape(-1, self.gdim)
                # magnitude_values = np.sqrt(np.sum(u_vectors**2, axis=1))
                # u_min = magnitude_values.min()
                # u_max = magnitude_values.max()
                # print(f"Global min/max displacement magnitude: {u_min:.2e} / {u_max:.2e} m")

                u_vec = u_n.x.array.reshape(-1, self.gdim)
                umag = np.sqrt((u_vec**2).sum(axis=1))
                print(f"Global min/max displacement magnitude: {umag.min():.2e} / {umag.max():.2e} m")
                return True

        print("\n[WARNING] Staggered solver did not converge. Using last iteration state.")
        # Update the main state with the last (non-converged) result
        self.T.x.array[:] = T_n.x.array
        self.u.x.array[:] = u_n.x.array

        if self.on.get("fracture", False):
            self.d_n.x.array[:] = self.d.x.array

        return False
