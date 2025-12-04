# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import numpy as np
import ufl
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from petsc4py import PETSc


class Solver:
    def __init__(self):
        print("[Solver] initializer")
        solver_settings = self.input_file.get("solver_settings", {})
        self.coupling = self.input_file.get("solver_settings", {}).get("coupling", "staggered")

        self.relax_T = float(solver_settings.get("relax_T", 0.9))
        self.relax_u = float(solver_settings.get("relax_u", 0.4))

        print("  Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")

        self.relax_adaptive = bool(solver_settings.get("relax_adaptive", False))
        self.relax_growth = float(solver_settings.get("relax_growth", 1.2))
        self.relax_shrink = float(solver_settings.get("relax_shrink", 0.5))
        self.relax_min = float(solver_settings.get("relax_min", 0.05))
        self.relax_max = float(solver_settings.get("relax_max", 1.0))

        if self.relax_adaptive:
            print("  Adaptive relaxation enabled")
            print(f"  → relax_growth  : {self.relax_growth}")
            print(f"  → relax_shrink : {self.relax_shrink}")
            print(f"  → relax_min  : {self.relax_min}")
            print(f"  → relax_max : {self.relax_max}")

        else:
            print("  Adaptive relaxation disabled")

        print("\n")

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
            ksp_type = "gmres"  # Non-symmetric matrix

        # 2. Choose preconditioner based on solver_type
        if solver_type == "direct_mumps":
            return {
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",  # MUMPS, direct solver
            }

        elif solver_type == "iterative_amg":
            coarse_eq_limit = 500 if physics == "thermal" else 1000
            return {
                "ksp_type": ksp_type,
                "pc_type": "gamg",  # PETSc's built-in Algebraic Multigrid (GAMG)
                "ksp_rtol": rtol,
                "pc_gamg_coarse_eq_limit": coarse_eq_limit,
                # "ksp_gmres_restart": 30 # Restart for GMRES can be useful for memory
            }

        elif solver_type == "iterative_hypre":
            return {
                "ksp_type": ksp_type,
                "pc_type": "hypre",
                "pc_hypre_type": "boomeramg",  # HYPRE's BoomerAMG, often more robust than GAMG
                "ksp_rtol": rtol,
            }

        else:
            raise ValueError(f"Unknown solver_type '{solver_type}'.")

    def solve_monolithic(self, tol=1e-8, max_iter=50):
        """
        Monolithic thermo-mechanical solver (MVP, Phase 1).

        Supported features
        ------------------
        - Dirichlet and Neumann boundary conditions for both thermal and mechanical problems (including Clamp_x/y/z constraints).
        - Material thermal conductivity ``k`` can be constant or a UFL symbolic expression (already stored in ``self.materials[*]["k"]``).

        Not included in this phase
        --------------------------
        - Robin-type gap conduction
        - Interface projection operators
        - "Gas" gap conductance ``h_gap``

        Notes
        -----
        - Slip_* boundary conditions and Robin conditions are not yet implemented here.
        - Dirichlet BCs are reconstructed on the mixed space ``W``

        """

        print("\n[INFO] Starting monolithic thermo-mechanical solve (Phase 1)")
        print(f"  → Tolerance (Newton)     : {tol}")
        print(f"  → Max iterations (Newton): {max_iter}")

        # Split (initialized) mixed unknowns into mechanical (u_m) and thermal (u_t) parts
        (u_m, u_t) = ufl.split(self.sol_mixed)
        (v_m, v_t) = ufl.TestFunctions(self.W)

        # Initialize total residual
        F = 0

        # --. VOLUME CONTRIBUTIONS: loop over materials --..
        for label, material in self.materials.items():
            tag = self.label_map[label]
            dx_local = ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag
            )

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
                sigma_tot = sigma_m + sigma_th
            else:
                sigma_tot = sigma_m

            if self.on["mechanical"]:
                F += (
                    ufl.inner(sigma_tot, ufl.sym(ufl.grad(v_m))) - ufl.dot(body_force, v_m)
                ) * dx_local

        # Thermal Neumann BCs (heat flux) — reuse lists from set_thermal_boundary_conditions
        for label in self.materials:
            for bc_info in self.neumann_thermal.get(label, []):
                ds_local = ufl.Measure(
                    "ds",
                    domain=self.mesh,
                    subdomain_data=self.facet_tags,
                    subdomain_id=bc_info["id"],
                )
                F -= bc_info["value"] * v_t * ds_local  # outward heat flux term

        # Mechanical traction BCs — reuse lists from set_mechanical_boundary_conditions
        for label in self.materials:
            for bc_info in self.traction.get(label, []):
                ds_local = ufl.Measure(
                    "ds",
                    domain=self.mesh,
                    subdomain_data=self.facet_tags,
                    subdomain_id=bc_info["id"],
                )
                # bc_info["value"] is already a traction vector (Constant * n), use directly
                F -= ufl.dot(bc_info["value"], v_m) * ds_local

        # --. DIRICHLET CONDITIONS on mixed space W.sub(1) (T) and W.sub(0) (u) --..
        bcs_mixed = []

        # Helper: locate DOFs on a subspace for a given facet region
        def locate_dofs_on_sub(subspace, region_id):
            facets = self.facet_tags.find(region_id)
            return dolfinx.fem.locate_dofs_topological(subspace, self.fdim, facets)

        # Thermal Dirichlet BCs reconstructed on W.sub(1)
        thermal_bcs_defs = self.boundary_conditions.get("thermal_bcs", {})
        for _, bc_list in thermal_bcs_defs.items():
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
        for _, bc_list in mech_bcs_defs.items():
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
                    vec = np.asarray(disp, dtype=dolfinx.default_scalar_type).reshape(self.tdim)
                    u_d_fun = dolfinx.fem.Function(V_u_sub, name="u_dirichlet_const")
                    u_d_fun.interpolate(
                        lambda x, v=vec: np.tile(v.reshape(self.tdim, 1), (1, x.shape[1]))
                    )
                    # Locate DOFs mapping mixed space → collapsed space
                    facets = self.facet_tags.find(region_id)
                    dofs_mixed = dolfinx.fem.locate_dofs_topological(
                        (self.W.sub(0), V_u_sub), self.fdim, facets
                    )
                    bcs_mixed.append(dolfinx.fem.dirichletbc(u_d_fun, dofs_mixed, self.W.sub(0)))

                elif bc_type == "Clamp_x":
                    dofs_x = locate_dofs_on_sub(self.W.sub(0).sub(0), region_id)
                    bcs_mixed.append(
                        dolfinx.fem.dirichletbc(
                            dolfinx.default_scalar_type(0.0), dofs_x, self.W.sub(0).sub(0)
                        )
                    )

                elif bc_type == "Clamp_y":
                    dofs_y = locate_dofs_on_sub(self.W.sub(0).sub(1), region_id)
                    bcs_mixed.append(
                        dolfinx.fem.dirichletbc(
                            dolfinx.default_scalar_type(0.0), dofs_y, self.W.sub(0).sub(1)
                        )
                    )

                elif bc_type == "Clamp_z":
                    val = float(bc.get("value", 0.0))  # allowing GPS
                    dofs_z = locate_dofs_on_sub(self.W.sub(0).sub(2), region_id)
                    bcs_mixed.append(
                        dolfinx.fem.dirichletbc(
                            dolfinx.default_scalar_type(val), dofs_z, self.W.sub(0).sub(2)
                        )
                    )

        # -- NONLINEAR PROBLEM & SOLVER (FEniCSx ≥ 0.10.0) --
        petsc_options = {
            "snes_type": "newtonls",
            "snes_linesearch_type": "none",
            "snes_monitor": None,
            "snes_atol": 1e-8,
            "snes_rtol": 1e-8,
            "snes_stol": 1e-8,
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        }
        problem = NonlinearProblem(
            F,
            self.sol_mixed,
            bcs=bcs_mixed,
            petsc_options=petsc_options,
            petsc_options_prefix="elasticity",
        )

        print("[INFO] Solving monolithic nonlinear problem...")
        problem.solve()

        snes = problem.solver
        iters = snes.getIterationNumber()
        converged = snes.getConvergedReason()

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
        u_vec = self.u.x.array.reshape(-1, self.tdim)
        umag = np.linalg.norm(u_vec, axis=1)
        print(f"Min/Max displacement magnitude: {umag.min():.2e} / {umag.max():.2e}")

    def solve_staggered(
        self, max_iter=20, stag_tol_th=1e-3, stag_tol_mech=1e-3, rtol_th=1e-6, rtol_mech=1e-6
    ):
        print(f"  → Max iterations              : {max_iter}")
        print(f"  → Staggering tolerance |ΔT|   : {stag_tol_th:.1e}")
        print(f"  → Staggering tolerance |Δu|   : {stag_tol_mech:.1e}")
        print(f"  → Relative tolerance th       : {rtol_th:.1e}")
        print(f"  → Relative tolerance mech     : {rtol_mech:.1e}")

        if self.on.get("thermal", False):
            T_i = dolfinx.fem.Function(self.V_t)
            T_i.x.array[:] = self.T.x.array  # Start with the initial/previous converged state
            T_old = dolfinx.fem.Function(
                self.V_t
            )  # Functions to store the previous state for convergence check

            bcs_t = []
            [bcs_t.extend(bc_list) for bc_list in self.dirichlet_thermal.values()]

        if self.on.get("mechanical", False):
            u_i = dolfinx.fem.Function(self.V_m)
            u_i.x.array[:] = self.u.x.array
            u_old = dolfinx.fem.Function(self.V_m)

            bcs_m = []
            [bcs_m.extend(bc_list) for bc_list in self.dirichlet_mechanical.values()]

        self.dx_tags = {
            tag: ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag
            )
            for tag in np.unique(self.cell_tags.values)
        }

        self.ds_tags = {
            id_: ufl.Measure(
                "ds", domain=self.mesh, subdomain_data=self.facet_tags, subdomain_id=id_
            )
            for id_ in np.unique(self.facet_tags.values)
        }

        prev_res_T = None
        prev_res_u = None

        # --- Staggering loop ---
        for iteration in range(max_iter):
            print(f"\n--- Staggering iteration {iteration+1}/{max_iter} ---")

            # 1) THERMAL
            if self.on.get("thermal", False):
                T_old.x.array[:] = T_i.x.array

                print("\n[INFO] Assembling thermal problem...")
                u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
                a_t = 0
                L_t = 0

                # --- Volume integrals ---
                for label, material in self.materials.items():
                    tag = self.label_map[label]
                    dx = self.dx_tags[tag]

                    print(f"\n  Building weak form, volume integrals (dx) for {label}, tag = {tag}")

                    # Thermal conductivity
                    k = material["k"]

                    # Output the average value of the thermal conductivity
                    # if isinstance(k, (int, float)):
                    #     k_avg = k
                    # else:
                    #     V0 = dolfinx.fem.functionspace(self.mesh, ("CG", 1))
                    #     k_func = dolfinx.fem.Function(V0, name=f"k_{label}")
                    #     expr = dolfinx.fem.Expression(k, V0.element.interpolation_points)
                    #     k_func.interpolate(expr)

                    #     num = dolfinx.fem.assemble_scalar(dolfinx.fem.form(k_func * dx))
                    #     den = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx))
                    #     k_avg = num / den

                    # print(f"  Average conductivity for {label} = {k_avg:.3f} W/m·K")

                    a_t += k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
                    L_t += self.q_third * v_t * dx

                    # q_third value in the volume "label" / "tag"
                    dofs = self.mgr.locate_domain_dofs(label=self.label_map[label], V=self.V_t)
                    q_vals = self.q_third.x.array[dofs]
                    print(
                        f"  → q_third[{label}](W/m3) min = {q_vals.min():.2e}, max = {q_vals.max():.2e}, mean = {q_vals.mean():.2e}"
                    )

                # Neumann
                for label in self.materials:
                    for bc_info in self.neumann_thermal[label]:
                        print(f"  Applying flux on subdomain id = {bc_info['id']}")
                        ds_neumann = self.ds_tags[bc_info["id"]]
                        L_t += bc_info["value"] * v_t * ds_neumann

                # Gap (Robin)
                h_gap = self.set_gap_conductance(T_i)

                for label in self.materials:
                    for bc_info in self.robin_thermal[label]:
                        print(f"  Applying thermal Robin BC on subdomain id = {bc_info['id']}")

                        region_id = bc_info["id"]  # e.g., interface 1 --> 2
                        pair_region = bc_info["pair"]  # e.g., interface 1 <-- 2

                        ds_interface = self.ds_tags[region_id]  # e.g., 1 interface surface

                        # 1) project
                        T_other = dolfinx.fem.Function(self.V_t)  # e.g., T_cyl_2
                        dofs_here = self.mgr.locate_facets_dofs(
                            region_id, self.V_t
                        )  # e.g., dofs of cyl_1_interface
                        dofs_other = self.mgr.locate_facets_dofs(
                            self.label_map[pair_region], self.V_t
                        )  # e.g., dofs of cyl_2_interface

                        T_other.x.array[dofs_here] = T_i.x.array[
                            dofs_other
                        ]  # e.g., T_cyl_2 projected on cyl_1_dofs

                        # 2) mean
                        # T_other = T_i.x.array[dofs_other].mean()

                        # -) bias
                        # if label == 'cyl_1':
                        #     print(f"Biasing cyl_2 temperature")
                        #     T_other.x.array[dofs_here[:200]] += 500

                        a_t += h_gap * u_t * v_t * ds_interface
                        L_t += h_gap * T_other * v_t * ds_interface

                        print(
                            f"  [INFO] Gap Robin BC between '{label}' (region={region_id}) and '{pair_region}' (region={self.label_map[pair_region]})"
                        )

                        # Print heat flux
                        # self.heat_flux(T_i)

                # --- Solve the system ---
                if self.thermal_options.get("solver") == "linear":
                    print(f"  Linear solver")
                    petsc_opts_thermal = self.get_solver_options(
                        solver_type=self.thermal_options.get("linear_solver", None),
                        physics="thermal",
                        rtol=rtol_th,
                    )
                    problem_t = LinearProblem(
                        a_t,
                        L_t,
                        bcs=bcs_t,
                        u=T_i,
                        petsc_options=petsc_opts_thermal,
                        petsc_options_prefix="thermal_",
                    )
                    problem_t.solve()

                elif self.thermal_options.get("solver") == "non-linear":
                    print(f"  Non-linear solver")
                    print(f"  [ERROR] Non-linear thermal solver not yet implemented.")

                # Check T_fo
                # dofs = self.locateFacetDofs(self.label_map["lateral"], self.V_t)
                # T = T_i.x.array[dofs].mean()
                # print(f"  T_lateral = {T:.2f} K")

                # Apply relaxation
                T_i.x.array[:] = self.relax_T * T_i.x.array + (1 - self.relax_T) * T_old.x.array

            # --. MECHANICAL SUB-PROBLEM --..
            if self.on.get("mechanical", False):
                u_old.x.array[:] = u_i.x.array

                print("\n[INFO] Assembling mechanical problem...")

                u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
                a_m, L_m = 0, 0  # linear formulation a(u,v) = F(v)
                F_m = 0  # non-linear formulation F(u,v) = 0

                # --- Volume integrals ---
                for label, material in self.materials.items():
                    tag = self.label_map[label]
                    print(f"  Building weak form, volume integrals (dx) for {label}, tag = {tag}")
                    dx = self.dx_tags[tag]

                    rho = dolfinx.default_scalar_type(material["rho"])
                    g = dolfinx.default_scalar_type(self.g)
                    body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -rho * g))

                    if self.mech_options.get("solver") == "linear":
                        stress_tensor = self.sigma_mech(u_m, material)  # in mechanical_model

                        a_m += ufl.inner(stress_tensor, self.epsilon(v_m)) * dx
                        L_m += ufl.dot(body_force, v_m) * dx
                        if self.on.get("thermal", False):
                            L_m -= ufl.inner(self.sigma_th(T_i, material), self.epsilon(v_m)) * dx

                    elif self.mech_options.get("solver") == "non-linear":
                        stress_tensor = self.sigma_mech(u_i, material)
                        F_m += (
                            ufl.inner(stress_tensor, self.epsilon(v_m)) * dx
                            - ufl.dot(body_force, v_m) * dx
                        )

                # --. Neumann BCs (traction) --..
                for label in self.materials:
                    for bc_info in self.traction[label]:
                        print(f"  Applying mechanical traction on subdomain id = {bc_info['id']}")
                        ds = self.ds_tags[bc_info["id"]]
                        if self.mech_options.get("solver") == "linear":
                            L_m += ufl.dot(bc_info["value"], v_m) * ds
                        elif self.mech_options.get("solver") == "non-linear":
                            F_m -= ufl.dot(bc_info["value"], v_m) * ds

                # --. Neumann BCs (penalties) --..
                for label in self.materials:
                    for bc_info in self.clamp_r[label]:
                        alpha = bc_info["penalty"]
                        val = bc_info["value"]
                        print(
                            f"  Applying Clamp_r (weak penalty) on region id = {bc_info["id"]} (α = {alpha:.2e})"
                        )

                        ds = self.ds_tags[bc_info["id"]]
                        n = ufl.FacetNormal(self.mesh)

                        if self.mech_options.get("solver") == "linear":
                            a_m += alpha * ufl.dot(u_m, n) * ufl.dot(v_m, n) * ds
                            if abs(val) > 1e-16:
                                L_m += alpha * val * ufl.dot(v_m, n) * ds

                        elif self.mech_options.get("solver") == "non-linear":
                            F_m += alpha * ufl.dot(u_m, n) * ufl.dot(v_m, n) * ds
                            if abs(val) > 1e-16:
                                F_m -= alpha * val * ufl.dot(v_m, n) * ds

                # --. check_/ debug mechanical_constraints --..
                if self.mech_debug:
                    if not hasattr(self, "_mech_diagnostic_done"):
                        print("[DEBUG] Mechanical diagnostic enabled (mechanical.debug=True)")
                        try:
                            from z3st.core.diagnostic import (
                                check_fixed_dofs,
                                check_mechanical_constraints,
                                debug_dirichlet_mechanical,
                            )

                            debug_dirichlet_mechanical(self)
                            check_fixed_dofs(self, "xmin")
                            check_mechanical_constraints(
                                self, eig_check=False, plot=False, tol=1e-5
                            )
                        except Exception as e:
                            print(f"[WARN] Mechanical diagnostic failed: {e}")
                        finally:
                            self._mech_diagnostic_done = True

                # --. solve --..
                if self.mech_options.get("solver") == "linear":
                    print(f"  Linear solver")
                    petsc_opts_mech = self.get_solver_options(
                        solver_type=self.mech_options.get("linear_solver", None),
                        physics="mechanical",
                        rtol=rtol_mech,
                    )
                    problem_m = LinearProblem(
                        a_m,
                        L_m,
                        bcs=bcs_m,
                        u=u_i,
                        petsc_options=petsc_opts_mech,
                        petsc_options_prefix="mechanical_",
                    )
                    problem_m.solve()

                elif self.mech_options.get("solver") == "non-linear":
                    print(f"  Non-linear solver")

                    petsc_opts_mech = self.get_solver_options(
                        solver_type=self.mech_options.get("linear_solver", None),
                        physics="mechanical",
                        rtol=rtol_mech,
                    )

                    problem_m = NonlinearProblem(
                        F_m,
                        u_i,
                        bcs=bcs_m,
                        petsc_options=petsc_opts_mech,
                        petsc_options_prefix="elasticity",
                    )

                    problem_m.solve()

                    snes = problem_m.solver
                    iters = snes.getIterationNumber()
                    converged = snes.getConvergedReason()

                    assert converged

                    if not converged:
                        raise RuntimeError(
                            f"[mechanical] NewtonSolver failed after {iters} iterations"
                        )

                # --. Apply relaxation --..
                u_i.x.array[:] = self.relax_u * u_i.x.array + (1 - self.relax_u) * u_old.x.array

            # --. CONVERGENCE CHECK and LOGGING --..
            print(f"\nConvergence check")

            rel_norm_dT = 0.0
            rel_norm_du = 0.0

            # Thermal convergence
            conv_th = True
            if self.on.get("thermal", False):
                T_i.x.scatter_forward()
                T_old.x.scatter_forward()

                vec_T_new = T_i.x.petsc_vec
                vec_T_old = T_old.x.petsc_vec

                # Compute difference vector: diff_T = T_now - T_old
                diff_T = vec_T_new.copy()
                diff_T.axpy(-1.0, vec_T_old)

                # Compute norms
                norm_dT = diff_T.norm(PETSc.NormType.NORM_2)
                norm_T = vec_T_new.norm(PETSc.NormType.NORM_2)
                rel_norm_dT = norm_dT / norm_T if norm_T > 1e-12 else norm_dT

                # print(f"  ||ΔT|| = {norm_dT:.3e}, ||ΔT||/||T|| = {rel_norm_dT:.3e}")

                if self.th_convergence == "norm":
                    print(f"  ||ΔT|| = {norm_dT:.3e}")
                    conv_th = norm_dT < stag_tol_th or not self.on.get("thermal", False)
                elif self.th_convergence == "rel_norm":
                    print(f"  ||ΔT||/||T|| = {rel_norm_dT:.3e}")
                    conv_th = rel_norm_dT < stag_tol_th or not self.on.get("thermal", False)

                # --. ADAPTIVE relax_T --..
                if self.relax_adaptive:
                    resT_curr = norm_dT if self.th_convergence == "norm" else rel_norm_dT
                    if prev_res_T is not None:
                        if resT_curr < prev_res_T:
                            self.relax_T = min(self.relax_T * self.relax_growth, self.relax_max)
                        else:
                            self.relax_T = max(self.relax_T * self.relax_shrink, self.relax_min)
                    prev_res_T = resT_curr

                    print(f"  [adaptive] relax_T={self.relax_T:.2f}")

            # Mechanical convergence
            conv_mech = True
            if self.on.get("mechanical", False):
                u_i.x.scatter_forward()
                u_old.x.scatter_forward()

                # Get the PETSc vectors
                vec_u_new = u_i.x.petsc_vec
                vec_u_old = u_old.x.petsc_vec

                # Compute difference vector: diff_u = u_now - u_old
                diff_u = vec_u_new.copy()
                diff_u.axpy(-1.0, vec_u_old)

                # Compute norms
                norm_du = diff_u.norm(PETSc.NormType.NORM_2)
                norm_u = vec_u_new.norm(PETSc.NormType.NORM_2)
                rel_norm_du = norm_du / norm_u if norm_u > 1e-12 else norm_du

                # print(f"  ||Δu|| = {norm_du:.3e}, ||Δu||/||u|| = {rel_norm_du:.3e}")

                if self.mech_convergence == "norm":
                    print(f"  ||Δu|| = {norm_du:.3e}")
                    conv_mech = norm_du < stag_tol_mech or not self.on.get("mechanical", False)
                elif self.mech_convergence == "rel_norm":
                    print(f"  ||Δu||/||u|| = {rel_norm_du:.3e}")
                    conv_mech = rel_norm_du < stag_tol_mech or not self.on.get("mechanical", False)

                # --. ADAPTIVE relax_u --..
                if self.relax_adaptive:
                    resU_curr = norm_du if self.mech_convergence == "norm" else rel_norm_du
                    if prev_res_u is not None:
                        if resU_curr < prev_res_u:
                            self.relax_u = min(self.relax_u * self.relax_growth, self.relax_max)
                        else:
                            self.relax_u = max(self.relax_u * self.relax_shrink, self.relax_min)
                    prev_res_u = resU_curr

                    print(f"  [adaptive] relax_u={self.relax_u:.2f}")

            if conv_mech and conv_th:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")
                if self.on.get("thermal", False):
                    self.T.x.array[:] = T_i.x.array
                    T_min, T_max = self.T.x.array.min(), self.T.x.array.max()
                    print(f"Global min/max temperature: {T_min:.2f} / {T_max:.2f} K")

                if self.on.get("mechanical", False):
                    self.u.x.array[:] = u_i.x.array
                    u_vec = u_i.x.array.reshape(-1, self.tdim)
                    umag = np.sqrt((u_vec**2).sum(axis=1))
                    print(
                        f"Global min/max displacement magnitude: {umag.min():.2e} / {umag.max():.2e} m"
                    )

                return True

        print("\n[WARNING] Staggered solver did not converge. Using last iteration state.")

        return False
