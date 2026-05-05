# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import numpy as np
import ufl
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI
from petsc4py import PETSc


class Solver:
    def __init__(self):
        print("[Solver] initializer")
        solver_settings = self.input_file.get("solver_settings", {})

        self.coupling = solver_settings.get("coupling", "staggered")

        self.relax_T = float(solver_settings.get("relax_T", 0.9))
        self.relax_u = float(solver_settings.get("relax_u", 0.4))
        self.relax_D = float(solver_settings.get("relax_D", 0.4))

        print("  Applied relaxation factor:")
        print(f"  → Temperature  : {self.relax_T}")
        print(f"  → Displacement : {self.relax_u}")
        print(f"  → Damage       : {self.relax_D}")

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

        physics: "thermal", "mechanical" or "damage".
        """
        if physics not in ["thermal", "mechanical", "damage"]:
            raise ValueError(
                f"Unknown physics '{physics}'. Must be 'thermal', 'mechanical' or 'damage'."
            )

        # 1) KSP type
        if physics in ["thermal", "damage"]:
            ksp_type = "cg"  # SPD
        else:  # mechanical
            ksp_type = "gmres"  # non-symmetric

        # 2) Preconditioner
        if solver_type == "direct_mumps":
            return {
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
            }

        elif solver_type == "iterative_amg":
            coarse_eq_limit = 500 if physics == "thermal" else 1000
            return {
                "ksp_type": ksp_type,
                "pc_type": "gamg",
                "ksp_rtol": rtol,
                "pc_gamg_coarse_eq_limit": coarse_eq_limit,
            }

        elif solver_type == "iterative_hypre":
            return {
                "ksp_type": ksp_type,
                "pc_type": "hypre",
                "pc_hypre_type": "boomeramg",
                "ksp_rtol": rtol,
            }

        else:
            raise ValueError(f"Unknown solver_type '{solver_type}'.")

    def _build_measures(self):
        """Build dx_tags and ds_tags measures with axisymmetric and cartesian support."""
        x = ufl.SpatialCoordinate(self.mesh)
        regime = self.regime
        
        # Integration weight logic:
        # - Axisymmetric: 2*pi*r
        # - Cartesian 2D: 1.0 (Area)
        # - 3D: 1.0 (Volume)

        if regime == "axisymmetric":
            self.weight = 2.0 * ufl.pi * x[0]
        elif regime == "2d":
            self.weight = 1.0
        else:
            self.weight = 1.0

        metadata = {}
        if getattr(self, "q_degree", None) is not None:
            metadata = {"quadrature_degree": self.q_degree, "quadrature_scheme": "default"}
            print(f"  [Solver] Using quadrature degree {self.q_degree} for integration measures.")

        self.dx_tags = {
            tag: ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag, metadata=metadata
            )
            for tag in np.unique(self.cell_tags.values)
        }

        self.ds_tags = {
            id_: ufl.Measure(
                "ds", domain=self.mesh, subdomain_data=self.facet_tags, subdomain_id=id_
            )
            for id_ in np.unique(self.facet_tags.values)
        }

    def _thermal_step(self, T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T):

        analysis = self.th_cfg.get("analysis", "stationary")

        # Transient mode with dt=0: preserve IC, only apply BCs
        if analysis == "transient" and self.dt <= 0:
            print("\n[INFO] Transient thermal: dt=0 → preserving initial condition (applying BCs only)")
            bcs_thermal_actual = [
                bc["value"] if isinstance(bc, dict) else bc
                for _, bc_list in self.dirichlet_thermal.items()
                for bc in bc_list
            ]
            dolfinx.fem.set_bc(T_new.x.array, bcs_thermal_actual)
            T_new.x.scatter_forward()
            return True, 0.0, 0.0, prev_res_T

        T_old.x.array[:] = T_new.x.array

        print("\n[INFO] Assembling thermal problem...")
        u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
        a_t = 0
        L_t = 0

        w = self.weight
        dt = self.dt
        transient = analysis == "transient"

        # Volume integrals
        for label, material in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]

            print(f"\n  Building weak form, volume integrals (dx) for {label}, tag = {tag}")
            k = material["k"]
            rho = material["rho"]
            cp = material["cp"]
            rho_cp = rho * cp

            # Diffusion + source (always present)
            a_t += w * k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
            L_t += w * self.q_third * v_t * dx

            # Mass term (backward Euler, only for transient)
            # self.T holds the converged temperature from the previous time step (T^n)
            # T_old is used for stagger relaxation only, not for backward Euler
            if transient:
                rho_cp_dt = rho_cp / dt
                a_t += w * rho_cp_dt * u_t * v_t * dx
                L_t += w * rho_cp_dt * self.T * v_t * dx

            dofs = self.mgr.locate_domain_dofs(label=self.label_map[label], V=self.V_t)
            q_vals = self.q_third.x.array[dofs]
            print(
                f"  → q_third[{label}](W/m3) min = {q_vals.min():.2e}, "
                f"max = {q_vals.max():.2e}, mean = {q_vals.mean():.2e}"
            )

        # Neumann
        for label in self.materials:
            for bc_info in self.neumann_thermal[label]:
                print(f"  Applying flux on subdomain id = {bc_info['id']}")
                ds_neumann = self.ds_tags[bc_info["id"]]
                L_t += w * (-bc_info["value"]) * v_t * ds_neumann

        # Robin BCs (gap or convective)
        h_gap = self.set_gap_conductance(T_new)

        for label in self.materials:
            for bc_info in self.robin_thermal[label]:
                region_id = bc_info["id"]
                ds_robin = self.ds_tags[region_id]

                if "pair" in bc_info:
                    # Gap mode: h from gap model, T_ext from paired subdomain
                    pair_region = bc_info["pair"]
                    T_other = dolfinx.fem.Function(self.V_t)
                    dofs_here = self.mgr.locate_facets_dofs(region_id, self.V_t)
                    dofs_other = self.mgr.locate_facets_dofs(self.label_map[pair_region], self.V_t)
                    T_other.x.array[dofs_here] = T_new.x.array[dofs_other]

                    a_t += w * h_gap * u_t * v_t * ds_robin
                    L_t += w * h_gap * T_other * v_t * ds_robin
                    print(f"  Robin (gap) BC on region {region_id}, paired with '{pair_region}'")

                else:
                    # Convective mode: fixed h_conv and T_ext
                    h_conv = bc_info["h_conv"]
                    T_ext = bc_info["T_ext"]
                    a_t += w * h_conv * u_t * v_t * ds_robin
                    L_t += w * h_conv * T_ext * v_t * ds_robin
                    print(f"  Robin (convective) BC on region {region_id}: h={h_conv:.1f} W/(m²·K), T_ext={T_ext:.1f} K")

        # Extract actual DirichletBC objects (handles both dict and direct BCs)
        bcs_thermal_actual = [
            bc["value"] if isinstance(bc, dict) else bc
            for _, bc_list in self.dirichlet_thermal.items()
            for bc in bc_list
        ]

        # Solve
        if self.th_cfg["solver"] == "linear":
            print("  Linear solver")
            petsc_opts_thermal = self.get_solver_options(
                solver_type=self.th_cfg["linear_solver"],
                physics="thermal",
                rtol=rtol_th,
            )
            problem_t = dolfinx.fem.petsc.LinearProblem(
                a_t,
                L_t,
                bcs=bcs_t,
                u=T_new,
                petsc_options=petsc_opts_thermal,
                petsc_options_prefix="thermal_",
            )
            dolfinx.fem.set_bc(T_new.x.array, bcs_t)
            problem_t.solve()
            print(f"  T_new: min={T_new.x.array.min():.2f} K, max={T_new.x.array.max():.2f} K, mean={T_new.x.array.mean():.2f} K")
            if transient:
                print(f"  T^n (self.T): min={self.T.x.array.min():.2f} K, max={self.T.x.array.max():.2f} K")
        else:
            print("  [ERROR] Non-linear thermal solver not yet implemented.")

        # Relax
        T_new.x.array[:] = self.relax_T * T_new.x.array + (1.0 - self.relax_T) * T_old.x.array
        dolfinx.fem.set_bc(T_new.x.array, bcs_thermal_actual)

        # Convergenza (norma o rel_norm)
        T_new.x.scatter_forward()
        T_old.x.scatter_forward()

        vec_T_new = T_new.x.petsc_vec
        vec_T_old = T_old.x.petsc_vec

        diff_T = vec_T_new.copy()
        diff_T.axpy(-1.0, vec_T_old)

        norm_dT = diff_T.norm(PETSc.NormType.NORM_2)
        norm_T = vec_T_new.norm(PETSc.NormType.NORM_2)
        rel_norm_dT = norm_dT / norm_T if norm_T > 1e-15 else norm_dT

        if self.th_cfg["convergence"] == "norm":
            print(f"  ||ΔT|| = {norm_dT:.3e}")
            conv_th = norm_dT < stag_tol_th
            res_curr = norm_dT
        else:
            print(f"  ||ΔT||/||T|| = {rel_norm_dT:.3e}")
            conv_th = rel_norm_dT < stag_tol_th
            res_curr = rel_norm_dT

        if self.relax_adaptive:
            if prev_res_T is not None:
                ema_alpha = 0.3
                ema_T = ema_alpha * res_curr + (1 - ema_alpha) * prev_res_T
                if res_curr < ema_T:
                    self.relax_T = min(self.relax_T * self.relax_growth, self.relax_max)
                else:
                    self.relax_T = max(self.relax_T * self.relax_shrink, self.relax_min)
                prev_res_T = ema_T
            else:
                prev_res_T = res_curr
            print(f"  [adaptive] relax_T={self.relax_T:.2f}")

        return conv_th, norm_dT, rel_norm_dT, prev_res_T

    def _mechanical_step(self, u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current):

        u_old.x.array[:] = u_new.x.array
        print("\n[INFO] Assembling mechanical problem...")

        # --- update step-dependent displacement ---
        for _, bc_list in self.dirichlet_mechanical.items():
            for bc in bc_list:
                # Skip BCs that are Clamp, Slip, etc. (not yet step-dependent)
                if not isinstance(bc, dict):
                    continue

                raw = bc.get("raw", None)
                if isinstance(raw, list):
                    idx = min(self.current_step, len(raw) - 1)
                    val = raw[idx]
                    bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)
                    print(f"  [INFO] Updating Displacement Dirichlet on region {bc['id']} → {val}")

        # --- update step-dependent tractions ---
        for _, bc_list in self.traction.items():
            for bc in bc_list:
                raw = bc.get("raw", None)
                
                if isinstance(raw, list):
                    idx = min(self.current_step, len(raw) - 1)
                    val = raw[idx]
                elif isinstance(raw, (int, float)):
                    val = raw
                else:
                    raise RuntimeError("Invalid traction 'raw' format")

                bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)
                print(f"  [INFO] Updating traction on region {bc['id']} → {val} Pa")

                regime = self.regime
                if regime in ["axisymmetric", "2d"]:
                    n_vec = ufl.as_vector([self.normal[0], self.normal[1]])
                else:
                    n_vec = self.normal
                
                bc["value"] = bc["const"] * n_vec

        u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
        a_m, L_m = 0, 0
        F_m = 0
        w = self.weight

        for label, material in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]
            print(f"  Building weak form, volume integrals (dx) for {label}, tag = {tag}")

            rho = dolfinx.default_scalar_type(material["rho"])
            g = dolfinx.default_scalar_type(self.g)

            regime = self.regime
            if regime == "axisymmetric" or regime == "2d":
                # 2D: (F_r, F_z) or (F_x, F_y)
                body_force = dolfinx.fem.Constant(self.mesh, (0.0, -rho * g))
            else:
                # 3D: (F_x, F_y, F_z)
                body_force = dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -rho * g))

            if self.mech_cfg["solver"] == "linear":
                sigma = self.sigma_mech(u_m, material)
                a_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx
                L_m += w * ufl.dot(body_force, v_m) * dx
                if self.on.get("thermal", False):
                    L_m -= w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx
            else:
                mode = material.get("constitutive_mode", "lame")
                if mode == "hyperelastic":
                    F_m += self.hyperelastic_residual(u_new, v_m, material, dx, w)
                    F_m -= w * ufl.dot(body_force, v_m) * dx
                    if self.on.get("thermal", False):
                        F_m += w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx
                else:
                    sigma = self.sigma_mech(u_new, material)
                    F_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx - w * ufl.dot(body_force, v_m) * dx

        # Traction BCs
        for label in self.materials:
            for bc_info in self.traction[label]:
                print(f"  Applying mechanical traction on subdomain id = {bc_info['id']}")
                ds = self.ds_tags[bc_info["id"]]
                if self.mech_cfg["solver"] == "linear":
                    L_m += w * ufl.dot(bc_info["value"], v_m) * ds
                else:
                    F_m -= w * ufl.dot(bc_info["value"], v_m) * ds

        # --- Extract actual DirichletBC objects (handles both dict and direct BCs) ---
        bcs_mech = [
            bc["value"] if isinstance(bc, dict) else bc
            for _, bc_list in self.dirichlet_mechanical.items()
            for bc in bc_list
        ]

        # Solve
        if self.mech_cfg["solver"] == "linear":
            print("  Linear solver")
            petsc_opts_mech = self.get_solver_options(
                solver_type=self.mech_cfg["linear_solver"],
                physics="mechanical",
                rtol=rtol_mech,
            )
            problem_m = dolfinx.fem.petsc.LinearProblem(
                a_m,
                L_m,
                bcs=bcs_mech,
                u=u_new,
                petsc_options=petsc_opts_mech,
                petsc_options_prefix="mechanical_",
            )
            dolfinx.fem.set_bc(u_new.x.array, bcs_mech)
            problem_m.solve()
        else:
            print("  Non-linear solver (SNES Newton)")
            linear_solver = self.mech_cfg.get("linear_solver", "direct_mumps")

            # SNES Newton options + inner linear solver
            if linear_solver == "direct_mumps":
                petsc_opts_mech = {
                    "snes_type": "newtonls",
                    "snes_linesearch_type": "basic",
                    "snes_atol": rtol_mech,
                    "snes_rtol": rtol_mech,
                    "snes_max_it": int(self.mech_cfg.get("snes_max_it", 50)),
                    "ksp_type": "preonly",
                    "pc_type": "lu",
                    "pc_factor_mat_solver_type": "mumps",
                }
            else:
                # Iterative inner solver (AMG / HYPRE)
                ksp_opts = self.get_solver_options(
                    solver_type=linear_solver,
                    physics="mechanical",
                    rtol=rtol_mech,
                )
                petsc_opts_mech = {
                    "snes_type": "newtonls",
                    "snes_linesearch_type": "bt",
                    "snes_atol": rtol_mech,
                    "snes_rtol": rtol_mech,
                    "snes_max_it": 100,
                    "snes_divergence_tolerance": 1e10,
                    **ksp_opts,
                }

            problem_m = NonlinearProblem(
                F_m,
                u_new,
                bcs=bcs_mech,
                petsc_options=petsc_opts_mech,
                petsc_options_prefix="elasticity_",
            )
            dolfinx.fem.set_bc(u_new.x.array, bcs_mech)
            problem_m.solve()

        # Relax
        u_new.x.array[:] = self.relax_u * u_new.x.array + (1 - self.relax_u) * u_old.x.array
        dolfinx.fem.set_bc(u_new.x.array, bcs_mech)

        # Convergence
        u_new.x.scatter_forward()
        u_old.x.scatter_forward()

        vec_u_new = u_new.x.petsc_vec
        vec_u_old = u_old.x.petsc_vec

        diff_u = vec_u_new.copy()
        diff_u.axpy(-1.0, vec_u_old)

        norm_du = diff_u.norm(PETSc.NormType.NORM_2)
        norm_u = vec_u_new.norm(PETSc.NormType.NORM_2)
        rel_norm_du = norm_du / norm_u if norm_u > 1e-15 else norm_du

        if self.mech_cfg["convergence"] == "norm":
            print(f"  ||Δu|| = {norm_du:.3e}")
            conv_mech = norm_du < stag_tol_mech
            res_curr = norm_du
        else:
            print(f"  ||Δu||/||u|| = {rel_norm_du:.3e}")
            conv_mech = rel_norm_du < stag_tol_mech
            res_curr = rel_norm_du

        if self.relax_adaptive:
            if prev_res_u is not None:
                ema_alpha = 0.3
                ema_u = ema_alpha * res_curr + (1 - ema_alpha) * prev_res_u
                if res_curr < ema_u:
                    self.relax_u = min(self.relax_u * self.relax_growth, self.relax_max)
                else:
                    self.relax_u = max(self.relax_u * self.relax_shrink, self.relax_min)
                prev_res_u = ema_u
            else:
                prev_res_u = res_curr
            print(f"  [adaptive] relax_u={self.relax_u:.2f}")

        return conv_mech, norm_du, rel_norm_du, prev_res_u

    def _damage_step(self, D_new, D_old, rtol_dmg, stag_tol_dmg, prev_res_D):

        D_old.x.array[:] = D_new.x.array
        
        w = self.weight
        
        lc = float(self.dmg_cfg["lc"])
        damage_type = self.dmg_cfg["type"]

        print(f"\n[INFO] Assembling damage ({damage_type}) problem...")

        u_d, v_d = ufl.TrialFunction(self.V_d), ufl.TestFunction(self.V_d)
        a_d, L_d = 0, 0

        bcs_d = []
        for mat_name in self.materials:
            for bc_entry in self.dirichlet_damage.get(mat_name, []):
                bcs_d.append(bc_entry["value"])

        for label, material in self.materials.items():
            print(
                f"Solving damage problem for '{label}' material"
            )
            
            tag = self.label_map[label]
            dx = self.dx_tags[tag]

            Gc = material["Gc"]
            sigma_c = material["sigma_c"]
            E = material["E"]

            if damage_type == "AT2":
                
                a_d += w*((self.H + 1.0) * u_d * v_d 
                + lc**2 * ufl.inner(ufl.grad(u_d), ufl.grad(v_d))) * dx
                L_d += w*self.H * v_d * dx

            elif damage_type == "AT1":

                # AT1 surface energy density (Pham-Marigo-Bourdin convention):
                #   E_s = (3*Gc/8) * integral( D/lc + lc * |grad D|^2 ) dx
                # Variation w.r.t. D yields the strong form
                #   -(3*Gc*lc/4) * Laplacian(D) + 2*H*D = 2*H - 3*Gc/(8*lc)
                # The constant -3*Gc/(8*lc) on the RHS produces the sharp
                # AT1 elastic threshold: no damage until 2H > 3*Gc/(8*lc),
                # i.e. sigma > sigma_c.
                #
                # Irreversibility (D_n+1 >= D_n) is enforced after the
                # linear solve via np.maximum(D_new, D_old). A symmetric
                # (D - D_old)^2 penalty in the weak form would freeze
                # D ~= D_old whenever the penalty coefficient dominates the
                # driving force 2H, which is the typical regime in
                # thermal-shock cases where 2H is small. The post-solve
                # max-projection is the correct one-sided enforcement.
                #
                # A small Tikhonov shift (diag_shift) stabilises the linear
                # system in cells where H = 0 (otherwise the LHS bilinear
                # form would be just the Laplacian, which is positive
                # semi-definite with a constant nullspace under natural BCs).
                cw = 8.0 / 3.0
                pref = Gc / cw                       # = 3*Gc/8  (surface density coeff)
                grad_coeff = 2.0 * pref * lc         # = 3*Gc*lc/4 (bilinear form coeff)
                diag_shift = 1.0e-8 * (Gc / lc)

                print(f"  - Material '{label}': AT1 solve. Gc={Gc:.2e}, sigma_c={sigma_c:.2e}")

                a_d += w*(2.0 * self.H + diag_shift) * u_d * v_d * dx + \
                    grad_coeff * ufl.inner(ufl.grad(u_d), ufl.grad(v_d)) * dx

                L_d += w*(2.0 * self.H - (pref / lc)) * v_d * dx

        petsc_opts_damage = self.get_solver_options(
            physics="damage",
            solver_type=self.dmg_cfg["linear_solver"],
            rtol=rtol_dmg,
        )

        problem_d = dolfinx.fem.petsc.LinearProblem(
            a_d,
            L_d,
            bcs=bcs_d,
            u=D_new,
            petsc_options=petsc_opts_damage,
            petsc_options_prefix="damage_",
        )
        problem_d.solve()
        
        if damage_type == "AT1":
            # Clipping:
            D_new.x.array[:] = np.clip(D_new.x.array, 0.0, 1.0)

        # Relax
        D_new.x.array[:] = self.relax_D * D_new.x.array + (1 - self.relax_D) * D_old.x.array

        # irreversibility
        D_new.x.array[:] = np.maximum(D_new.x.array, D_old.x.array)
        D_new.x.array[:] = np.clip(D_new.x.array, 0.0, 1.0)

        # Convergence
        D_new.x.scatter_forward()
        D_old.x.scatter_forward()

        vec_D_new = D_new.x.petsc_vec
        vec_D_old = D_old.x.petsc_vec

        diff_D = vec_D_new.copy()
        diff_D.axpy(-1.0, vec_D_old)

        norm_dD = diff_D.norm(PETSc.NormType.NORM_2)
        norm_D = vec_D_new.norm(PETSc.NormType.NORM_2)
        rel_norm_dD = norm_dD / norm_D if norm_D > 1e-15 else norm_dD

        if self.dmg_cfg["convergence"] == "norm":
            print(f"  ||ΔD|| = {norm_dD:.3e}")
            conv_damage = norm_dD < stag_tol_dmg
            res_curr = norm_dD
        else:
            print(f"  ||ΔD||/||D|| = {rel_norm_dD:.3e}")
            conv_damage = rel_norm_dD < stag_tol_dmg
            res_curr = rel_norm_dD

        if self.relax_adaptive:
            if prev_res_D is not None:
                ema_alpha = 0.3
                ema_D = ema_alpha * res_curr + (1 - ema_alpha) * prev_res_D
                if res_curr < ema_D:
                    self.relax_D = min(self.relax_D * self.relax_growth, self.relax_max)
                else:
                    self.relax_D = max(self.relax_D * self.relax_shrink, self.relax_min)
                prev_res_D = ema_D
            else:
                prev_res_D = res_curr
            print(f"  [adaptive] relax_D={self.relax_D:.2f}")

        # Residual in L_inf norm
        res_D_inf = np.linalg.norm(D_new.x.array - D_old.x.array, ord=np.inf)
        print(f"  |ΔD|_∞ = {res_D_inf:.3e}")

        return conv_damage, norm_dD, rel_norm_dD, prev_res_D

    def _cluster_step(self, c_new, c_old, dt):
        """
        Solve the cluster dynamics step with mass conservation using DG.
        
        Solves: ∂c/∂t = -v ∂c/∂n + D ∂²c/∂n²

        Case v > 0: The clusters grow. The distribution moves to the right (larger n).
        Case v < 0: The clusters shrink (evaporation/dissolution). The distribution moves to the left (towards n=1).

        C_tot = ∫ c·n dn = constant
        
        DG formulation:
        - Upwind for advection
        - Symmetric Interior Penalty (SIPG) for diffusion
        """
        c_old.x.array[:] = c_new.x.array
        
        u_c, v_c = ufl.TrialFunction(self.V_c), ufl.TestFunction(self.V_c)
        
        # Parameters
        v_vel = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.v_cluster))
        D_diff = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.D_cluster))
        dt_c = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(dt))
        
        # Geometric info
        n = ufl.FacetNormal(self.mesh)
        h = ufl.CellDiameter(self.mesh)
        h_avg = (h('+') + h('-')) / 2.0
        
        # Penalty parameter for SIPG (diffusion)
        # For P1 elements, gamma = 10.0 is usually sufficient.
        gamma = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(10.0))

        # Péclet number diagnostics
        num_cells = self.mesh.topology.index_map(self.mesh.topology.dim).size_global
        local_coords = self.mesh.geometry.x[:, 0]
        if local_coords.size > 0:
            x_min_local = float(local_coords.min())
            x_max_local = float(local_coords.max())
        else:
            x_min_local = float("inf")
            x_max_local = float("-inf")
        x_min = self.mesh.comm.allreduce(x_min_local, op=MPI.MIN)
        x_max = self.mesh.comm.allreduce(x_max_local, op=MPI.MAX)
        L_domain = x_max - x_min
        h_cell = L_domain / num_cells
        v = abs(self.v_cluster)
        D = self.D_cluster
        pe = (v * h_cell) / (2 * D) if D > 0 else float('inf')

        print(f"  [Cluster DG] Peclet number Pe: {pe:.4e}")
        if pe > 1:
            print(f"    [INFO] Advection-dominated system. DG Upwind will provide stability.")

        # Variational form (Implicit uler + DG)
        # Mass matrix (time derivative)
        a = (u_c / dt_c) * v_c * ufl.dx
        L = (c_old / dt_c) * v_c * ufl.dx

        if dt > 0:
            # Advection term (Upwind)
            # Volume term
            a += - u_c * v_vel * v_c.dx(0) * ufl.dx
            
            # Interior facets
            v_n = v_vel * n[0]
            a += (ufl.avg(u_c * v_vel * n[0]) * ufl.jump(v_c) \
                 + 0.5 * abs(v_n('+')) * ufl.jump(u_c) * ufl.jump(v_c)) * ufl.dS
            
            # Boundary facets (outflow/inflow)
            a += ufl.conditional(v_n > 0, v_n * u_c * v_c, 0.0) * ufl.ds

            # Diffusion term (SIPG)
            a += D_diff * u_c.dx(0) * v_c.dx(0) * ufl.dx
            
            # Consistency and symmetry terms on interior facets
            a += - D_diff * ufl.avg(u_c.dx(0)) * ufl.jump(v_c, n[0]) * ufl.dS
            a += - D_diff * ufl.avg(v_c.dx(0)) * ufl.jump(u_c, n[0]) * ufl.dS
            
            # Penalty term on interior facets
            a += D_diff * (gamma / h_avg) * ufl.jump(u_c) * ufl.jump(v_c) * ufl.dS

            # Solve
            petsc_options = {
                "ksp_type": "gmres",
                "pc_type": "ilu",
                "ksp_rtol": 1e-12,
                "ksp_atol": 1e-15,
                "ksp_max_it": 1000,
            }

            problem = dolfinx.fem.petsc.LinearProblem(
                a, L, u=c_new, 
                petsc_options=petsc_options,
                petsc_options_prefix="cluster_"
            )

            problem.solve()
        else:
            print("  [Cluster] dt=0: skipping PDE solve.")
            c_new.x.array[:] = c_old.x.array

        # Mass conservation
        x = ufl.SpatialCoordinate(self.mesh)
        n_coord = x[0]
        
        C_tot_curr_new = self.mesh.comm.allreduce(
            dolfinx.fem.assemble_scalar(dolfinx.fem.form(c_new * n_coord * ufl.dx)),
            op=MPI.SUM,
        )
        
        if self.C_tot_target is not None:
            if abs(C_tot_curr_new) > 0.0:
                renorm_factor = self.C_tot_target / C_tot_curr_new
                c_new.x.array[:] *= renorm_factor
                
                print(f"  [Cluster] Mass conservation: target = {self.C_tot_target:.6e}")
                print(f"  [Cluster] Mass conservation: before = {C_tot_curr_new:.6e}")
                print(f"  [Cluster] Mass conservation: factor = {renorm_factor:.8f}")
            else:
                print("  [Cluster] Warning: mass is zero, cannot renormalize.")

        c_max_local = float(np.max(c_new.x.array)) if c_new.x.array.size > 0 else float("-inf")
        c_max = self.mesh.comm.allreduce(c_max_local, op=MPI.MAX)
        print(f"   [Diagnostics] Max density c_max: {c_max:.2f}")

    def solve_staggered(
        self,
        max_iter=20,
        dt=0.0,
        stag_tol_th=1e-3,
        stag_tol_mech=1e-3,
        stag_tol_dmg=1e-3,
        rtol_th=1e-6,
        rtol_mech=1e-6,
        rtol_dmg=1e-5,
    ):
        # Store dt as instance attribute for access in material models
        self.dt = dt

        print(f"  → Max iterations              : {max_iter}")
        print(f"  → Staggering tolerance |ΔT|   : {stag_tol_th:.1e}")
        print(f"  → Staggering tolerance |Δu|   : {stag_tol_mech:.1e}")
        print(f"  → Staggering tolerance |ΔD|   : {stag_tol_dmg:.1e}")
        print(f"  → Relative tolerance th       : {rtol_th:.1e}")
        print(f"  → Relative tolerance mech     : {rtol_mech:.1e}")
        print(f"  → Relative tolerance dmg      : {rtol_dmg:.1e}")

        # Build measures once
        self._build_measures()

        # Allocate local fields
        if self.on.get("thermal", False):
            T_new = dolfinx.fem.Function(self.V_t)
            T_new.x.array[:] = self.T.x.array
            T_old = dolfinx.fem.Function(self.V_t)

            bcs_t = []
            [bcs_t.extend(bc_list) for bc_list in self.dirichlet_thermal.values()]
        else:
            T_new = T_old = None
            bcs_t = []

        if self.on.get("mechanical", False):
            u_new = dolfinx.fem.Function(self.V_m)
            u_new.x.array[:] = self.u.x.array
            u_old = dolfinx.fem.Function(self.V_m)

            bcs_m = []
            [bcs_m.extend(bc_list) for bc_list in self.dirichlet_mechanical.values()]
        else:
            u_new = u_old = None
            bcs_m = []

        if self.on.get("damage", False):
            D_new = dolfinx.fem.Function(self.V_d)
            D_new.x.array[:] = self.D.x.array
            D_old = dolfinx.fem.Function(self.V_d)
        else:
            D_new = D_old = None

        if self.on.get("cluster", False):            
            c_new = dolfinx.fem.Function(self.V_c)
            c_new.x.array[:] = self.c.x.array
            self.c_n.x.array[:] = self.c.x.array
        else:
            c_new = None

        prev_res_T = None
        prev_res_u = None
        prev_res_D = None

        for iteration in range(max_iter):
            print(f"\n--- Staggering iteration {iteration+1}/{max_iter} ---")

            # Defaults
            conv_th = True
            conv_mech = True
            conv_damage = True

            # --. THERMAL STEP --..
            if self.on.get("thermal", False):
                conv_th, _, _, prev_res_T = self._thermal_step(
                    T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T
                )

            # --. MECHANICAL STEP --..
            if self.on.get("mechanical", False):
                conv_mech, _, _, prev_res_u = self._mechanical_step(
                    u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current=T_new
                )

            # --. DAMAGE STEP --..
            if self.on.get("damage", False):
                self.update_history(u_new)
                conv_damage, _, _, prev_res_D = self._damage_step(
                    D_new,
                    D_old,
                    rtol_dmg,
                    stag_tol_dmg,
                    prev_res_D,
                )
            
            # --. CLUSTER STEP --..
            if self.on.get("cluster", False):
                self._cluster_step(c_new, self.c_n, dt)

            # --.. GLOBAL CONVERGENCE --..
            print("\nConvergence check")

            if conv_th and conv_mech and conv_damage:
                print(f"\n[SUCCESS] Staggered solver converged in {iteration+1} iterations.")

                if self.on.get("thermal", False):
                    self.T.x.array[:] = T_new.x.array

                if self.on.get("mechanical", False):
                    self.u.x.array[:] = u_new.x.array

                if self.on.get("damage", False):
                    self.D.x.array[:] = D_new.x.array
                
                if self.on.get("cluster", False):
                    self.c.x.array[:] = c_new.x.array

                if self.on.get("mechanical", False) and self.on.get("plasticity", False):
                    self.update_plastic_history(u_new)

                return True

        # --. IF NOT CONVERGED --..
        print("\n[WARNING] Staggered solver did not converge. Using last iteration state.")

        if self.on.get("thermal", False):
            self.T.x.array[:] = T_new.x.array
        if self.on.get("mechanical", False):
            self.u.x.array[:] = u_new.x.array
        if self.on.get("damage", False):
            self.D.x.array[:] = D_new.x.array
        if self.on.get("cluster", False):
            self.c.x.array[:] = c_new.x.array

        return False
