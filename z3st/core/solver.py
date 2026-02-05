# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import numpy as np
import ufl
from dolfinx.fem.petsc import NonlinearProblem
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
        regime = self.mech_cfg.get("mechanical_regime", "3d").lower()

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

    def _thermal_step(self, T_new, T_old, bcs_t, rtol_th, stag_tol_th, prev_res_T):

        T_old.x.array[:] = T_new.x.array

        print("\n[INFO] Assembling thermal problem...")
        u_t, v_t = ufl.TrialFunction(self.V_t), ufl.TestFunction(self.V_t)
        a_t = 0
        L_t = 0

        w = self.weight

        # Volume integrals
        for label, material in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]

            print(f"\n  Building weak form, volume integrals (dx) for {label}, tag = {tag}")
            k = material["k"]

            a_t += w * k * ufl.inner(ufl.grad(u_t), ufl.grad(v_t)) * dx
            L_t += w * self.q_third * v_t * dx

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

        # Gap (Robin)
        h_gap = self.set_gap_conductance(T_new)

        for label in self.materials:
            for bc_info in self.robin_thermal[label]:
                print(f"  Applying thermal Robin BC on subdomain id = {bc_info['id']}")
                region_id = bc_info["id"]
                pair_region = bc_info["pair"]
                ds_interface = self.ds_tags[region_id]

                T_other = dolfinx.fem.Function(self.V_t)
                dofs_here = self.mgr.locate_facets_dofs(region_id, self.V_t)
                dofs_other = self.mgr.locate_facets_dofs(self.label_map[pair_region], self.V_t)
                T_other.x.array[dofs_here] = T_new.x.array[dofs_other]

                a_t += w * h_gap * u_t * v_t * ds_interface
                L_t += w * h_gap * T_other * v_t * ds_interface

                print(
                    f"  [INFO] Gap Robin BC between '{label}' "
                    f"(region={region_id}) and '{pair_region}' "
                    f"(region={self.label_map[pair_region]})"
                )

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
            problem_t.solve()
        else:
            print("  [ERROR] Non-linear thermal solver not yet implemented.")

        # Relax
        T_new.x.array[:] = self.relax_T * T_new.x.array + (1 - self.relax_T) * T_old.x.array

        # Convergenza (norma o rel_norm)
        T_new.x.scatter_forward()
        T_old.x.scatter_forward()

        vec_T_new = T_new.x.petsc_vec
        vec_T_old = T_old.x.petsc_vec

        diff_T = vec_T_new.copy()
        diff_T.axpy(-1.0, vec_T_old)

        norm_dT = diff_T.norm(PETSc.NormType.NORM_2)
        norm_T = vec_T_new.norm(PETSc.NormType.NORM_2)
        rel_norm_dT = norm_dT / norm_T if norm_T > 1e-12 else norm_dT

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
                if res_curr < prev_res_T:
                    self.relax_T = min(self.relax_T * self.relax_growth, self.relax_max)
                else:
                    self.relax_T = max(self.relax_T * self.relax_shrink, self.relax_min)
            prev_res_T = res_curr
            print(f"  [adaptive] relax_T={self.relax_T:.2f}")

        return conv_th, norm_dT, rel_norm_dT, prev_res_T

    def _mechanical_step(
        self, u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current
    ):

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
                    bc["const"].value = raw[idx]
                    print(f"  [INFO] Updating Dirichlet on region {bc['id']} → {raw[idx]}")

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

                regime = self.mech_cfg.get("mechanical_regime").lower()
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

            regime = self.mech_cfg.get("mechanical_regime").lower()
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
                sigma = self.sigma_mech(u_new, material)
                F_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx - ufl.dot(body_force, v_m) * dx

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
            problem_m.solve()
        else:
            print("  Non-linear solver")
            petsc_opts_mech = self.get_solver_options(
                solver_type=self.mech_cfg["linear_solver"],
                physics="mechanical",
                rtol=rtol_mech,
            )
            problem_m = NonlinearProblem(
                F_m,
                u_new,
                bcs=bcs_mech,
                petsc_options=petsc_opts_mech,
                petsc_options_prefix="elasticity",
            )
            problem_m.solve()

        # Relax
        u_new.x.array[:] = self.relax_u * u_new.x.array + (1 - self.relax_u) * u_old.x.array

        # Convergence
        u_new.x.scatter_forward()
        u_old.x.scatter_forward()

        vec_u_new = u_new.x.petsc_vec
        vec_u_old = u_old.x.petsc_vec

        diff_u = vec_u_new.copy()
        diff_u.axpy(-1.0, vec_u_old)

        norm_du = diff_u.norm(PETSc.NormType.NORM_2)
        norm_u = vec_u_new.norm(PETSc.NormType.NORM_2)
        rel_norm_du = norm_du / norm_u if norm_u > 1e-12 else norm_du

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
                if res_curr < prev_res_u:
                    self.relax_u = min(self.relax_u * self.relax_growth, self.relax_max)
                else:
                    self.relax_u = max(self.relax_u * self.relax_shrink, self.relax_min)
            prev_res_u = res_curr
            print(f"  [adaptive] relax_u={self.relax_u:.2f}")

        return conv_mech, norm_du, rel_norm_du, prev_res_u

    def _damage_step(self, D_new, D_old, rtol_dmg, stag_tol_dmg, prev_res_D, u_current):

        D_old.x.array[:] = D_new.x.array
        
        lc = float(self.dmg_cfg["lc"])
        damage_type = self.dmg_cfg["type"]

        print(f"\n[INFO] Assembling damage ({damage_type}) problem...")

        u_d, v_d = ufl.TrialFunction(self.V_d), ufl.TestFunction(self.V_d)
        a_d, L_d = 0, 0

        bcs_d = []
        for mat_name in self.materials:
            for bc_entry in self.dirichlet_damage.get(mat_name, []):
                bcs_d.append(bc_entry["value"])

        # Updating H field (according to damage model)
        self.update_history(u_current) 

        for label, material in self.materials.items():
            print(
                f"Solving damage problem for '{label}' material, with sigma_c = {material['sigma_c']*1e-6} MPa"
            )
            
            tag = self.label_map[label]
            dx = self.dx_tags[tag]

            if damage_type == "AT2":
                print(f"  - Material '{label}': Stress-based solve.")
                a_d += (1.0 + self.H) * u_d * v_d * dx + lc**2 * ufl.inner(
                    ufl.grad(u_d), ufl.grad(v_d)
                ) * dx
                L_d += self.H * v_d * dx

            elif damage_type == "AT1":

                E_mod = material["E"]
                sigma_c = material["sigma_c"]

                cw = 8.0 / 3.0 
                Gc = (8.0 * lc * sigma_c**2) / (3.0 * E_mod)
                pref = Gc / cw
                
                tol_ir = 0.05
                gamma_penalty = (Gc / lc) * (27.0 / (64.0 * tol_ir**2))
                w_act = (3.0 * Gc) / (8.0 * lc)

                print(f"  - Material '{label}': AT1 solve. Gc={Gc:.2e}, Gamma={gamma_penalty:.2e}")

                a_d += (2.0 * self.H + gamma_penalty) * u_d * v_d * dx + \
                    (pref * lc) * ufl.inner(ufl.grad(u_d), ufl.grad(v_d)) * dx
                                                
                L_d += (2.0 * self.H - (pref / lc) + gamma_penalty * D_old) * v_d * dx

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
        rel_norm_dD = norm_dD / norm_D if norm_D > 1e-12 else norm_dD

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
                if res_curr < prev_res_D:
                    self.relax_D = min(self.relax_D * self.relax_growth, self.relax_max)
                else:
                    self.relax_D = max(self.relax_D * self.relax_shrink, self.relax_min)
            prev_res_D = res_curr
            print(f"  [adaptive] relax_D={self.relax_D:.2f}")

        # Residual in L_inf norm
        res_D_inf = np.linalg.norm(D_new.x.array - D_old.x.array, ord=np.inf)
        print(f"  |ΔD|_∞ = {res_D_inf:.3e}")

        return conv_damage, norm_dD, rel_norm_dD, prev_res_D

    def solve_staggered(
        self,
        max_iter=20,
        stag_tol_th=1e-3,
        stag_tol_mech=1e-3,
        stag_tol_dmg=1e-3,
        rtol_th=1e-6,
        rtol_mech=1e-6,
        rtol_dmg=1e-5,
    ):
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

            # --. DAMAGE STEP --..
            if self.on.get("damage", False):
                conv_damage, _, _, prev_res_D = self._damage_step(
                    D_new,
                    D_old,
                    rtol_dmg,
                    stag_tol_dmg,
                    prev_res_D,
                    u_current=u_new,
                )

            # --. MECHANICAL STEP --..
            if self.on.get("mechanical", False):
                conv_mech, _, _, prev_res_u = self._mechanical_step(
                    u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current=T_new
                )

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

                return True

        # --- IF NOT CONVERGED ---
        print("\n[WARNING] Staggered solver did not converge. Using last iteration state.")

        if self.on.get("thermal", False):
            self.T.x.array[:] = T_new.x.array
        if self.on.get("mechanical", False):
            self.u.x.array[:] = u_new.x.array
        if self.on.get("damage", False):
            self.D.x.array[:] = D_new.x.array

        return False

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
        thermal_bcs_defs = self.boundary_conditions.get("thermal", {})
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
        mech_bcs_defs = self.boundary_conditions.get("mechanical", {})
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
