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
from mpi4py import MPI
from dolfinx import mesh, fem
import numpy as np
import yaml



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

    def _mechanical_step(self, u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current=None):
        """Étape mécanique standard mono-domaine pour le mixin Solver."""
        u_old.x.array[:] = u_new.x.array

        print("\n[INFO] Assembling mechanical problem...")

        # Mise à jour des conditions aux limites dépendantes du temps
        for _, bc_list in self.dirichlet_mechanical.items():
            for bc in bc_list:
                if isinstance(bc, dict) and "raw" in bc:
                    idx = min(self.current_step, len(bc["raw"]) - 1)
                    bc["const"].value = np.array(bc["raw"][idx], dtype=dolfinx.default_scalar_type)

        for _, bc_list in self.traction.items():
            for bc in bc_list:
                if "raw" in bc:
                    idx = min(self.current_step, len(bc["raw"]) - 1)
                    bc["const"].value = np.array(bc["raw"][idx], dtype=dolfinx.default_scalar_type)

        # Assemblage des formes
        w = self.weight
        v_m = ufl.TestFunction(self.V_m)
        a_m, L_m = 0, 0
        
        for label, material in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]
            
            a_m += w * ufl.inner(self.sigma_mech(ufl.TrialFunction(self.V_m), material), self.epsilon(v_m)) * dx
            
            if T_current is not None:
                L_m -= w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx

            rho = dolfinx.default_scalar_type(material.get("rho", 0.0))
            g = dolfinx.default_scalar_type(9.81)
            body_f = fem.Constant(self.mesh, (0.0, -rho*g) if self.regime != "3d" else (0.0, 0.0, -rho*g))
            L_m += w * ufl.dot(body_f, v_m) * dx

        for label in self.materials:
            for bc in self.traction[label]:
                L_m += w * ufl.dot(bc["value"], v_m) * self.ds_tags[bc["id"]]

        # Résolution
        petsc_opts = self.get_solver_options("mechanical", rtol=rtol_mech)
        problem = dolfinx.fem.petsc.LinearProblem(a_m, L_m, bcs=bcs_m, u=u_new, petsc_options=petsc_opts, petsc_options_prefix="mechanical_")
        problem.solve()

        # Relaxation et Convergence
        u_new.x.array[:] = self.relax_u * u_new.x.array + (1.0 - self.relax_u) * u_old.x.array
        dolfinx.fem.set_bc(u_new.x.array, bcs_m)
        u_new.x.scatter_forward()
        u_old.x.scatter_forward()
        
        diff = u_new.x.petsc_vec.copy()
        diff.axpy(-1.0, u_old.x.petsc_vec)
        norm_du = diff.norm(PETSc.NormType.NORM_2)
        norm_u = u_new.x.petsc_vec.norm(PETSc.NormType.NORM_2)
        rel_norm_du = norm_du / norm_u if norm_u > 1e-15 else norm_du
        
        res_curr = rel_norm_du if self.mech_cfg.get("convergence") == "rel_norm" else norm_du
        if self.relax_adaptive and prev_res_u is not None:
            self.relax_u = min(self.relax_u * self.relax_growth, self.relax_max) if res_curr < prev_res_u else max(self.relax_u * self.relax_shrink, self.relax_min)
        
        return res_curr < stag_tol_mech, norm_du, rel_norm_du, res_curr

    def _damage_step(self, D_new, D_old, rtol_dmg, stag_tol_dmg, prev_res_D):
        D_old.x.array[:] = D_new.x.array
        w = self.weight
        lc = float(self.dmg_cfg["lc"])
        damage_type = self.dmg_cfg["type"]
        print(f"\n[INFO] Assembling damage ({damage_type}) problem...")

        u_d, v_d = ufl.TrialFunction(self.V_d), ufl.TestFunction(self.V_d)
        a_d, L_d = 0, 0
        bcs_d = [bc["value"] for mat in self.materials for bc in self.dirichlet_damage.get(mat, [])]

        for label, material in self.materials.items():
            dx = self.dx_tags[self.label_map[label]]
            Gc = material["Gc"]
            if damage_type == "AT2":
                a_d += w*((self.H + 1.0) * u_d * v_d + lc**2 * ufl.inner(ufl.grad(u_d), ufl.grad(v_d))) * dx
                L_d += w*self.H * v_d * dx
            elif damage_type == "AT1":
                cw = 8.0/3.0; pref = Gc/cw; grad_c = 2.0*pref*lc; d_shift = 1e-8*(Gc/lc)
                a_d += w*(2.0*self.H + d_shift)*u_d*v_d*dx + w*grad_c*ufl.inner(ufl.grad(u_d), ufl.grad(v_d))*dx
                L_d += w*(2.0*self.H - (pref/lc))*v_d*dx

        opts = self.get_solver_options("damage", solver_type=self.dmg_cfg["linear_solver"], rtol=rtol_dmg)
        dolfinx.fem.petsc.LinearProblem(a_d, L_d, bcs=bcs_d, u=D_new, petsc_options=opts, petsc_options_prefix="damage_").solve()
        
        D_new.x.array[:] = self.relax_D * D_new.x.array + (1 - self.relax_D) * D_old.x.array
        D_new.x.array[:] = np.clip(np.maximum(D_new.x.array, D_old.x.array), 0.0, 1.0)
        D_new.x.scatter_forward(); D_old.x.scatter_forward()
        
        diff = D_new.x.petsc_vec.copy(); diff.axpy(-1.0, D_old.x.petsc_vec)
        ndD = diff.norm(PETSc.NormType.NORM_2); nD = D_new.x.petsc_vec.norm(PETSc.NormType.NORM_2)
        rndD = ndD / nD if nD > 1e-15 else ndD
        res = rndD if self.dmg_cfg.get("convergence") == "rel_norm" else ndD
        if self.relax_adaptive and prev_res_D is not None:
            self.relax_D = min(self.relax_D * self.relax_growth, self.relax_max) if res < prev_res_D else max(self.relax_D * self.relax_shrink, self.relax_min)
        return res < stag_tol_dmg, ndD, rndD, res

    def _cluster_step(self, c_new, c_old, dt):
        c_old.x.array[:] = c_new.x.array
        u_c, v_c = ufl.TrialFunction(self.V_c), ufl.TestFunction(self.V_c)
        v_v, D_v, dt_v = fem.Constant(self.mesh, self.v_cluster), fem.Constant(self.mesh, self.D_cluster), fem.Constant(self.mesh, dt)
        n, h = ufl.FacetNormal(self.mesh), ufl.CellDiameter(self.mesh)
        a = (u_c/dt_v)*v_c*ufl.dx; L = (c_old/dt_v)*v_c*ufl.dx
        if dt > 0:
            a += -u_c*v_v*v_c.dx(0)*ufl.dx + (ufl.avg(u_c*v_v*n[0])*ufl.jump(v_c) + 0.5*abs(v_v*n[0]('+'))*ufl.jump(u_c)*ufl.jump(v_c))*ufl.dS
            a += ufl.conditional(v_v*n[0]>0, v_v*n[0]*u_c*v_c, 0.0)*ufl.ds + D_v*u_c.dx(0)*v_c.dx(0)*ufl.dx
            a += -D_v*ufl.avg(u_c.dx(0))*ufl.jump(v_c, n[0])*ufl.dS - D_v*ufl.avg(v_c.dx(0))*ufl.jump(u_c, n[0])*ufl.dS
            a += D_v*(10.0/((h('+')+h('-'))/2.0))*ufl.jump(u_c)*ufl.jump(v_c)*ufl.dS
            dolfinx.fem.petsc.LinearProblem(a, L, u=c_new, petsc_options={"ksp_type":"gmres","pc_type":"ilu"}).solve()
        else: c_new.x.array[:] = c_old.x.array
        ct = self.mesh.comm.allreduce(fem.assemble_scalar(fem.form(c_new*ufl.SpatialCoordinate(self.mesh)[0]*ufl.dx)), MPI.SUM)
        if getattr(self, "C_tot_target", None) and abs(ct)>0: c_new.x.array[:] *= self.C_tot_target/ct

    def solve_staggered(self, max_iter=20, dt=0.0, stag_tol_th=1e-3, stag_tol_mech=1e-3, stag_tol_dmg=1e-3, rtol_th=1e-6, rtol_mech=1e-6, rtol_dmg=1e-5):
        self.dt = dt
        self._build_measures()

        # Initialize variables to avoid UnboundLocalError when models are disabled
        T_n = None
        u_n = getattr(self, "u", None)

        if self.on.get("thermal"):
            T_n = fem.Function(self.V_t); T_n.x.array[:] = self.T.x.array; T_o = fem.Function(self.V_t)
            b_t = [bc for l in self.dirichlet_thermal.values() for bc in l]
        if self.on.get("mechanical"):
            u_n = fem.Function(self.V_m); u_n.x.array[:] = self.u.x.array; u_o = fem.Function(self.V_m)
            b_m = [bc for l in self.dirichlet_mechanical.values() for bc in l]
        if self.on.get("damage"): D_n = fem.Function(self.V_d); D_n.x.array[:] = self.D.x.array; D_o = fem.Function(self.V_d)
        if self.on.get("cluster"): c_n = fem.Function(self.V_c); c_n.x.array[:] = self.c.x.array; self.c_n.x.array[:] = self.c.x.array
        rT = ru = rD = None
        for i in range(max_iter):
            print(f"\n--- Staggering iteration {i+1}/{max_iter} ---")
            cth = cm = cd = True
            if self.on.get("thermal"): cth, _, _, rT = self._thermal_step(T_n, T_o, b_t, rtol_th, stag_tol_th, rT)
            if self.on.get("mechanical"): cm, _, _, ru = self._mechanical_step(u_n, u_o, b_m, rtol_mech, stag_tol_mech, ru, T_current=T_n)
            if self.on.get("damage"): self.update_history(u_n, T=T_n); cd, _, _, rD = self._damage_step(D_n, D_o, rtol_dmg, stag_tol_dmg, rD)
            if self.on.get("cluster"): self._cluster_step(c_n, self.c_n, dt)
            if cth and cm and cd:
                print(f"[SUCCESS] Converged in {i+1} iterations.")
                if self.on.get("thermal"): self.T.x.array[:] = T_n.x.array
                if self.on.get("mechanical"): self.u.x.array[:] = u_n.x.array
                if self.on.get("damage"): self.D.x.array[:] = D_n.x.array
                if self.on.get("cluster"): self.c.x.array[:] = c_n.x.array
                if self.on.get("mechanical") and self.on.get("plasticity"): self.update_plastic_history(u_n)
                return True
        print("[WARNING] Did not converge."); return False


class MechanicalSolver:
    def __init__(self, config):
        self.comm = MPI.COMM_WORLD
        self.config = config

        # 1. Charger les deux maillages
        self.meshes = self._load_meshes()

        # 2. Initialiser les espaces fonctionnels
        self.V_m = {
            "steel_in": fem.VectorFunctionSpace(self.meshes["steel_in"], ("CG", 1)),
            "steel_o": fem.VectorFunctionSpace(self.meshes["steel_o"], ("CG", 1))
        }

        # 3. Charger les matériaux
        self.materials = self._load_materials(config["materials"])

        # 4. Charger les conditions aux limites
        self.dirichlet_mechanical, self.traction = self._load_boundary_conditions(config["boundary_conditions_path"])

        # 5. Initialiser les champs de déplacement
        self.u_old = {name: fem.Function(self.V_m[name]) for name in self.meshes}
        self.u_new = {name: fem.Function(self.V_m[name]) for name in self.meshes}

        # 6. Définir les surfaces de contact
        self.contact_surfaces = self._find_contact_surfaces()

        # 7. Initialiser le problème de contact
        self.contact = self._init_contact()

        # 8. Autres initialisations
        self.current_step = 0
        self.relax_u = config["solver_settings"]["relax_u"]

    def _load_meshes(self):
        # Charger les deux maillages depuis les fichiers .msh
        meshes = {}
        for name in ["steel_in", "steel_o"]:
            with dolfinx.io.XDMFFile(self.comm, f"{name}.xdmf", "r") as xdmf:
                meshes[name] = xdmf.read_mesh()
        return meshes

    def _find_contact_surfaces(self):
        # Trouver les facettes de contact pour chaque maillage
        contact_surfaces = {}
        for name, mesh_data in self.meshes.items():
            if name == "steel_in":
                # Surface de contact côté steel_in : steel_in_right (tag 4)
                contact_surfaces[name] = mesh.find_facets(
                    mesh_data,
                    lambda x: np.isclose(x[0], self.config["geometry"]["R_mid"] - self.config["geometry"]["gap"]/2)
                )
            elif name == "steel_o":
                # Surface de contact côté steel_o : steel_o_left (tag 8)
                contact_surfaces[name] = mesh.find_facets(
                    mesh_data,
                    lambda x: np.isclose(x[0], self.config["geometry"]["R_mid"] + self.config["geometry"]["gap"]/2)
                )
        return contact_surfaces

    def _init_contact(self):
        # Initialiser le problème de contact
        mesh_in = self.meshes["steel_in"]
        mesh_o = self.meshes["steel_o"]
        facets_in = self.contact_surfaces["steel_in"]
        facets_o = self.contact_surfaces["steel_o"]

        return PenaltyContact(
            slave=(mesh_in, facets_in),
            master=(mesh_o, facets_o),
            penalty=1e8,  # Coefficient de pénalité pour l'acier
            gap=0.0       # Jeu initial (0 pour contact parfait)
        )

    def _load_materials(self, materials_path):
        materials = {}
        for name, path in materials_path.items():
            with open(path, "r") as f:
                materials[name] = yaml.safe_load(f)
        return materials

    def _load_boundary_conditions(self, bc_path):
        with open(bc_path, "r") as f:
            bc_data = yaml.safe_load(f)

        dirichlet_mechanical = {}
        traction = {}

        for domain, bcs in bc_data["mechanical"].items():
            dirichlet_mechanical[domain] = []
            traction[domain] = []

            for bc in bcs:
                if bc["type"] == "Clamp_x":
                    # Condition de Dirichlet (encastrement)
                    mesh = self.meshes[domain]
                    dofs = fem.locate_dofs_geometrical(
                        self.V_m[domain],
                        lambda x: self._get_region_function(bc["region"])(x)
                    )
                    bc_value = fem.Constant(mesh, [bc["value"], 0.0, 0.0])  # Clamp en x seulement
                    dirichlet_bc = fem.DirichletBC(bc_value, dofs)
                    dirichlet_mechanical[domain].append(dirichlet_bc)

                elif bc["type"] == "Neumann":
                    # Condition de Neumann (traction)
                    traction[domain].append({
                        "region": bc["region"],
                        "value": bc["traction"],
                        "const": fem.Constant(self.meshes[domain], bc["traction"])
                    })

        return dirichlet_mechanical, traction

    def _get_region_function(self, region_name):
        # Retourne une fonction pour identifier une région (ex : steel_in_left)
        if region_name == "steel_in_left":
            return lambda x: np.isclose(x[0], self.config["geometry"]["R_int"])
        elif region_name == "steel_o_right":
            return lambda x: np.isclose(x[0], self.config["geometry"]["R_o"])
        elif region_name == "steel_in_right":
            return lambda x: np.isclose(x[0], self.config["geometry"]["R_mid"] - self.config["geometry"]["gap"]/2)
        elif region_name == "steel_o_left":
            return lambda x: np.isclose(x[0], self.config["geometry"]["R_mid"] + self.config["geometry"]["gap"]/2)
        else:
            raise ValueError(f"Unknown region: {region_name}")
