# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class PorosityMigrationModel:
    """
    Porosity migration model representing thermal-gradient-driven pore transport.
    """

    def __init__(self):
        print("[PorosityMigrationModel] initializer")

        self.porosity_cfg = self.input_file.get("porosity", {})

        # Sens (1972) pore-velocity correlation, J. Nucl. Mater. 43, 293, Eqs. 20-24,
        # as used by Novascone et al. (2018) Eq. 4 and Barani et al. (2022) Eq. 9
        # (the benchmark reproduced here), in SI units:
        #   v = v0 (c1 + c2 T + c3 T^2 + c4 T^3) T^-2.5 exp(-Hs/RT) grad T   [m/s]
        # - (c1..c4) are the UO2 molecular-volume polynomial from the Christensen
        #   density data, Sens Eq. 24: rho_T = 10.97/(c1 + c2 T + c3 T^2 + c4 T^3).
        # - Hs and the vapour-pressure pre-exponential ln(p0) = 33.7085 [dyn/cm^2]
        #   are Sens Table 1 (Ohse data, 1500-2000 K range): the Arrhenius term
        #   exp(-Hs/RT) is the UO2 vapour pressure p = p0 exp(-Hs/RT), Sens Eq. 18.
        # - v0 is the single lumped prefactor (Sens Eq. 21 geometric/kinetic factor
        #   4 p0 Omega_T DeltaH ... / (3 pi sigma^2 P N ...), which Novascone calls
        #   "the product of many other constants"). v0 = 1.303427e8 reproduces
        #   Sens' own Table 2 disc-pore velocity (4.2e-11 m/s at T=2000 K,
        #   gradT=2500 K/cm) to ~8%, and yields the published central-void radius
        #   (~0.2 r/Ro).
        self.v0 = float(self.porosity_cfg.get("v0", 1.303427e8))
        self.c1 = float(self.porosity_cfg.get("c1", 0.988))
        self.c2 = float(self.porosity_cfg.get("c2", 6.395e-6))
        self.c3 = float(self.porosity_cfg.get("c3", 3.543e-9))
        self.c4 = float(self.porosity_cfg.get("c4", 3.0e-12))
        self.H_s = float(self.porosity_cfg.get("Hs", 5.98e5))   # heat of vaporisation [J/mol] (Sens Table 1)

        print("[PorosityMigrationModel] options loaded:")
        for key, val in [
            ("v0", self.v0), ("c1", self.c1), ("c2", self.c2),
            ("c3", self.c3), ("c4", self.c4), ("Hs", self.H_s),
        ]:
            print(f"  {key:<20}: {val}")

    def set_porosity_initial_conditions(self):
        """
        Sets the initial condition for the porosity field using material cards.
        """
        print("\n[PorosityMigration] Setting initial porosity conditions...")
        self.porosity.x.array[:] = 0.0

        for name, mat in self.materials.items():
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_p)
            p_init = float(mat.get("initial_porosity", 0.0))
            self.porosity.x.array[dofs] = p_init
            print(f"  [INFO] Subdomain '{name}': set initial porosity to {p_init:.4f} ({len(dofs)} DOFs)")

        self.porosity.x.scatter_forward()
        self.porosity_n.x.array[:] = self.porosity.x.array
        self.porosity_n.x.scatter_forward()

    def update_porosity_dependent_properties(self, T_eval, p_eval):
        """
        Dynamically evaluate and update porosity-dependent material properties.
        This modifies mat["k"] in place so the staggered conduction step sees the
        updated field.
        """
        for name, mat in self.materials.items():
            if mat.get("thermal_conductivity_model") == "kato_porosity":
                T_vals = T_eval.x.array
                p_vals = p_eval.x.array

                x = float(mat.get("stoichiometry_deviation", 0.025))
                k_He = float(mat.get("helium_conductivity", 0.69))

                # Kato temperature dependent conductivity k(T)
                A = 2.713 * x + 1.595e-2
                B = (2.493 - 2.625 * x) * 1.0e-4
                C = 1.541e11
                D = 1.522e4

                T_clamped = np.clip(T_vals, 298.0, 4000.0)
                k_T = 1.0 / (A + B * T_clamped) + (C / (T_clamped**2.5)) * np.exp(-D / T_clamped)

                # Maxwell-Eucken correction for porosity p
                p_clamped = np.clip(p_vals, 0.0, 1.0)
                num = k_He + 2.0 * k_T - 2.0 * p_clamped * (k_T - k_He)
                den = k_He + 2.0 * k_T + p_clamped * (k_T - k_He)

                k_Tp = k_T * (num / den)

                mat["k"].x.array[:] = k_Tp
                mat["k"].x.scatter_forward()

    def _porosity_step(self, p_new, p_n, dt, T_current, stag_tol_p, prev_res_p):
        """
        Solve one pore-advection sub-step (Barani et al. 2022, Eq. 1b) with SU
        stabilisation, inside the staggered T-p fixed-point loop.

        p_new : current porosity iterate, updated in place.
        p_n   : porosity at the start of the time step (t^n). Held FIXED here and
                used only for the backward-Euler time term — it must NOT be
                overwritten, otherwise each staggered iteration would advect a
                full dt forward from the previous iterate.
        """
        if dt <= 0.0:
            print("  Porosity step: dt <= 0.0 -> preserving porosity (converged: True)")
            return True, 0.0, 0.0, prev_res_p

        # Previous staggered iterate: reference for the Eq. (2) convergence check
        # and for under-relaxation. Distinct from p_n (the t^n time level).
        p_prev_iter = p_new.x.array.copy()

        u_p, v_test = ufl.TrialFunction(self.V_p), ufl.TestFunction(self.V_p)
        dt_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(dt))

        R = 8.314  # universal gas constant [J/(mol K)]
        c1, c2, c3, c4 = self.c1, self.c2, self.c3, self.c4
        Hs = self.H_s

        # Pore velocity (Novascone 2018 Eq. 4 / Barani 2022 Eq. 9), SI. The lumped
        # prefactor v0 bundles Sens' c0, the heat-of-vaporisation prefactor and the
        # vapour-pressure pre-exponential; the Arrhenius term uses Hs [J/mol] and R
        # [SI]. v points up the thermal gradient (inward, toward the hot centre).
        poly = c1 + c2 * T_current + c3 * T_current**2 + c4 * T_current**3
        mobility = self.v0 * poly * (T_current**(-2.5)) * ufl.exp(-Hs / (R * T_current))
        v_vec = mobility * ufl.grad(T_current)

        # SU stabilisation: streamline diffusion K = (h_e / 2|v|) v⊗v (Barani Eq. 3-4),
        # contributing ∫ (h_e/2|v|) (v·∇p)(v·∇w) dx to the bilinear form.
        h_e = ufl.CellDiameter(self.mesh)
        v_mag = ufl.sqrt(ufl.dot(v_vec, v_vec) + 1e-30)
        K_term = (h_e / (2.0 * v_mag)) * ufl.dot(v_vec, ufl.grad(u_p)) * ufl.dot(v_vec, ufl.grad(v_test))

        # Conservative advection integrated by parts:
        #   ∫ ∇·(v p) w dx = -∫ p v·∇w dx + ∫_∂Ω p (v·n) w ds
        # Default: drop the boundary term (homogeneous natural BC, as in Barani
        # Sec. 4.1 for SU). v·n = 0 on the symmetry faces and is negligible at the
        # cold outer rim, so the rim keeps its fabrication porosity without an
        # imposed value. Set porosity.rim_inflow_porosity to force a rim value.
        a_p = (u_p / dt_const) * v_test * ufl.dx
        a_p -= u_p * ufl.dot(v_vec, ufl.grad(v_test)) * ufl.dx
        a_p += K_term * ufl.dx

        L_p = (p_n / dt_const) * v_test * ufl.dx

        rim_inflow = self.porosity_cfg.get("rim_inflow_porosity", None)
        if rim_inflow is not None:
            n_vec = ufl.FacetNormal(self.mesh)
            ds_outer = self.ds_tags[self.label_map["outer"]]
            p_inflow_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(float(rim_inflow)))
            L_p -= p_inflow_const * ufl.dot(v_vec, n_vec) * v_test * ds_outer

        petsc_opts = self.get_solver_options(
            solver_type=self.porosity_cfg.get("linear_solver", "direct_mumps"),
            physics="porosity",
            rtol=self.porosity_cfg.get("rtol", 1.0e-8),
        )
        problem_p = dolfinx.fem.petsc.LinearProblem(
            a_p, L_p, bcs=[], u=p_new,
            petsc_options=petsc_opts,
            petsc_options_prefix="porosity_",
        )
        problem_p.solve()

        # Clamp to the physical bounds [0, 1].
        p_new.x.array[:] = np.clip(p_new.x.array, 0.0, 1.0)
        p_new.x.scatter_forward()

        # Under-relaxation against the previous staggered iterate.
        relax_p = float(self.porosity_cfg.get("relax", 1.0))
        if relax_p < 1.0:
            p_new.x.array[:] = relax_p * p_new.x.array + (1.0 - relax_p) * p_prev_iter
            p_new.x.scatter_forward()

        # Convergence: mixed relative/absolute criterion between successive
        # staggered iterates (Barani Eq. 2), reduced over owned DOFs across ranks.
        eps_rel = float(self.porosity_cfg.get("stag_tol_rel", 1.0e-6))
        eps_abs = float(self.porosity_cfg.get("stag_tol_abs", 1.0e-8))
        n_owned = self.V_p.dofmap.index_map.size_local
        p_now = p_new.x.array[:n_owned]
        diff = np.abs(p_now - p_prev_iter[:n_owned])
        check_vals = diff - np.abs(p_now) * eps_rel - eps_abs
        local_max = float(np.max(check_vals)) if check_vals.size > 0 else -1.0
        max_val = self.mesh.comm.allreduce(local_max, op=MPI.MAX)
        conv_p = max_val < 0.0

        print(f"  Porosity step: max(diff - rel - abs) = {max_val:.3e} (converged: {conv_p})")
        print(f"  p_new: min={p_new.x.array.min():.4f}, max={p_new.x.array.max():.4f}, mean={p_new.x.array.mean():.4f}")

        return conv_p, 0.0, 0.0, prev_res_p
