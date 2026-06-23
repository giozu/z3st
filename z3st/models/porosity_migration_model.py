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
        
        # Load physical model constants (Sens 1972)
        # Default parameters matching Barani et al. (2022) / Novascone (2018)
        self.c0 = float(self.porosity_cfg.get("c0", 5.0e-16))
        self.c1 = float(self.porosity_cfg.get("c1", 0.988))
        self.c2 = float(self.porosity_cfg.get("c2", 6.395e-6))
        self.c3 = float(self.porosity_cfg.get("c3", 3.543e-9))
        self.c4 = float(self.porosity_cfg.get("c4", 3.0e-12))
        self.P0_s = float(self.porosity_cfg.get("P0", 33708.5))
        self.H_s = float(self.porosity_cfg.get("Hs", 5.98e5))

        print("[PorosityMigrationModel] options loaded:")
        for key in ["c0", "c1", "c2", "c3", "c4", "P0", "Hs"]:
            val = self.porosity_cfg.get(key)
            if val is not None:
                print(f"  {key:<20}: {val}")

    def set_porosity_initial_conditions(self):
        """
        Sets the initial condition for the porosity field using material cards.
        """
        print("\n[PorosityMigration] Setting initial porosity conditions...")
        self.p.x.array[:] = 0.0

        for name, mat in self.materials.items():
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_p)
            p_init = float(mat.get("initial_porosity", 0.0))
            self.p.x.array[dofs] = p_init
            print(f"  [INFO] Subdomain '{name}': set initial porosity to {p_init:.4f} ({len(dofs)} DOFs)")

        self.p.x.scatter_forward()
        self.p_n.x.array[:] = self.p.x.array
        self.p_n.x.scatter_forward()

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

    def _porosity_step(self, p_new, p_old, dt, T_current, stag_tol_p, prev_res_p):
        """
        Solve the pore advection step with SU stabilization.
        """
        if dt <= 0.0:
            print("  Porosity step: dt <= 0.0 -> preserving porosity (converged: True)")
            return True, 0.0, 0.0, prev_res_p

        p_old.x.array[:] = p_new.x.array
        p_old.x.scatter_forward()


        u_p, v_test = ufl.TrialFunction(self.V_p), ufl.TestFunction(self.V_p)
        dt_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(dt))

        # Pore velocity model parameters
        R = 8.314
        c0 = self.c0
        c1 = self.c1
        c2 = self.c2
        c3 = self.c3
        c4 = self.c4
        P0 = self.P0_s
        Hs = self.H_s

        # Mobility function
        # Based on the Sens (1972) pore velocity model in cgs units converted to SI:
        # Pre-exponential factor in dyn/cm2: p0_cgs = exp(33.7085)
        # Activation energy in erg/mol: Hs_cgs = Hs * 1e7
        # Conversion factor from cgs to SI: 1e-4
        p0_cgs = ufl.exp(33.7085)
        Hs_cgs = Hs * 1.0e7
        poly = c1 + c2 * T_current + c3 * T_current**2 + c4 * T_current**3
        mobility = c0 * poly * (T_current**(-2.5)) * Hs_cgs * p0_cgs * ufl.exp(-Hs / (R * T_current)) * 1.0e-4
        v_vec = mobility * ufl.grad(T_current)

        # SU stabilization terms
        h_e = ufl.CellDiameter(self.mesh)
        v_mag = ufl.sqrt(ufl.dot(v_vec, v_vec) + 1e-20)
        K_term = (h_e / (2.0 * v_mag)) * ufl.dot(v_vec, ufl.grad(u_p)) * ufl.dot(v_vec, ufl.grad(v_test))

        # Outward normal vector
        n_vec = ufl.FacetNormal(self.mesh)
        
        # Outermost boundary tag and measure for inflow boundary
        outer_id = self.label_map["outer"]
        ds_outer = self.ds_tags[outer_id]
        
        # Prescribed inflow porosity at the outer boundary (where flow is inwards, v · n < 0)
        p_inflow = float(self.materials["fuel"].get("initial_porosity", 0.15))
        p_inflow_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(p_inflow))

        # Bilinear and linear forms using conservative weak formulation:
        # We integrate the divergence term div(v * p) * test by parts:
        #   \int \nabla \cdot (p v) \psi dx = - \int p v \cdot \nabla \psi dx + \int_{\partial \Omega} p (v · n) \psi ds
        # Since v · n < 0 on the outer boundary (inflow) and v · n = 0 on symmetry boundaries (top/bottom),
        # the boundary term only needs to be integrated over the outer boundary where p is set to p_inflow.
        # This term goes to the RHS (L_p) since it doesn't contain the trial function u_p.
        
        a_p = (u_p / dt_const) * v_test * ufl.dx
        a_p -= u_p * ufl.dot(v_vec, ufl.grad(v_test)) * ufl.dx
        a_p += K_term * ufl.dx

        L_p = (p_old / dt_const) * v_test * ufl.dx
        L_p -= p_inflow_const * ufl.dot(v_vec, n_vec) * v_test * ds_outer

        # PETSc options for porosity linear solver
        petsc_opts = self.get_solver_options(
            solver_type=self.porosity_cfg.get("linear_solver", "direct_mumps"),
            physics="porosity",
            rtol=self.porosity_cfg.get("rtol", 1.0e-8),
        )

        problem_p = dolfinx.fem.petsc.LinearProblem(
            a_p, L_p, bcs=[], u=p_new,
            petsc_options=petsc_opts,
            petsc_options_prefix="porosity_"
        )
        problem_p.solve()

        # Clamp porosity to physical bounds: [0.0, 1.0]
        p_new.x.array[:] = np.clip(p_new.x.array, 0.0, 1.0)
        p_new.x.scatter_forward()

        # Under-relaxation
        relax_p = float(self.porosity_cfg.get("relax", 1.0))
        if relax_p < 1.0:
            p_new.x.array[:] = relax_p * p_new.x.array + (1.0 - relax_p) * p_old.x.array
            p_new.x.scatter_forward()

        # Convergence check: mixed relative/absolute criterion (Eq 2 in the paper)
        eps_rel = float(self.porosity_cfg.get("stag_tol_rel", 1.0e-6))
        eps_abs = float(self.porosity_cfg.get("stag_tol_abs", 1.0e-8))

        diff = np.abs(p_new.x.array - p_old.x.array)
        check_vals = diff - np.abs(p_new.x.array) * eps_rel - eps_abs
        max_val = np.max(check_vals) if check_vals.size > 0 else -1.0

        conv_p = max_val < 0.0

        print(f"  Porosity step: max(diff - rel - abs) = {max_val:.3e} (converged: {conv_p})")
        print(f"  p_new: min={p_new.x.array.min():.4f}, max={p_new.x.array.max():.4f}, mean={p_new.x.array.mean():.4f}")

        return conv_p, 0.0, 0.0, prev_res_p
