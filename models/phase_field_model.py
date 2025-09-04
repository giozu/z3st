# models/phase_field_model.py
import numpy as np
import ufl
import dolfinx
from petsc4py import PETSc
import yaml

class PhaseFieldModel:
    """
    Phase-field fracture (AT-2 type formulation) based on the paper,
    with volumetric-deviatoric split (Amor et al.) and
    irreversibility enforced via a history field H(x) = max_tau Psi_e^+(tau).

    Energy functional (paper's (2.2) with H+ and g(s)=s^2+eta):
        E_l(u, s) = ∫_Omega (g(s) H+ + Psi_el.0^-(epsilon(u))) dOmega + Gc ∫_Omega ((1-s)^2/(4l) + l |∇s|^2) dOmega
    
    The weak form for the phase-field problem (s-problem) is derived by taking the
    variational derivative of the energy functional (2.2) with respect to s,
    using g(s) = s^2 + eta:
        ∫_Omega ( (2s * H+ - Gc/(2l) * (1-s)) * q + Gc * l * ∇s ⋅ ∇q ) dOmega = 0
    where q is the test function for s.

    PUBLIC API used by the solver:
      - setup_phasefield_space()
      - update_history(u)           # updates H from current u
      - build_damage_forms()        # builds the non-linear residual F_s
      - assemble_solve_damage()     # solves for s with NewtonSolver
      - sigma_degraded(u_trial, s_field, material) # degraded stress tensor (paper's (2.7))
      - _material_constants(E, nu) # Helper to get kappa, mu from E, nu
    """

    # -------- init & parameters --------
    def __init__(self, input_file="input.yaml"):
        print("__PhaseFieldModel initializer__")

        with open(input_file, "r") as f:
            data = yaml.safe_load(f)

        pf = data.get("fracture", {})
        if not pf:
            raise ValueError("[PhaseFieldModel] 'fracture' missing in input.yaml.")

        self.Gc = pf.get("Gc", 0.01) # N/mm
        self.ell = pf.get("ell", 0.6241) # mm
        self.eta = pf.get("eta", 1e-8) # Small parameter for g(s) to avoid division by zero

        print("\n[PhaseFieldModel] parameters loaded from input.yaml:")
        print(f"  Gc               : {self.Gc} N/mm")
        print(f"  ell              : {self.ell} mm")
        print(f"  eta              : {self.eta} (small regularization parameter)\n")

        # FE containers
        self.dirichlet_s = []     # optional Dirichlet BCs on s
        self.F_s_residual = None  # Non-linear residual for s-problem

    # -------- helpers --------
    def _I(self): return ufl.Identity(self.gdim)

    def _pos(self, x):  # <x>_+ = max(x, 0)
        return ufl.max_value(x, 0.0)

    def _neg(self, x):  # <x>_- = min(x, 0)
        return ufl.min_value(x, 0.0)

    def _dev(self, eps):
        return eps - (ufl.tr(eps)/self.gdim) * self._I()

    # def g_s(self, s_val):
    #     """Degradation function for phase-field, g(s) = s^2 + eta (paper Eq. 2.5)."""
    #     return s_val**2 + self.eta

    def g_s(self, s_val):
        return s_val**2 + self.eta

    def g_prime(self, s_val):
        """Derivative of degradation function, g'(s) = 2s."""
        return 2 * s_val

    def setup(self):
        self.s   = dolfinx.fem.Function(self.V_s, name="PhaseField")
        self.s.x.array[:] = 0.0

        self.s_n = dolfinx.fem.Function(self.V_s, name="PhaseField_n") # For irreversibility
        self.s_n.x.array[:] = 0.0

        self.H   = dolfinx.fem.Function(self.V_s, name="History")
        self.H.x.array[:] = 0.0

        self.tr_eps_pos_n = dolfinx.fem.Function(self.V_s, name="TrEpsPos_n")
        self.tr_eps_pos_n.x.array[:] = 0.0

        self.tr_eps_neg_n = dolfinx.fem.Function(self.V_s, name="TrEpsNeg_n")
        self.tr_eps_neg_n.x.array[:] = 0.0

    # -------- elastic split (Amor et al.) --------
    def psi_positive_elastic_energy(self, u, material):
        """
        Psi_el.0^+(epsilon(u)) = 0.5 * K * (tr_eps)_+^2 + G * ||eps_dev||^2 (paper Eq. 2.4)
        """
        eps = ufl.sym(ufl.grad(u))
        tr_eps = ufl.tr(eps)
        eps_dev = self._dev(eps)

        kappa = material["bulk_modulus"]
        G = material["G"]

        tr_eps_pos = self._pos(tr_eps)
        return 0.5 * kappa * tr_eps_pos**2 + G * ufl.inner(eps_dev, eps_dev)

    def psi_negative_elastic_energy(self, u, material):
        """
        Psi_el.0^-(epsilon(u)) = 0.5 * K * (tr_eps)_-^2 (paper Eq. 2.4)
        """
        eps = ufl.sym(ufl.grad(u))
        tr_eps = ufl.tr(eps)

        kappa = material["bulk_modulus"]

        tr_eps_neg = self._neg(tr_eps)
        return 0.5 * kappa * tr_eps_neg**2

    def sigma_degraded(self, u, s, material):
        g = self.g_s(s)
        sigma = self.sigma_mech(u, material)
        return g * sigma

    def update_history(self, u):
        """
        H(x) := max( H(x), Psi_e^+(u)(x) ) over all P1 DOFs.
        Assumes self.materials is set up (e.g., in main solver class).
        """
        V1 = self.V_s # P1 space for interpolation

        H_new = dolfinx.fem.Function(V1)
        H_new.x.array[:] = self.H.x.array # Start with current history field

        # If multiple materials are defined, update H for each region.
        # Otherwise, assume a single material covering the entire domain.
        if self.materials and self.volume_tags:
            for name, mat_params in self.materials.items():
                tag = self.label_map[name] # Assuming label_map is also set (e.g., in main)
                cells = self.volume_tags.find(tag)

                # Create a cell-wise indicator function
                chi0 = dolfinx.fem.Function(dolfinx.fem.functionspace(self.mesh, ("DG", 0)))
                chi0.x.array[:] = 0.0
                chi0.x.array[cells] = 1.0

                # Interpolate the indicator to P1 nodes
                chi1 = dolfinx.fem.Function(V1)
                chi1.interpolate(dolfinx.fem.Expression(chi0, V1.element.interpolation_points()))

                # Compute Psi_e^+ on P1 nodes within this material region
                psi_p = self.psi_positive_elastic_energy(u, mat_params)
                psi1  = dolfinx.fem.Function(V1)
                psi1.interpolate(dolfinx.fem.Expression(psi_p, V1.element.interpolation_points()))

                # Update H using the mask
                mask = (chi1.x.array > 0.5)
                H_arr   = H_new.x.array
                psi_arr = psi1.x.array
                np.maximum(H_arr[mask], psi_arr[mask], out=H_arr[mask])
        else: # Single material case
            # Assumes self.materials has one entry, or default material params are used directly
            # For this example, we need *some* material_params, let's use the first one if available
            # or default to a dummy one.
            if self.materials:
                mat_params = list(self.materials.values())[0] # Use first material if defined
            else:
                # Fallback: create dummy material params if not set up externally
                kappa, mu = self._material_constants(18.0e3, 0.2) # Example values from paper Table 2
                mat_params = {"kappa": kappa, "mu": mu}

            psi_p = self.psi_positive_elastic_energy(u, mat_params)
            psi1  = dolfinx.fem.Function(V1)
            psi1.interpolate(dolfinx.fem.Expression(psi_p, V1.element.interpolation_points()))
            np.maximum(self.H.x.array, psi1.x.array, out=self.H.x.array) # Update directly

        self.H.x.array[:] = H_new.x.array if (self.materials and self.volume_tags) else self.H.x.array
        # Ensure H is not NaN if any part of the domain has zero elastic energy (e.g. initial state)
        # This isn't strictly necessary with np.maximum, but good for robustness
        self.H.x.array[np.isnan(self.H.x.array)] = 0.0 

    # -------- damage forms (for non-linear solve) --------
    def build_damage_forms(self):
        """
        Build the non-linear residual form for the s-problem (phase-field).
        F_s = ∫_Omega ( (g'(s) H - Gc/(2l) * (1-s)) * q + Gc * l * ∇s ⋅ ∇q ) dOmega = 0
        """
        dx = ufl.Measure("dx", domain=self.mesh)

        # Residual from variational derivative of Energy Functional (2.2) wrt s
        # using g(s) = s^2 + eta, so g'(s) = 2s
        F_s = (self.g_prime(self.s) * self.H - (self.Gc / (2 * self.ell)) * (1 - self.s)) * self.q_test * dx \
              + self.Gc * self.ell * ufl.inner(ufl.grad(self.s), ufl.grad(self.q_test)) * dx
        
        self.F_s_residual = F_s
        print("[PF] Non-linear phase-field residual form built.")

    # -------- non-linear solve for s --------
    def assemble_solve_damage(self):
        """
        Solve the non-linear system for s using a NewtonSolver.
        Applies (optional) Dirichlet BCs in self.dirichlet_s and enforces irreversibility.
        """
        from dolfinx.fem.petsc import NonlinearProblem
        from dolfinx.nls.petsc import NewtonSolver

        problem_s = NonlinearProblem(self.F_s_residual, self.s, bcs=self.dirichlet_s)
        solver_s = NewtonSolver(self.mesh.comm, problem_s)
        solver_s.convergence_criterion = "incremental" # Use "incremental" for delta_u norm
        solver_s.rtol = 1e-6
        solver_s.atol = 1e-8
        solver_s.max_it = 20
        solver_s.report = False
        
        # solver_s.view() # Uncomment to see solver options

        num_iters, converged = solver_s.solve(self.s)
        if not converged:
            print(f"WARNING: s-solver did not converge in {num_iters} iterations.")
            # Depending on robustness needs, might raise an error or try a smaller step

        # Enforce physical constraints: 0 <= s <= 1
        np.clip(self.s.x.array, 0.0, 1.0, out=self.s.x.array)

        # Enforce irreversibility: s >= s_n
        if self.s_n is not None:
            np.maximum(self.s.x.array, self.s_n.x.array, out=self.s.x.array)
            self.s_n.x.array[:] = self.s.x.array # Update s_n for the next step

        return num_iters, converged