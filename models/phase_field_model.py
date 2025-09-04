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
    
    def g_s(self, dmg):
        """Degradation function for phase-field"""
        return dmg**2 + self.eta

    def g_prime(self, dmg):
        """Derivative of degradation function"""
        return 2 * dmg

    def setup(self):
        self.s   = dolfinx.fem.Function(self.V_pf, name="PhaseField")
        self.s.x.array[:] = 0.0

        self.s_n = dolfinx.fem.Function(self.V_pf, name="PhaseField_n") # For irreversibility
        self.s_n.x.array[:] = 0.0

        self.H   = dolfinx.fem.Function(self.V_pf, name="History")
        self.H.x.array[:] = 0.0

        self.tr_eps_pos_n = dolfinx.fem.Function(self.V_pf, name="TrEpsPos_n")
        self.tr_eps_pos_n.x.array[:] = 0.0

        self.tr_eps_neg_n = dolfinx.fem.Function(self.V_pf, name="TrEpsNeg_n")
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
    
    def sigma_pos(self, s, u, material):

        K = material["bulk_modulus"]
        G = material["G"]

        return self.g_s(s) * (2*G*self._dev(self.epsilon(u))) + K*self._pos(ufl.tr(self.epsilon(u)))*self._I

    def sigma_neg(self, s, u, material):

        K = material["bulk_modulus"]

        return K*self._neg(ufl.tr(self.epsilon(u)))*self._I

    def psi_split(self, u, material):
        tr_eps = ufl.tr(ufl.sym(ufl.grad(u)))
        eps_dev = self._dev(ufl.sym(ufl.grad(u)))
        
        kappa = material["bulk_modulus"]
        G = material["G"]

        tr_eps_pos = self._pos(tr_eps)
        tr_eps_neg = self._neg(tr_eps)

        psi_plus  = 0.5 * kappa * tr_eps_pos**2 + G * ufl.inner(eps_dev, eps_dev)
        psi_minus = 0.5 * kappa * tr_eps_neg**2
        return psi_plus, psi_minus


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