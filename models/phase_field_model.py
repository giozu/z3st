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
        print(f"  Gc               : {self.Gc} N/m")
        print(f"  ell              : {self.ell} m")
        print(f"  eta              : {self.eta} (small regularization parameter)\n")

        # FE containers
        self.dirichlet_s = []     # optional Dirichlet BCs on s
        self.F_s_residual = None  # Non-linear residual for s-problem

    # -------- helpers --------
    def _I(self): return ufl.Identity(self.gdim)

    def _pos(self, f):
        # f_+ = f if f>0
        return ufl.max_value(f, 0.0)

    def _neg(self, f):
        # f_- = -f if f<0
        return ufl.min_value(f, 0.0)

    def _dev(self, eps):
        return eps - (ufl.tr(eps)/self.gdim) * self._I()
    
    def g_s(self, dmg):
        """Degradation function for phase-field"""
        # return (1.0 - self.eta) * dmg**2 + self.eta
        return dmg**2 + self.eta

    def g_prime(self, dmg):
        """Derivative of degradation function"""
        # return  (1.0 - self.eta) * 2 * dmg
        return  2 * dmg

    def setup(self):
        
        self.d   = dolfinx.fem.Function(self.V_pf, name="PhaseField")
        self.d.x.array[:] = 0.0

        self.H   = dolfinx.fem.Function(self.V_pf, name="History")
        self.H.x.array[:] = 0.0

    def elastic_energy_pos(self, u, material):
        """
        Psi_el.0^+(epsilon(u)) = 0.5 * K * (tr_eps)_+^2 + G * ||eps_dev||^2
        """
        eps = ufl.sym(ufl.grad(u))
        tr_eps = ufl.tr(eps)
        eps_dev = self._dev(eps)

        kappa = material["bulk_modulus"]
        G = material["G"]

        tr_eps_pos = self._pos(tr_eps)
        return 0.5 * kappa * tr_eps_pos**2 + G * ufl.inner(eps_dev, eps_dev)

    def elastic_energy_neg(self, u, material):
        """
        Psi_el.0^-(epsilon(u)) = 0.5 * K * (tr_eps)_-^2
        """
        eps = ufl.sym(ufl.grad(u))
        tr_eps = ufl.tr(eps)

        kappa = material["bulk_modulus"]

        tr_eps_neg = self._neg(tr_eps)
        return 0.5 * kappa * tr_eps_neg**2

    def sigma_pos(self, u, s, material):
        K = material["bulk_modulus"]
        G = material["G"]
        eps = self.epsilon(u)
        return ( self.g_s(s) * ( 2*G*self._dev(eps) ) + K * self._pos(ufl.tr(eps)) * self._I() )

    def sigma_neg(self, u, material):
        K = material["bulk_modulus"]
        tr_neg = self._neg(ufl.tr(self.epsilon(u)))
        return K * tr_neg * self._I()

    # def sigma_degraded(self, u, s, material):
    #     return self.sigma_pos(u, s, material) + self.sigma_neg(u, material)
    
    def sigma_degraded(self, u, s, material):
        lmbda = material["lmbda"]
        mu = material["mu"]
        eps = self.epsilon(u)                # sym(grad(u))
        tr_eps = ufl.tr(eps)

        # split volumetric strain
        tr_pos = (tr_eps + ufl.abs(tr_eps)) / 2
        tr_neg = (tr_eps - ufl.abs(tr_eps)) / 2

        dev = eps - (1.0/3.0)*tr_eps*ufl.Identity(self.gdim)

        return (
            ( (1-s)**2 + self.eta ) * (2.0*mu*dev + lmbda*tr_pos*ufl.Identity(self.gdim))
            + lmbda*tr_neg*ufl.Identity(self.gdim)
        )
