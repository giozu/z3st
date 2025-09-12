# models/phase_field_model.py
import numpy as np
import ufl
import dolfinx
import yaml

class PhaseFieldModel:
    """
    Phase-field fracture

    Convention:
        d = 0 → intact material
        d = 1 → fully broken
    """

    def __init__(self, input_file="input.yaml"):
        print("__PhaseFieldModel initializer__")

        with open(input_file, "r") as f:
            data = yaml.safe_load(f)

        pf = data.get("fracture", {})
        if not pf:
            raise ValueError("[PhaseFieldModel] 'fracture' missing in input.yaml.")

        self.Gc = pf.get("Gc", 0.01)   # Fracture toughness (N/m)
        self.ell = pf.get("ell", 0.001) # Regularization length (m)
        self.eta = pf.get("eta", 1e-8)  # Stability parameter

        print("\n[PhaseFieldModel] parameters loaded from input.yaml:")
        print(f"  Gc   : {self.Gc}")
        print(f"  ell  : {self.ell}")
        print(f"  eta  : {self.eta}\n")

    # --- helpers ---
    def _I(self): 
        return ufl.Identity(self.gdim)

    def _pos(self, f):  # positive part
        return ufl.max_value(f, 0.0)

    def _neg(self, f):  # negative part
        return ufl.min_value(f, 0.0)

    def _dev(self, eps):
        return eps - (ufl.tr(eps)/self.gdim) * self._I()

    def g_d(self, d):
        """Degradation function"""
        return (1.0 - d)**2 + self.eta

    def g_d_prime(self, d):
        """dg(d)/dd"""
        return -2.0*(1.0-d)

    def setup(self):
        """Initialize PF variables"""
        self.d = dolfinx.fem.Function(self.V_pf, name="PhaseField")
        self.d.x.array[:] = 0.0
        # self.H = dolfinx.fem.Function(self.V_pf, name="History")
        # self.H.x.array[:] = 0.0

    def psi_pos(self, u, material):
        """Positive elastic energy (tensile)"""
        eps = self.epsilon(u)
        tr_eps = ufl.tr(eps)
        tr_pos = self._pos(tr_eps)
        dev = self._dev(eps)
        K = material["bulk_modulus"]
        G = material["G"]
        return 0.5*K*tr_pos**2 + G*ufl.inner(dev, dev)

    def psi_neg(self, u, material):
        """Negative elastic energy (compressive)"""
        eps = self.epsilon(u)
        tr_eps = ufl.tr(eps)
        tr_neg = self._neg(tr_eps)
        K = material["bulk_modulus"]
        return 0.5*K*tr_neg**2

    def sigma_degraded(self, u, d, material):
        eps = ufl.variable(self.epsilon(u))
        psi_pos = self.psi_pos(u, material)
        psi_neg = self.psi_neg(u, material)
        sigma = self.g_d(d) * ufl.diff(psi_pos, eps) + ufl.diff(psi_neg, eps)
        return sigma
    
    def sigma_degraded(self, u, d, material):
        eps = ufl.variable(self.epsilon(u))

        tr_eps = ufl.tr(eps)
        tr_pos = self._pos(tr_eps)
        tr_neg = self._neg(tr_eps)
        dev = self._dev(eps)

        K = material["bulk_modulus"]
        G = material["G"]

        psi_pos = 0.5*K*tr_pos**2 + G*ufl.inner(dev, dev)
        psi_neg = 0.5*K*tr_neg**2

        sigma = self.g_d(d) * ufl.diff(psi_pos, eps) + ufl.diff(psi_neg, eps)
        return sigma
