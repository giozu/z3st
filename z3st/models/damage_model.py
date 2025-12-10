# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import numpy as np
import ufl


class DamageModel:
    def __init__(self):
        print("[DamageModel] initializer")

        data = self.input_file

        self.damage_options = data.get("damage", {})
        if not self.damage_options:
            raise ValueError("[DamageModel] 'damage' missing in input.yaml.")

        print("[DamageModel] options loaded from input.yaml:")
        for key, value in self.damage_options.items():
            print(f"  {key:<20}: {value}")

    @staticmethod
    def degradation_function(D, K=1e-4):
        """
        Stress degradation function g_d(D) = (1-D)**2 + K

        :param D: scalar field variable. D=0 corresponds to the undamaged state, D=1 to the fully damaged state.
        :param K: numerical regularization parameter
        """
        return (1 - D) ** 2 + K

    def effective_strain_energy(self, u, material):
        eps = ufl.sym(ufl.grad(u))
        lam = material["lmbda"]
        G = material["G"]
        E = material["E"]

        sigma_eff = lam * ufl.tr(eps) * ufl.Identity(3) + 2 * G * eps
        tr_pos = 0.5 * (ufl.tr(sigma_eff) + abs(ufl.tr(sigma_eff)))  # ⟨tr(σ)⟩
        psi_eff = 0.5 / E * tr_pos**2
        return psi_eff

    def crack_driving_force(self, u, material):
        """
        Crack-driving quantity H(u) in forma energetica positiva.
        Usa le costanti di un materiale (E, λ, G, σ_c, ζ).
        """
        E = material["E"]
        lam = material["lmbda"]
        G = material["G"]
        sigma_c = material["sigma_c"]
        zeta = material["zeta"]

        eps = self.epsilon(u)
        tr_eps = ufl.tr(eps)
        tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
        psi_pos = 0.5 * lam * tr_eps_pos**2 + G * ufl.inner(ufl.dev(eps), ufl.dev(eps))

        psi_c = sigma_c**2 / (2 * E)

        H = zeta * ufl.conditional(
            ufl.gt(psi_pos / psi_c, 1.0),
            psi_pos / psi_c - 1.0,
            0.0,
        )
        return H

    @staticmethod
    def gamma_density(D, grad_D, lc):
        """Fracture energy density γ(D, ∇D) = 0.5*(D**2/lc + lc*|∇D|²)."""
        return 0.5 * (D**2 / lc + lc * ufl.dot(grad_D, grad_D))

    def update_history(self, u):
        """
        Aggiorna il campo di storia H in modo irreversibile:
            H(x) = max(H_old(x), H_new(x))

        Proiezione L2 su DG0 tramite Expression + interpolate.
        Per semplicità usa i parametri del primo materiale.
        """
        if not hasattr(self, "H"):
            raise RuntimeError("[DamageModel] self.H non inizializzato.")

        mat_name, mat = next(iter(self.materials.items()))
        H_expr = self.crack_driving_force(u, mat)

        Q = self.Q
        H_new = dolfinx.fem.Function(Q, name="CrackDrivingForce_new")
        expr = dolfinx.fem.Expression(H_expr, Q.element.interpolation_points)
        H_new.interpolate(expr)

        H_old_arr = self.H.x.array
        H_new_arr = H_new.x.array
        self.H.x.array[:] = np.maximum(H_old_arr, H_new_arr)
