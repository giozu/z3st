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

        # --. Damage model options --..
        self.dmg_cfg = self.input_file.get("damage", {})

        if not self.dmg_cfg:
            raise ValueError("[DamageModel] 'damage' missing in input.yaml.")

        print("[DamageModel] options loaded from input.yaml:")
        for key, value in self.dmg_cfg.items():
            print(f"  {key:<20}: {value}")

    @staticmethod
    def degradation_function(D, K=1e-4):
        """
        Stress degradation function g_d(D) = (1-D)**2 + K

        :param D: scalar field variable. D=0 corresponds to the undamaged state, D=1 to the fully damaged state.
        :param K: numerical regularization parameter
        """
        return (1 - D) ** 2 + K

    def crack_driving_force(self, u, material):
        """
        Crack-driving force H(u)
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
        Ensure the irreversibility for the field H:
            H(x) = max(H_old(x), H_new(x))
        """

        Q = self.Q
        H_new = dolfinx.fem.Function(Q)

        # iterate over materials
        for label, material in self.materials.items():
            tag = self.label_map[label]

            # get cells belonging to this material
            cells = np.where(self.cell_tags.values == tag)[0]
            if len(cells) == 0:
                continue

            # compute crack driving force for this material
            H_expr = self.crack_driving_force(u, material)
            tmp = dolfinx.fem.Function(Q)
            expr = dolfinx.fem.Expression(H_expr, Q.element.interpolation_points)
            tmp.interpolate(expr)

            # assign tmp values to H_new only on these cells
            for cell in cells:
                dof = Q.dofmap.cell_dofs(cell)[0]  # DG0 --> 1 dof per cell
                H_new.x.array[dof] = tmp.x.array[dof]

        # irreversibility
        self.H.x.array[:] = np.maximum(self.H.x.array, H_new.x.array)
