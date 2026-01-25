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
        print("DamageModel initializer")

        # --. Damage model options --..
        self.dmg_cfg = self.input_file.get("damage", {})

        if not self.dmg_cfg:
            raise ValueError("'damage' entry missing in input.yaml.")

        print("Options loaded from input.yaml:")
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
        Crack-driving force
        """

        # Material properties
        E = material["E"]
        lam = material["lmbda"]
        G = material["G"]
        sigma_c = material["sigma_c"]
        zeta = material["zeta"]

        # Definitions
        eps = self.epsilon(u)
        tr_eps = ufl.tr(eps)
        tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
        # No cracking under pure hydrostatic compression
        psi_pos = 0.5 * lam * tr_eps_pos**2 + G * ufl.inner(ufl.dev(eps), ufl.dev(eps))

        psi_c = sigma_c**2 / (2 * E)

        # The damage initiates once the energy density ratio
        # psi_pos / psi_c > 1.0
        # microscopic fracture energy <---> macroscopic stress sigma_c
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
        Vectorized function (faster), to update the field H (crack driving force)
        Handles irreversibility:  H_n+1 = max(H_old, H_new) without loops over cells.
        """
        Q = self.Q
        # Temporary array, storing H of this step
        H_new_array = np.zeros_like(self.H.x.array)

        for label, material in self.materials.items():
            tag = self.label_map[label]
            
            # Finding cells of this material
            entities = self.cell_tags.find(tag)
            if len(entities) == 0:
                continue
            
            H_expr = self.crack_driving_force(u, material)
            expr = dolfinx.fem.Expression(H_expr, Q.element.interpolation_points)

            tmp_func = dolfinx.fem.Function(Q)
            tmp_func.interpolate(expr, entities)
            
            dofs = dolfinx.fem.locate_dofs_topological(Q, self.mesh.topology.dim, entities)
            
            H_new_array[dofs] = tmp_func.x.array[dofs]

        self.H.x.array[:] = np.maximum(self.H.x.array, H_new_array)
        self.H.x.scatter_forward()



    def get_damage_residual(self, D, D_test, material):
        lc = material["lc"]
        H = self.H
        
        res = ( (1 + H) * D * D_test + lc**2 * ufl.dot(ufl.grad(D), ufl.grad(D_test)) - H * D_test ) * ufl.dx
        return res

    def compute_energy_balance(self):
        """
        Compute the energy component to verify the conservation law.
        """
        E_el = 0.0
        E_frac = 0.0
        
        for label, mat in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]
            
            # 1. Elastic energy (already degraded by damage)
            E_el += dolfinx.fem.assemble_scalar(dolfinx.fem.form(self.elastic_energy_density(self.u, mat) * dx))
            
            # 2. Fracture energy
            # Gamma = Gc/2 * (d^2/lc + lc * grad(d)^2)
            if self.on.get("damage"):
                lc = self.dmg_cfg.get("lc") # characteristic length
                Gc = (3.0 * mat["sigma_c"]**2 * self.dmg_cfg["lc"]) / (2.0 * mat["E"]) # Gc, fracture energy
                print(f"Fracture energy: Gc = 3*sigma_c**2/(lc*2*E) = {Gc} J")
                
                gamma = (Gc / 2.0) * ((self.D**2 / lc) + lc * ufl.dot(ufl.grad(self.D), ufl.grad(self.D)))
                E_frac += dolfinx.fem.assemble_scalar(dolfinx.fem.form(gamma * dx))
                
        return E_el, E_frac