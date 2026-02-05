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

        self.dirichlet_damage = {}

    @staticmethod
    def degradation_function(D, K=1e-6):
        """
        Stress degradation function g_d(D) = (1-D)**2 + K

        :param D: scalar field variable. D=0 corresponds to the undamaged state, D=1 to the fully damaged state.
        :param K: numerical regularization parameter
        """
        return (1 - D) ** 2 + K

    def crack_driving_force(self, u, material):
        """
        Compute the dimensionless driving force H for the damage formulations
        """

        # Material properties
        lam = material["lmbda"]
        G = material["G"]
        E = material["E"]

        lc = self.dmg_cfg.get("lc")
        sigma_c = material["sigma_c"]
        Gc = 256/27 * lc * sigma_c**2 / E

        zeta = material["zeta"]

        damage_type = self.dmg_cfg["type"]

        if damage_type == "AT2":

            # Spectral
            psi_pos = self.psi_miehe_spectral(u, material)
            H = (2.0 * lc / Gc) * psi_pos

            # Volumetric
            # lam = material["lmbda"]
            # G = material["G"]
            # eps = self.epsilon(u)
            # tr_eps = ufl.tr(eps)
            # tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
            # psi_pos = 0.5 * lam * tr_eps_pos**2 + G * ufl.inner(ufl.dev(eps), ufl.dev(eps))

            return H

        elif damage_type == "AT1":
            
            # Definitions
            eps = self.epsilon(u)
            
            # Volumetric-Deviatoric Split (Amor et al.)
            tr_eps = ufl.tr(eps)
            tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps)) # Macauley bracket pos
            
            # No cracking under pure hydrostatic compression
            psi_pos = 0.5 * lam * tr_eps_pos**2 + G * ufl.inner(ufl.dev(eps), ufl.dev(eps))

            return psi_pos

    @staticmethod
    def gamma_density(D, grad_D, lc):

        damage_type = self.dmg_cfg["type"]

        if damage_type == "AT2":
            """Fracture energy density gamma(D, ∇D) = 0.5*(D**2/lc + lc*|∇D|²)."""
            return 0.5 * (D**2 / lc + lc * ufl.dot(grad_D, grad_D))

        elif damage_type == "AT1":
            """
            Fracture energy density: gamma(D, ∇D) = (1 / 4*cw) * (w(D)/lc + lc*|∇D|²)
            """
            
            return D / lc + lc * ufl.dot(grad_D, grad_D)

    def psi_miehe_spectral(self, u, material):

        lmbda = material["lmbda"]
        mu = material["G"]
        
        eps = self.epsilon(u)
                
        eps_xx = eps[0, 0]
        eps_yy = eps[1, 1]
        eps_xy = eps[0, 1]
        
        tr_eps = eps_xx + eps_yy
        
        tol = 1.0e-16
        R = ufl.sqrt(((eps_xx - eps_yy) / 2.0)**2 + eps_xy**2 + tol)
        
        eig1 = tr_eps / 2.0 + R
        eig2 = tr_eps / 2.0 - R
        
        # <x>_+ = (x + |x|) / 2
        eig1_pos = 0.5 * (eig1 + abs(eig1))
        eig2_pos = 0.5 * (eig2 + abs(eig2))
        
        tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
                
        psi_pos = (0.5 * lmbda * tr_eps_pos**2) + (mu * (eig1_pos**2 + eig2_pos**2))
        
        return psi_pos

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
            if len(entities) == 0: continue
            
            H_expr = self.crack_driving_force(u, material)
            expr = dolfinx.fem.Expression(H_expr, Q.element.interpolation_points)

            tmp_func = dolfinx.fem.Function(Q)
            tmp_func.interpolate(expr, entities)
            
            dofs = dolfinx.fem.locate_dofs_topological(Q, self.mesh.topology.dim, entities)
            
            H_new_array[dofs] = tmp_func.x.array[dofs]

        damage_type = self.dmg_cfg["type"]
        
        if damage_type == "AT1":
            self.H.x.array[:] = H_new_array
        elif damage_type == "AT2":
            self.H.x.array[:] = np.maximum(self.H.x.array, H_new_array)
        
        self.H.x.scatter_forward()

    def compute_energy_balance(self, u):
        """
        Compute the energy component to verify the conservation law.
        """
        E_el = 0.0
        E_frac = 0.0
        damage_type = self.dmg_cfg["type"]

        for label, mat in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]
            
            # 1. Elastic energy (already degraded by damage)
            E_el += dolfinx.fem.assemble_scalar(dolfinx.fem.form(self.elastic_energy_density(u, mat) * dx))
            
            # 2. Fracture energy
            lc = self.dmg_cfg.get("lc") # characteristic length
            
            if damage_type == "AT1":
                Gc = (8.0 * mat["sigma_c"]**2 * self.dmg_cfg["lc"]) / (3.0 * mat["E"])
            elif damage_type == "AT2":
                Gc = (256.0 * lc * mat["sigma_c"]**2) / (27.0 * mat["E"])
            
            print(f"Fracture energy ({self.dmg_cfg['type']}): Gc = {Gc} J")

            gamma = (Gc / 2.0) * ((self.D**2 / lc) + lc * ufl.dot(ufl.grad(self.D), ufl.grad(self.D)))
            E_frac += dolfinx.fem.assemble_scalar(dolfinx.fem.form(gamma * dx))
                
        return E_el, E_frac

    def set_damage_boundary_conditions(self, V_damage):
            print("\nSetting damage boundary conditions...")
            damage_bcs_defs = self.boundary_conditions.get("damage", {})

            for name in self.materials:
                self.dirichlet_damage[name] = []

            for mat_type, bc_list in damage_bcs_defs.items():
                for bc_info in bc_list:
                    region_name = bc_info.get("region")
                    bc_type = bc_info.get("type")
                    
                    val_d = bc_info.get("value")

                    if region_name is None or bc_type is None or val_d is None:
                        print(f"  [ERROR] Incomplete damage BC definition for '{mat_type}'.")
                        continue

                    region_id = self.label_map.get(region_name)
                    if region_id is None:
                        print(f"  [ERROR] Region '{region_name}' not found in label_map.")
                        continue

                    facets = self.facet_tags.find(region_id)
                    
                    if bc_type == "Dirichlet":
                        d_const = dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(val_d))
                        
                        dofs = dolfinx.fem.locate_dofs_topological(V_damage, self.fdim, facets)
                        
                        bc = dolfinx.fem.dirichletbc(d_const, dofs, V_damage)

                        self.dirichlet_damage[mat_type].append({
                            "id": region_id,
                            "value": bc,
                            "const": d_const
                        })

                        print(f"  [INFO] Dirichlet damage BC on '{mat_type}' → D = {val_d} at region '{region_name}'")