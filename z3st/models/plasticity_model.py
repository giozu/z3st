# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import ufl
import numpy as np
from petsc4py import PETSc

class PlasticityModel:
    def __init__(self):        
        """
        Initialize internal variables for plasticity model.
        """

        print("[PlasticityModel] initializer")
        self.plasticity_cfg = self.input_file.get("plasticity", {})
        
        # Scalar function for cumulative plastic strain
        self.p = dolfinx.fem.Function(self.Q_pl, name="CumulativePlasticStrain")
        self.p.x.array[:] = 0.0

        # Tensor function for plastic strain tensor
        self.ep = dolfinx.fem.Function(self.V_pl, name="PlasticStrainTensor")
        self.ep.x.array[:] = 0.0
        
        # Scalar function for cumulative plastic strain (old)
        self.p_n = dolfinx.fem.Function(self.Q_pl, name="CumulativePlasticStrain_old")
        self.p_n.x.array[:] = 0.0
        
        # Tensor function for plastic strain tensor (old)
        self.ep_n = dolfinx.fem.Function(self.V_pl, name="PlasticStrainTensor_old")
        self.ep_n.x.array[:] = 0.0
        
    def sigma_plastic(self, u, material):
        """
        Calculate stress and update internal variables using return mapping (J2 plasticity).
        
        Returns:
            sigma (ufl.Tensor): The stress tensor expression.
            eff_stress (ufl.Scalar): Von Mises stress.
            
        Note:
            This function defines the constitutive relation in UFL for the Newton solver.
            Ideally, we should return the stress `sigma` that depends on `u` and history `ep_n`, `p_n`.
            
            J2 Plasticity with isotropic hardening.
        """
        
        # Elastic predictor
        eps = self.epsilon(u) # Total strain
        eps_el_tr = eps - self.ep_n # Trial elastic strain
        
        lmbda = material["lmbda"]
        mu = material["G"]
        
        # Sigma trial
        dim = eps.ufl_shape[0]
        I = ufl.Identity(dim)
        
        sigma_tr = lmbda * ufl.tr(eps_el_tr) * I + 2 * mu * eps_el_tr
        
        # Yield condition
        # Deviatoric stress (trial)
        s_dev_tr = sigma_tr - (1./3.) * ufl.tr(sigma_tr) * I
        
        # Von Mises equivalent stress (trial)
        sigma_eq_tr = ufl.sqrt(1.5 * ufl.inner(s_dev_tr, s_dev_tr))
        
        # Yield stress
        # sigma_y = sigma_0 + H * p # Linearly plastic J2, yield strength increases with plastic strain
        sigma_y = material["yield_strength"]
        H = material["hardening_modulus"]
        sigma_y_n = sigma_y + H * self.p_n

        # Check yield condition: f = sigma_eq_tr - sigma_y(p_n) <= 0
        f_val = sigma_eq_tr - sigma_y_n
        
        # Return mapping
        dp = ufl.conditional(ufl.gt(f_val, 0), f_val / (3*mu + H), 0.0)
        
        # Update
        n_flow = ufl.conditional(ufl.gt(sigma_eq_tr, 0), s_dev_tr / sigma_eq_tr, s_dev_tr*0)
        
        sigma = sigma_tr - 3*mu*dp*n_flow
        
        return sigma

    def get_plastic_internal_variables(self, u, material):
        """
        Returns UFL expressions for the updated internal variables (ep_new, p_new) for the current displacement u.
        Used for updating history variables after convergence.

        Supports both standard J2 plasticity and custom crystal plasticity models.
        """
        # Check if custom plasticity mode
        mode = self.plasticity_cfg.get("mode", "j2")

        if mode == "custom":
            # Import custom function for crystal plasticity
            import importlib
            stress_func_path = material.get("stress_function", "")
            module_path = ".".join(stress_func_path.split(".")[:-1])
            module = importlib.import_module(module_path)

            if hasattr(module, 'get_cp_internal_variables'):
                get_vars_func = getattr(module, 'get_cp_internal_variables')
                T_field = self.T if self.on.get("thermal", False) else None
                ep_new = get_vars_func(u, T_field, material, model=self)

                # Calculate p_new from ep_new
                # NOTE: For J2 plasticity, p = sqrt(3/2 * ep:ep) is exact.
                #       For crystal plasticity, this gives p ≈ 0.866*γ instead of p = γ.
                #       This is acceptable if p is not used in the constitutive law.
                #       For slip-system hardening, p should be computed as cumulative slip.
                p_new = ufl.sqrt(1.5 * ufl.inner(ep_new, ep_new))
                return ep_new, p_new
            else:
                raise AttributeError(f"Module {module_path} does not have 'get_cp_internal_variables' function")

        # Standard J2 plasticity
        eps = self.epsilon(u)
        eps_el_tr = eps - self.ep_n

        lmbda = material["lmbda"]
        mu = material["G"]

        dim = eps.ufl_shape[0]
        I = ufl.Identity(dim)
        sigma_tr = lmbda * ufl.tr(eps_el_tr) * I + 2 * mu * eps_el_tr

        s_tr = sigma_tr - (1./3.) * ufl.tr(sigma_tr) * I
        sigma_eq_tr = ufl.sqrt(1.5 * ufl.inner(s_tr, s_tr))

        sigma_y = material["yield_strength"]
        H = material["hardening_modulus"]

        sigma_y_n = sigma_y + H * self.p_n
        f_val = sigma_eq_tr - sigma_y_n

        dp = ufl.conditional(ufl.gt(f_val, 0), f_val / (3*mu + H), 0.0)
        n_flow = ufl.conditional(ufl.gt(sigma_eq_tr, 0), s_tr / sigma_eq_tr, s_tr*0)

        # Updated internal variables
        p_new = self.p_n + dp
        ep_new = self.ep_n + 1.5 * dp * n_flow # plastic strain increment (J2) = 3/2 * dp * n_flow

        return ep_new, p_new

    def update_plastic_history(self, u):
        """
        Update the history variables (ep_n, p_n) with the converged values.
        """
        print("[PlasticityModel] Updating plastic history...")
        
        for name, mat in self.materials.items():
            ep_expr, p_expr = self.get_plastic_internal_variables(u, mat)
            
            tag = self.label_map[name]
            cells = self.cell_tags.find(tag)
            
            V_ep = self.ep.function_space
            V_p = self.p.function_space
            
            expr_ep = dolfinx.fem.Expression(ep_expr, V_ep.element.interpolation_points)
            expr_p = dolfinx.fem.Expression(p_expr, V_p.element.interpolation_points)
            
            # Update current state variables (using old _n values)
            self.ep.interpolate(expr_ep, cells)
            self.p.interpolate(expr_p, cells)

            # Then update history variables from current
            self.ep_n.x.array[:] = self.ep.x.array[:]
            self.p_n.x.array[:] = self.p.x.array[:]
                                
        # Sync ghost
        self.ep_n.x.scatter_forward()
        self.p_n.x.scatter_forward()
        self.ep.x.scatter_forward()
        self.p.x.scatter_forward()