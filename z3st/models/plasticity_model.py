# SPDX-License-Identifier: Apache-2.0
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
        
    def _j2_return_map(self, u, material):
        """Trial state and radial return for J2 with linear isotropic hardening.

        Returns ``(sigma_tr, dp, n_flow, mu)``. Both consumers need the same
        predictor and the same plastic multiplier and finish differently:
        sigma_plastic wants the corrected stress, get_plastic_internal_variables
        wants the increments. The two used to compute all of this line for line
        twice.
        """
        eps = self.epsilon(u)
        eps_el_tr = eps - self.ep_n                     # trial elastic strain

        lmbda = material["lmbda"]
        mu = material["G"]
        H = material["hardening_modulus"]

        I = ufl.Identity(eps.ufl_shape[0])
        sigma_tr = lmbda * ufl.tr(eps_el_tr) * I + 2 * mu * eps_el_tr
        s_tr = sigma_tr - (1. / 3.) * ufl.tr(sigma_tr) * I          # deviator
        sigma_eq_tr = ufl.sqrt(1.5 * ufl.inner(s_tr, s_tr))         # von Mises

        # Yield stress grows linearly with the accumulated plastic strain.
        f_val = sigma_eq_tr - (material["yield_strength"] + H * self.p_n)
        dp = ufl.conditional(ufl.gt(f_val, 0), f_val / (3 * mu + H), 0.0)
        n_flow = ufl.conditional(ufl.gt(sigma_eq_tr, 0), s_tr / sigma_eq_tr, s_tr * 0)

        return sigma_tr, dp, n_flow, mu

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
        
        sigma_tr, dp, n_flow, mu = self._j2_return_map(u, material)
        return sigma_tr - 3 * mu * dp * n_flow

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
        _, dp, n_flow, _ = self._j2_return_map(u, material)
        # J2 plastic strain increment = 3/2 * dp * n_flow
        return self.ep_n + 1.5 * dp * n_flow, self.p_n + dp

    def update_plastic_history(self, u):
        """
        Update the history variables (ep_n, p_n) with the converged values.
        """
        print("[PlasticityModel] Updating plastic history...")

        mode = self.plasticity_cfg.get("mode", "j2")
        for name, mat in self.materials.items():
            # Mixed runs: skip materials that carry no plastic law (e.g. an
            # elastic clad next to a plastic fuel) — the unconditional
            # yield_strength read crashed such cases with a KeyError.
            if mode != "custom" and "yield_strength" not in mat:
                continue
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

        # Refresh the history copies ONCE, after every material's cells have
        # been interpolated: the per-material expressions reference ep_n/p_n,
        # so overwriting them inside the loop would hand later materials an
        # already-advanced trial state if cell sets ever overlapped.
        self.ep_n.x.array[:] = self.ep.x.array[:]
        self.p_n.x.array[:] = self.p.x.array[:]

        # Sync ghost
        self.ep_n.x.scatter_forward()
        self.p_n.x.scatter_forward()
        self.ep.x.scatter_forward()
        self.p.x.scatter_forward()