# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import numpy as np
import ufl
from mpi4py import MPI


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

    def _thermal_eigenstrain(self, dim, material, T):
        """Return the thermal eigenstrain tensor used in the damage driver.

        Returns None if T is missing or the material has no thermal-expansion
        properties. The damage driving force must be evaluated on the
        ELASTIC (mechanical) strain eps_el = eps(u) - eth, not on the total
        strain eps(u). Otherwise uniform thermal expansion in the bulk
        produces a spurious psi_pos that drives damage everywhere.

        Regime-dependent z-component handling:

        - In 3D and axisymmetric, every diagonal strain component is
          dynamic (computed from u). The full isotropic eigenstrain
          alpha*(T-T_ref)*I is correct.

        - In 2D Cartesian PLANE STRAIN (regime: 2d in z3st), eps_zz is
          forced to 0 by the geometric constraint. Subtracting a non-zero
          eth_zz = alpha*(T-T_ref) would create eps_el_zz = -alpha*(T-T_ref),
          whose deviatoric component injects a spurious psi_pos
          (= (2/3) * G * alpha^2 * dT^2) throughout the bulk -- about
          5.6 MJ/m^3 for UO2 at dT = 760 K, right at the AT1 threshold for
          sigma_c = 2 GPa. This is a 2D-approximation artifact (the
          uniform sigma_zz from blocked z-thermal-expansion is a real
          stress in an infinite cylinder but cannot drive radial cracks).
          We therefore zero out the z-component of eth in plane strain.

        - PLANE STRESS (regime: plane_stress) has sigma_zz = 0 by
          construction, so the z-thermal-expansion creates no real stress
          either. Same treatment.
        """
        if T is None or "alpha" not in material or "T_ref" not in material:
            return None
        factor = material["alpha"] * (T - material["T_ref"])

        regime = str(getattr(self, "regime", "3d")).lower()
        if regime in ("2d", "plane_stress") and dim == 3:
            # Eigenstrain active in the in-plane components only.
            I_inplane = ufl.as_tensor([
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
            ])
            return factor * I_inplane

        return factor * ufl.Identity(dim)

    def crack_driving_force(self, u, material, T=None):
        """
        Compute the driving force H stored in self.H for the active damage model.

        AT2: H = (2*lc/Gc) * psi_pos     (dimensionless, normalized form)
        AT1: H = psi_pos                 (physical energy density, J/m^3)

        If T is provided, the elastic strain (total strain minus thermal
        eigenstrain) is used in the split, so thermal expansion does not
        appear as a damage driver.
        """

        lc = float(self.dmg_cfg["lc"])
        Gc = material["Gc"]
        damage_type = self.dmg_cfg["type"]

        psi_pos, _ = self.psi_split(u, material, T=T)

        if damage_type == "AT2":
            return (2.0 * lc / Gc) * psi_pos
        elif damage_type == "AT1":
            return psi_pos

    def gamma_density(self, D, grad_D, lc):

        damage_type = self.dmg_cfg["type"]

        if damage_type == "AT2":
            # Fracture energy density gamma(D, ∇D) = 0.5*(D**2/lc + lc*|∇D|²)
            return 0.5 * (D**2 / lc + lc * ufl.dot(grad_D, grad_D))

        elif damage_type == "AT1":
            # Fracture energy density: gamma(D, ∇D) = (1 / 4*cw) * (w(D)/lc + lc*|∇D|²)
            return D / lc + lc * ufl.dot(grad_D, grad_D)

    def psi_miehe_spectral(self, u, material, T=None):
        """
        Miehe spectral split of the elastic strain energy density.

        Returns the tuple (psi_pos, psi_neg) with
            psi_pos = (lambda/2) <tr eps_el>_+^2 + mu * sum_i <eps_el_i>_+^2
            psi_neg = (lambda/2) <tr eps_el>_-^2 + mu * sum_i <eps_el_i>_-^2
        where <.>_+ and <.>_- are the positive and negative parts and
        eps_el = eps(u) - alpha*(T - T_ref)*I is the elastic (mechanical)
        strain (the total strain when T is None or thermal props are absent).
        """

        lmbda = material["lmbda"]
        mu = material["G"]

        eps = self.epsilon(u)
        dim = eps.ufl_shape[0]
        tol = 1.0e-16

        # Subtract thermal eigenstrain so the split sees the elastic strain only.
        eth = self._thermal_eigenstrain(dim, material, T)
        if eth is not None:
            eps = eps - eth

        if dim == 2:
            eps_xx = eps[0, 0]
            eps_yy = eps[1, 1]
            eps_xy = eps[0, 1]

            tr_eps = eps_xx + eps_yy

            R = ufl.sqrt(((eps_xx - eps_yy) / 2.0)**2 + eps_xy**2 + tol)

            eig1 = tr_eps / 2.0 + R
            eig2 = tr_eps / 2.0 - R

            eig1_pos = 0.5 * (eig1 + abs(eig1))
            eig2_pos = 0.5 * (eig2 + abs(eig2))
            eig1_neg = 0.5 * (eig1 - abs(eig1))
            eig2_neg = 0.5 * (eig2 - abs(eig2))

            tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
            tr_eps_neg = 0.5 * (tr_eps - abs(tr_eps))

            psi_pos = (0.5 * lmbda * tr_eps_pos**2) + (mu * (eig1_pos**2 + eig2_pos**2))
            psi_neg = (0.5 * lmbda * tr_eps_neg**2) + (mu * (eig1_neg**2 + eig2_neg**2))

        else:
            # 3D: Cardano's formula for eigenvalues of symmetric 3x3 tensor
            e00 = eps[0, 0]; e11 = eps[1, 1]; e22 = eps[2, 2]
            e01 = eps[0, 1]; e02 = eps[0, 2]; e12 = eps[1, 2]

            tr_eps = e00 + e11 + e22

            # Invariants of eps
            I1 = tr_eps
            I2 = (e00*e11 + e11*e22 + e00*e22
                   - e01**2 - e02**2 - e12**2)
            I3 = (e00*(e11*e22 - e12**2)
                   - e01*(e01*e22 - e12*e02)
                   + e02*(e01*e12 - e11*e02))

            # Deviatoric invariants for Cardano
            p = (I1**2 - 3.0*I2) / 9.0
            q = (2.0*I1**3 - 9.0*I1*I2 + 27.0*I3) / 54.0

            s = ufl.sqrt(p + tol)

            # cos(theta) = q / s^3, clamped
            arg = q / (s**3 + tol)
            # Smooth clamp via conditional
            arg_c = ufl.conditional(ufl.gt(arg, 1.0 - tol), 1.0 - tol,
                     ufl.conditional(ufl.lt(arg, -1.0 + tol), -1.0 + tol, arg))
            theta = ufl.acos(arg_c) / 3.0

            eig1 = I1 / 3.0 + 2.0 * s * ufl.cos(theta)
            eig2 = I1 / 3.0 + 2.0 * s * ufl.cos(theta - 2.0 * ufl.pi / 3.0)
            eig3 = I1 / 3.0 + 2.0 * s * ufl.cos(theta + 2.0 * ufl.pi / 3.0)

            eig1_pos = 0.5 * (eig1 + abs(eig1))
            eig2_pos = 0.5 * (eig2 + abs(eig2))
            eig3_pos = 0.5 * (eig3 + abs(eig3))
            eig1_neg = 0.5 * (eig1 - abs(eig1))
            eig2_neg = 0.5 * (eig2 - abs(eig2))
            eig3_neg = 0.5 * (eig3 - abs(eig3))

            tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
            tr_eps_neg = 0.5 * (tr_eps - abs(tr_eps))

            psi_pos = (0.5 * lmbda * tr_eps_pos**2) + (mu * (eig1_pos**2 + eig2_pos**2 + eig3_pos**2))
            psi_neg = (0.5 * lmbda * tr_eps_neg**2) + (mu * (eig1_neg**2 + eig2_neg**2 + eig3_neg**2))

        return psi_pos, psi_neg

    def psi_amor_split(self, u, material, T=None):
        """
        Amor (volumetric/deviatoric) split of the elastic strain energy density.

        Returns (psi_pos, psi_neg) with
            psi_pos = (lambda/2) <tr eps_el>_+^2 + mu * dev(eps_el):dev(eps_el)
            psi_neg = (lambda/2) <tr eps_el>_-^2
        where eps_el = eps(u) - alpha*(T - T_ref)*I is the elastic strain
        (total strain when T is None). Note: this uses lambda rather than
        the n-dimensional bulk modulus K_n. Same convention as the existing
        AT1 driving force in z3st.
        """
        lam = material["lmbda"]
        G = material["G"]
        eps = self.epsilon(u)
        dim = eps.ufl_shape[0]

        # Subtract thermal eigenstrain so the split sees the elastic strain only.
        eth = self._thermal_eigenstrain(dim, material, T)
        if eth is not None:
            eps = eps - eth

        tr_eps = ufl.tr(eps)
        tr_eps_pos = 0.5 * (tr_eps + abs(tr_eps))
        tr_eps_neg = 0.5 * (tr_eps - abs(tr_eps))

        psi_pos = 0.5 * lam * tr_eps_pos**2 + G * ufl.inner(ufl.dev(eps), ufl.dev(eps))
        psi_neg = 0.5 * lam * tr_eps_neg**2

        return psi_pos, psi_neg

    def psi_split(self, u, material, T=None):
        """
        Dispatch to the elastic energy split matching the active damage model:
        Miehe spectral for AT2, Amor volumetric/deviatoric for AT1.
        Returns (psi_pos, psi_neg) in physical units (J/m^3).

        T (optional): the current temperature field. When provided, the split
        is evaluated on the elastic (mechanical) strain so thermal expansion
        does not contribute spuriously.
        """
        damage_type = self.dmg_cfg["type"]
        if damage_type == "AT2":
            return self.psi_miehe_spectral(u, material, T=T)
        elif damage_type == "AT1":
            return self.psi_amor_split(u, material, T=T)
        else:
            raise ValueError(f"Unknown damage type '{damage_type}'")

    def update_history(self, u, T=None):
        """
        Vectorized update of the crack driving force history field self.H.

        AT2 enforces irreversibility:  H_{n+1} = max(H_n, H_new).
        AT1 stores the current value (irreversibility is enforced on D itself).

        Hybrid (Ambati et al. 2015) constraint:
            if dmg_cfg["hybrid_constraint"] is True (default), then in cells where
            psi_neg > psi_pos the new contribution to H is set to zero, so the
            damage field cannot grow in compression-dominated regions. For
            monotonic loads (e.g. thermal shock) this matches Eq. (27c) of the
            paper. It is conservative under unloading: it suppresses crack
            re-opening rather than enforcing pointwise crack closure.
        """

        Q = self.Q
        H_new_array = np.zeros_like(self.H.x.array)

        use_hybrid_constraint = self.dmg_cfg.get("hybrid_constraint", True)

        psi_pos_fn = dolfinx.fem.Function(Q) if use_hybrid_constraint else None
        psi_neg_fn = dolfinx.fem.Function(Q) if use_hybrid_constraint else None
        tmp_func = dolfinx.fem.Function(Q)

        for label, material in self.materials.items():
            tag = self.label_map[label]

            entities = self.cell_tags.find(tag)
            if len(entities) == 0:
                continue

            H_expr = self.crack_driving_force(u, material, T=T)
            expr = dolfinx.fem.Expression(H_expr, Q.element.interpolation_points)
            tmp_func.interpolate(expr, entities)

            dofs = dolfinx.fem.locate_dofs_topological(Q, self.mesh.topology.dim, entities)
            H_new_array[dofs] = tmp_func.x.array[dofs]

            if use_hybrid_constraint:
                psi_pos_expr, psi_neg_expr = self.psi_split(u, material, T=T)

                expr_pos = dolfinx.fem.Expression(psi_pos_expr, Q.element.interpolation_points)
                expr_neg = dolfinx.fem.Expression(psi_neg_expr, Q.element.interpolation_points)

                psi_pos_fn.interpolate(expr_pos, entities)
                psi_neg_fn.interpolate(expr_neg, entities)

                # Cells where the negative-energy part dominates: suppress damage growth.
                local_dofs = dofs
                neg_mask = psi_neg_fn.x.array[local_dofs] > psi_pos_fn.x.array[local_dofs]
                H_new_array[local_dofs[neg_mask]] = 0.0

        damage_type = self.dmg_cfg["type"]

        if damage_type == "AT1":
            self.H.x.array[:] = H_new_array
        elif damage_type == "AT2":
            self.H.x.array[:] = np.maximum(self.H.x.array, H_new_array)

        self.H.x.scatter_forward()

    def compute_energy_balance(self, u):
        """
        Compute total elastic and fracture energies for the conservation diagnostic.

        Surface (fracture) energy density:
            AT2:  gamma = (Gc/2) * (D^2/lc + lc * |grad D|^2)
            AT1:  gamma = (3*Gc/8) * (D/lc + lc * |grad D|^2)

        Both integrals use the regime-dependent integration weight
        (`self.weight = 2*pi*r` for axisymmetric, 1 otherwise) so that the
        reported energies are true volume integrals and not per-(r,z)-area
        proxies. Without this factor, an axisymmetric run reports energies
        that are off by ~pi*R compared to the real 3D-volume value.

        The elastic-energy density is evaluated on the elastic strain
        (eps(u) - alpha*(T - T_ref)*I), so uniform thermal expansion in the
        bulk does not produce a spurious E_el offset.
        """
        E_el = 0.0
        E_frac = 0.0
        damage_type = self.dmg_cfg["type"]
        lc = self.dmg_cfg.get("lc")

        # Integration weight: 2*pi*r for axisymmetric, 1 elsewhere.
        # _build_measures sets self.weight when solve_staggered runs; fall
        # back to 1.0 if compute_energy_balance is called before any solve.
        w = getattr(self, "weight", 1.0)

        grad_D2 = ufl.dot(ufl.grad(self.D), ufl.grad(self.D))

        T_field = getattr(self, "T", None)

        for label, mat in self.materials.items():
            tag = self.label_map[label]
            dx = self.dx_tags[tag]

            # 1. Elastic energy (degraded by damage), evaluated on eps_el.
            E_el += dolfinx.fem.assemble_scalar(
                dolfinx.fem.form(self.elastic_energy_density(u, mat, T=T_field) * w * dx)
            )

            # 2. Fracture energy.
            Gc = mat["Gc"]
            if damage_type == "AT2":
                gamma = (Gc / 2.0) * (self.D**2 / lc + lc * grad_D2)
            elif damage_type == "AT1":
                gamma = (3.0 * Gc / 8.0) * (self.D / lc + lc * grad_D2)
            else:
                raise ValueError(f"Unknown damage type '{damage_type}'")

            E_frac += dolfinx.fem.assemble_scalar(dolfinx.fem.form(gamma * w * dx))

        # assemble_scalar returns per-rank partial sums; reduce so every
        # rank returns the same global value and downstream callers see
        # the physical energy (not the rank-0 share).
        comm = self.mesh.comm
        E_el = comm.allreduce(E_el, op=MPI.SUM)
        E_frac = comm.allreduce(E_frac, op=MPI.SUM)

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