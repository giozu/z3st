# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import sys

import dolfinx
import dolfinx.fem.petsc
import numpy as np
import ufl
from dolfinx.fem.petsc import NonlinearProblem
from petsc4py import PETSc

from z3st.core.solver import (
    aitken_omega,
    as_bool,
    build_constrained_rigid_nullspace,
    build_rigid_body_nullspace,
)

class MechanicalModel:
    def __init__(self):
        print("[MechanicalModel] initializer")
        self.traction = {}
        self.dirichlet_mechanical = {}

        # --. Mechanical model options --..
        self.mech_cfg = self.input_file.get("mechanical", {})
        # Default linear-solver backend; user can override in input.yaml.
        self.mech_cfg.setdefault("linear_solver", "iterative_hypre")

        print("[MechanicalModel] options loaded from input.yaml:")
        for key, value in self.mech_cfg.items():
            print(f"  {key:<20}: {value}")

    # BC types handled per displacement component (single sub-space).
    _COMPONENT_BC_TYPES = (
        "Dirichlet_x", "Dirichlet_y", "Dirichlet_z",
        "Clamp_x", "Clamp_y", "Clamp_z",
    )
    # BC types that free one axis and block the other two (frictionless slip).
    _SLIP_BC_TYPES = ("Slip_x", "Slip_y", "Slip_z")

    def _regime_normal(self):
        """Facet normal restricted to the regime's displacement components.

        On a 1-D mesh V_m has a single component, so the test function is rank-1
        of size 1 and the traction vector must match: self.normal still carries
        gdim = 3 entries because gmsh stores 3-D node coordinates, and only the
        axial one is physically meaningful on a line. Axisymmetric and 2-D keep
        (r, z) / (x, y); 3-D uses the normal as is.

        The staggered traction update in solver.py needs the same vector.
        """
        if self.mgr.tdim == 1:
            return ufl.as_vector([self.normal[0]])
        if self.regime in ["axisymmetric", "2d"]:
            return ufl.as_vector([self.normal[0], self.normal[1]])
        return self.normal

    def set_mechanical_boundary_conditions(self, V_u):
        """
        Apply mechanical boundary conditions.

        Each definition in ``boundary_conditions["mechanical"]`` is validated
        (region present and known, type present) and then dispatched to the
        handler for its type, which appends to ``self.dirichlet_mechanical`` or
        ``self.traction``. See the ``_add_*`` helpers for the per-type logic.

        Parameters:
            V_u: FunctionSpace (displacement space)
        """
        mechanical_bcs_defs = self.boundary_conditions.get("mechanical", {})

        self._warn_duplicate_mechanical_bcs(mechanical_bcs_defs)

        for name in self.materials:
            self.traction[name] = []
            self.dirichlet_mechanical[name] = []

        for mat_type, bc_list in mechanical_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                if region_name is None or bc_type is None:
                    print(f"  [ERROR] Incomplete mechanical BC definition for '{mat_type}'.")
                    sys.exit(1)

                region_id = self.label_map.get(region_name)
                if region_id is None:
                    print(
                        f"  [ERROR] Region '{region_name}' not found in label_map for mechanical BC."
                    )
                    sys.exit(1)

                facets = self.facet_tags.find(region_id)

                if bc_type == "Dirichlet":
                    self._add_dirichlet_bc(V_u, mat_type, bc_info, region_name, region_id, facets)
                elif bc_type == "Neumann":
                    self._add_neumann_bc(V_u, mat_type, bc_info, region_name, region_id, facets)
                elif bc_type in self._COMPONENT_BC_TYPES:
                    self._add_component_bc(V_u, mat_type, bc_info, region_name, region_id, facets, bc_type)
                elif bc_type in self._SLIP_BC_TYPES:
                    self._add_slip_bc(V_u, mat_type, region_name, facets, bc_type)
                else:
                    print(f"  [ERROR] Unknown mechanical BC type '{bc_type}' for '{mat_type}'.")
                    print(f"  Available are: Dirichlet, Dirichlet_x/y/z, Neumann, Clamp_x/y/z, Slip_x/y/z.")
                    sys.exit(1)

    @staticmethod
    def _warn_duplicate_mechanical_bcs(mechanical_bcs_defs):
        """Warn (but do not abort) on duplicate (region, type) BC definitions."""
        seen_regions = {}
        for mat_type, bc_list in mechanical_bcs_defs.items():
            for bc_info in bc_list:
                key = (bc_info.get("region"), bc_info.get("type"))
                if key in seen_regions:
                    print(
                        f"[WARNING] Duplicate mechanical BC of type '{key[1]}' defined for region '{key[0]}' (previously in '{seen_regions[key]}', now in '{mat_type}')."
                    )
                else:
                    seen_regions[key] = mat_type

    def _add_dirichlet_bc(self, V_u, mat_type, bc_info, region_name, region_id, facets):
        """Full-vector Dirichlet displacement BC (constant or step-dependent)."""
        displacement = bc_info.get("displacement")
        if displacement is None:
            print(
                f"[ERROR] Dirichlet BC on '{mat_type}' for region '{region_name}' has no displacement."
            )
            sys.exit(1)

        # --- interpret displacement input ---
        dim = self.tdim

        # Case 1: single vector [ux, uy, uz] or [u1, u2]
        if isinstance(displacement, (list, tuple)) and all(isinstance(x, (int, float)) for x in displacement):
            if len(displacement) != dim:
                print(f"[ERROR] Displacement vector length {len(displacement)} != mesh dimension {dim}.")
                sys.exit(1)
            raw_value = [displacement]
            print(f"  [INFO] Constant Dirichlet vector ({dim}D) → {displacement}")

        # Case 2: list of vectors over steps
        elif isinstance(displacement, list):
            if all(isinstance(v, (list, tuple)) and len(v) == dim for v in displacement):
                raw_value = displacement
                print(f"  [INFO] Step-dependent Dirichlet list ({dim}D), length {len(displacement)}")

                if len(displacement) != self.n_steps:
                    print(f"[ERROR] BC list length {len(displacement)} != n_steps {self.n_steps}.")
                    sys.exit(1)
            else:
                print(f"[ERROR] All vectors in the list must have length {dim} for a {dim}D mesh.")
                sys.exit(1)

        else:
            print(
                f"[ERROR] Invalid Dirichlet 'displacement' format for region '{region_name}'. "
                f"Must be [ux,uy,uz] or list of such vectors."
            )
            sys.exit(1)

        # --- create constant ---
        initial_disp = raw_value[0]
        # On a 1D mesh V_u is built with a blocked element of shape
        # (1,), which dolfinx treats as a scalar function space.
        # The Constant must then be rank-0, not a length-1 vector,
        # otherwise dirichletbc raises "Rank mismatch between
        # Constant and function space".
        if self.tdim == 1:
            disp_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(initial_disp[0]))
        else:
            disp_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(initial_disp))

        dofs = dolfinx.fem.locate_dofs_topological(V_u, self.fdim, facets)
        bc = dolfinx.fem.dirichletbc(disp_const, dofs, V_u)

        self.dirichlet_mechanical[mat_type].append(
            {
                "id": region_id,
                "value": bc,
                "const": disp_const,
                "raw": raw_value,
            }
        )

        print(
            f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {initial_disp} at region '{region_name}'"
        )

    def _add_neumann_bc(self, V_u, mat_type, bc_info, region_name, region_id, facets):
        """Neumann pressure/traction BC: stored weak-form term T = p * n."""
        traction_value = bc_info.get("traction")

        # Handling list
        if isinstance(traction_value, list):
            raw_value = traction_value
            if len(raw_value) != self.n_steps:
                print(
                    f"[ERROR] Neumann traction list length {len(raw_value)} "
                    f"!= n_steps {self.n_steps} for region '{region_name}'."
                )
                sys.exit(1)
            initial_val = float(traction_value[0])  # starting from step 0
        # Scalar
        else:
            raw_value = [traction_value] * self.n_steps
            initial_val = float(traction_value)

        # scalar constant (Pa)
        traction_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(initial_val))

        n_vec = self._regime_normal()

        # T = p * n
        traction_expr = traction_const * n_vec

        self.traction[mat_type].append({
            "id": region_id,
            "region_name": region_name,
            "value": traction_expr,  # --> in weak form
            "const": traction_const, # --> in the loop
            "raw": raw_value,        # --> list
        })

        print(f"  [INFO] Neumann mechanical BC on '{mat_type}' → {region_name}: {initial_val} Pa (list loaded)")

    def _add_component_bc(self, V_u, mat_type, bc_info, region_name, region_id, facets, bc_type):
        """Component-wise Dirichlet / Clamp on a single axis (x, y, z)."""
        # 2D/3D check for Clamp_z
        if bc_type == "Clamp_z" and self.regime in ["2d", "axisymmetric"]:
            raise ValueError(
                f"\n[ERROR] Boundary condition 'Clamp_z' is not allowed in 2D mode.\n"
                f"        In 2D axisymmetric regime, the axial/vertical component is Y.\n"
                f"        Please use 'Clamp_y' in your boundary_conditions.yaml for region '{region_name}'."
            )

        # Extract value: check 'displacement', then 'value', default to 0.0
        displacement = bc_info.get("displacement")
        if displacement is None:
            displacement = bc_info.get("value", 0.0)

        comp_map = {
            "Dirichlet_x": 0, "Dirichlet_y": 1, "Dirichlet_z": 2,
            "Clamp_x": 0, "Clamp_y": 1, "Clamp_z": 2
        }
        c_idx = comp_map[bc_type]

        if c_idx >= self.tdim:
            print(f"[ERROR] {bc_type} not valid for {self.tdim}D mesh (region '{region_name}').")
            sys.exit(1)

        # Interpret displacement (scalar or list of scalars)
        if isinstance(displacement, (int, float)):
            raw_value = [float(displacement)] * self.n_steps
        elif isinstance(displacement, list):
            if len(displacement) != self.n_steps:
                print(f"[ERROR] {bc_type} list length {len(displacement)} != n_steps {self.n_steps} (region '{region_name}').")
                sys.exit(1)
            raw_value = [float(v) for v in displacement]
        else:
            print(f"[ERROR] Invalid displacement format for {bc_type}. Must be scalar or list of scalars.")
            sys.exit(1)

        disp_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(raw_value[0]))
        # On a 1D mesh V_u has a single component (shape (1,)) and
        # dolfinx refuses V_u.sub(0) -- "Cannot extract subsystem...
        # no subsystems". Since the Clamp_z check above already
        # rejects Clamp_y / Clamp_z for tdim==1, the only valid
        # component is c_idx == 0, which is V_u itself.
        if self.tdim == 1:
            dofs = dolfinx.fem.locate_dofs_topological(V_u, self.fdim, facets)
            bc = dolfinx.fem.dirichletbc(disp_const, dofs, V_u)
        else:
            dofs = dolfinx.fem.locate_dofs_topological(V_u.sub(c_idx), self.fdim, facets)
            bc = dolfinx.fem.dirichletbc(disp_const, dofs, V_u.sub(c_idx))

        self.dirichlet_mechanical[mat_type].append(
            {
                "id": region_id,
                "value": bc,
                "const": disp_const,
                "raw": raw_value,
            }
        )
        print(f"  [INFO] {bc_type} mechanical BC on '{mat_type}' → {raw_value[0]} (first step) at region '{region_name}'")

    def _add_slip_bc(self, V_u, mat_type, region_name, facets, bc_type):
        """
        Frictionless slip: free the named axis and block the other two
        components (e.g. Slip_x leaves u_x free and sets u_y = u_z = 0).
        Each blocked component is appended as a bare DirichletBC, since these
        are not (yet) step-dependent.
        """
        axis_names = ("u_x", "u_y", "u_z")
        free_axis = {"Slip_x": 0, "Slip_y": 1, "Slip_z": 2}[bc_type]
        if free_axis >= self.tdim:
            # Mirrors the Clamp_z guard.
            raise ValueError(
                f"\n[ERROR] Boundary condition '{bc_type}' is not valid on a "
                f"{self.tdim}D mesh: the free axis does not exist, so every "
                f"in-plane component would be blocked (a full clamp). Use "
                f"Slip_x/Slip_y instead."
            )
        constrained = [i for i in (0, 1, 2) if i != free_axis and i < self.tdim]

        for i in constrained:
            V_m_sub = V_u.sub(i)
            boundary_dofs = dolfinx.fem.locate_dofs_topological(
                (V_u, V_m_sub), self.fdim, facets
            )
            boundary_dofs = np.concatenate(boundary_dofs).astype(np.int32)
            bc_i = dolfinx.fem.dirichletbc(
                dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(0.0)),
                boundary_dofs,
                V_m_sub,
            )
            self.dirichlet_mechanical[mat_type].append(bc_i)

        blocked = ", ".join(axis_names[i] for i in constrained)
        print(
            f"  [INFO] {bc_type} mechanical BC on '{mat_type}' → {blocked} = 0.0 at region '{region_name}'"
        )

    def epsilon(self, u):
        """
        Compute the infinitesimal strain tensor epsilon.

        This function supports:
        1. 'axisymmetric': A 2D formulation where the problem is symmetric with respect to the azimutal coordinate.
        The x-coordinate is treated as the radial component (r),
        and the y-coordinate as the axial component (z).
        3. '2D': x-y 2D formulation
        3. '3D' or other: Standard symmetric gradient of the displacement vector.

        Parameters:
            u: Displacement field.

        Returns:
            The 3x3 strain tensor.
        """
        regime = self.regime

        if regime == "axisymmetric":
            # u[0] is radial displacement (u_r), u[1] is axial displacement (u_z)
            r = ufl.SpatialCoordinate(self.mesh)[0]

            # Components of the strain tensor in cylindrical coordinates (r, theta, z)
            eps_rr = u[0].dx(0)  # Normal radial strain
            eps_tt = u[0] / r  # Hoop strain (tangential)
            eps_zz = u[1].dx(1)  # Normal axial strain
            eps_rz = 0.5 * (u[0].dx(1) + u[1].dx(0))  # Shear strain in the r-z plane

            # Return the 3x3 tensor.
            return ufl.as_tensor([[eps_rr, 0.0, eps_rz], [0.0, eps_tt, 0.0], [eps_rz, 0.0, eps_zz]])

        elif regime == "2d":
            # u[0] = x-displacement
            # u[1] = y-displacement
            eps_xx = u[0].dx(0)
            eps_yy = u[1].dx(1)
            eps_xy = 0.5 * (u[0].dx(1) + u[1].dx(0))

            return ufl.as_tensor([[eps_xx, eps_xy, 0.0], [eps_xy, eps_yy, 0.0], [0.0, 0.0, 0.0]])

        elif self.mgr.tdim == 1:
            # 1D mesh: u is a single-component vector (axial displacement only).
            # gmsh always stores 3D node coordinates so the mesh has gdim = 3
            # even when tdim = 1, and ufl.grad(u) returns a (1, 3) tensor that
            # cannot be symmetrized directly. The only meaningful component is
            # the axial derivative, computed explicitly. The full strain tensor
            # is returned padded to 3x3 (axial component in [0,0], zeros
            # elsewhere) so it is compatible with the rest of the pipeline:
            # sigma_mech downstream uses ufl.Identity(eps.ufl_shape[0]) = I_3
            # and the writer expects 3x3 tensors.
            eps_xx = u[0].dx(0)
            return ufl.as_tensor([[eps_xx, 0.0, 0.0],
                                  [0.0,    0.0, 0.0],
                                  [0.0,    0.0, 0.0]])

        else:
            # Default symmetric gradient: 0.5 * (grad(u) + grad(u).T)
            return ufl.sym(ufl.grad(u))

    def deformation_gradient(self, u):
        """
        Deformation gradient F = I + ∇u.

        - axisymmetric: 3x3 in cylindrical (r, θ, z)
        - 2d: 3x3 with F_zz = 1 (plane strain assumption)
        - 3d: F = I + grad(u)

        Parameters:
            u: Displacement field.

        Returns:
            F as a ufl.variable (for use with ufl.diff).
        """
        regime = self.regime

        if regime == "axisymmetric":
            r = ufl.SpatialCoordinate(self.mesh)[0]
            # F in cylindrical coordinates (r, θ, z)
            F_def = ufl.as_tensor([
                [1.0 + u[0].dx(0),  0.0,  u[0].dx(1)],
                [0.0,               1.0 + u[0] / r,  0.0],
                [u[1].dx(0),        0.0,  1.0 + u[1].dx(1)],
            ])

        elif regime == "2d":
            # Plane strain: F_zz = 1
            F_def = ufl.as_tensor([
                [1.0 + u[0].dx(0),  u[0].dx(1),  0.0],
                [u[1].dx(0),        1.0 + u[1].dx(1),  0.0],
                [0.0,               0.0,               1.0],
            ])

        else:
            # 3D: standard
            d = len(u)
            F_def = ufl.Identity(d) + ufl.grad(u)

        return ufl.variable(F_def)

    def sigma_mech(self, u, material):
        """
        Mechanical Cauchy stress σ(u).

        Mode selection (by ``material["constitutive_mode"]``): ``hyperelastic``,
        ``plasticity``, ``custom``, or the default isotropic Lamé
        σ = λ tr(ε) I + 2 G ε (``1d`` reduces to σ = E ε, uniaxial stress).

        ε = sym(∇u) is small strain.
        """

        # The constitutive_mode promotion (lame -> plasticity for materials
        # with yield_strength when plasticity is on) lives in
        # spine.py::load_materials so the material dict is deterministic at
        # load time.
        mode = material.get("constitutive_mode", "lame")
        regime = self.regime

        if mode == "hyperelastic":
            sigma = self.sigma_hyperelastic(u, material)

        elif mode == "plasticity":
            sigma = self.sigma_plastic(u, material)

        elif mode == "custom":
            # Custom material law
            if "stress_function" not in material:
                raise ValueError(f"Material with constitutive_mode='custom' requires 'stress_function' field")

            # Import and use the custom stress function
            import importlib
            stress_func_path = material["stress_function"]

            try:
                module_path, func_name = stress_func_path.rsplit(".", 1) # e.g., z3st.materials.single_crystal_law
                module = importlib.import_module(module_path)
                stress_func = getattr(module, func_name)

                # Call custom function: sigma = f(u, T, material, model)
                T_field = self.T if self.on.get("thermal", False) else None
                sigma = stress_func(u, T_field, material, model=self)

            except Exception as e:
                raise RuntimeError(f"Failed to load/execute custom stress function '{stress_func_path}': {e}")

        else:
            # True 1D structural element (line mesh): uniaxial *stress* state,
            # sigma_yy = sigma_zz = 0. The axial stiffness is therefore the
            # engineering Young's modulus E, giving sigma_11 = E eps_11.
            # epsilon() returns the 3x3 strain padded with only eps_xx != 0, so
            # E * eps has only sigma_xx = E eps_xx and zero transverse stress.
            if regime == "1d":
                eps = self.epsilon(u)
                sigma = material["E"] * eps

            # Default isotropic Lamé
            elif regime == "3d" or regime == "axisymmetric" or regime == "2d":
                eps = self.epsilon(u)
                dim = eps.ufl_shape[0]
                sigma = (
                    material["lmbda"] * ufl.tr(eps) * ufl.Identity(dim) + 2.0 * material["G"] * eps
                )

        if self.on.get("damage", False):
            g_d = self.degradation_function(self.D)
            sigma = g_d * sigma

        return sigma

    def psi_hyperelastic(self, u, material):
        """
        Neo-Hookean strain energy density for hyperelasticity.

        ψ = (μ/2)(I_C - 3) - μ ln(J) + (λ/2)(ln J)²

        Uses deformation_gradient() which handles axisymmetric, 2D, and 3D.

        Parameters:
            u: Displacement field.
            material (dict): Must contain 'lmbda' and 'G'.

        Returns:
            (psi, F_var): strain energy density and the deformation gradient
                          (as ufl.variable, needed for P = dψ/dF).
        """
        F_def = self.deformation_gradient(u)
        C = ufl.variable(F_def.T * F_def)
        Ic = ufl.variable(ufl.tr(C))
        J = ufl.variable(ufl.det(F_def))

        mu = material["G"]
        lmbda = material["lmbda"]

        # Always 3×3 tensor (axisymmetric, 2d, 3d all produce 3×3 F)
        psi = (mu / 2) * (Ic - 3) - mu * ufl.ln(J) + (lmbda / 2) * (ufl.ln(J)) ** 2
        return psi, F_def

    def hyperelastic_residual(self, u, v, material, dx, w):
        """
        Hyperelastic internal force residual via energy derivative.

        Uses S = 2 ∂ψ/∂C (second Piola-Kirchhoff) and E (Green-Lagrange):
            δΠ = ∫ S : δE dx

        Parameters:
            u: Displacement field (Function, not TrialFunction).
            v: TestFunction.
            material (dict): Material properties.
            dx: Integration measure.
            w: Integration weight (e.g. 2πr for axisymmetric).

        Returns:
            UFL form contribution for internal forces.
        """
        psi, F_def = self.psi_hyperelastic(u, material)
        P = ufl.diff(psi, F_def)

        regime = self.regime

        if regime == "axisymmetric":
            r = ufl.SpatialCoordinate(self.mesh)[0]
            # δF consistent with the 3×3 deformation gradient
            grad_v = ufl.as_tensor([
                [v[0].dx(0),  0.0,  v[0].dx(1)],
                [0.0,         v[0] / r,     0.0],
                [v[1].dx(0),  0.0,  v[1].dx(1)],
            ])
        elif regime == "2d":
            grad_v = ufl.as_tensor([
                [v[0].dx(0),  v[0].dx(1),  0.0],
                [v[1].dx(0),  v[1].dx(1),  0.0],
                [0.0,         0.0,         0.0],
            ])
        else:
            grad_v = ufl.grad(v)

        return w * ufl.inner(grad_v, P) * dx

    def sigma_hyperelastic(self, u, material):
        """
        Cauchy stress for the hyperelastic (Neo-Hookean) model.

        σ = (1/J) P Fᵀ

        where P = ∂ψ/∂F is the first Piola-Kirchhoff stress.
        Returns a 3x3 tensor for all regimes (axisymmetric, 2D, 3D).

        Parameters:
            u: Displacement field.
            material (dict): Material properties (needs 'lmbda', 'G').

        Returns:
            Cauchy stress tensor (UFL expression, 3x3).
        """
        psi, F_def = self.psi_hyperelastic(u, material)
        P = ufl.diff(psi, F_def)
        J = ufl.det(F_def)
        return (1.0 / J) * P * F_def.T

    def has_material_eigenstrain(self, material):
        """
        True if the material carries an eigenstrain of its own — swelling, a
        material eigenstrain callable, ... — i.e. one that is *not* the thermal
        eigenstrain.

        (The thermal expansion α(T-T_ref)I is itself an eigenstrain too, but it
        is supplied by the thermal block: it exists only when thermal is active
        and needs the temperature field. This predicate covers the material's
        *own* eigenstrains, which apply regardless of the thermal block.)
        """
        return ("swelling" in material) or ("_eigenstrain_func" in material)

    def applies_eigenstress(self, material):
        """
        Whether the eigenstress -C:ε* must be assembled for this material: when
        the thermal block contributes a thermal eigenstrain, or the material
        carries one of its own. Single predicate so the solver states the intent
        once.
        """
        return self.on.get("thermal", False) or self.has_material_eigenstrain(material)

    def eigenstrain(self, T, material):
        """
        Total inelastic eigenstrain tensor ε* of a material, σ = C : (ε - ε*)

        Every material with a thermal-expansion coefficient contributes the
        thermal eigenstrain α(T - T_ref) I. A material may add further inelastic
        contributions (swelling and densification, creep, ... )
        by exposing an ``eigenstrain`` callable in its card, which
        ``spine.load_materials`` resolves to ``_eigenstrain_func``.
        The result is a UFL tensor, so the Newton tangent stays automatic.
        """
        regime = self.regime
        dim = 1 if regime == "1d" else (3 if regime in ["axisymmetric", "3d", "2d"] else self.tdim)
        I = ufl.Identity(dim)
        eps_star = 0.0 * I  # zero eigenstrain tensor; contributions add below

        # Thermal eigenstrain α(T − T_ref) I
        if T is not None and "alpha" in material and "T_ref" in material:
            eps_star = eps_star + material["alpha"] * (T - material["T_ref"]) * I

        # Constant volumetric swelling ΔV/V from a scalar material field
        # Isotropic, the linear eigenstrain per direction is (ΔV/V)/3 since tr(ε*) = ΔV/V
        s = material.get("swelling")
        if s:
            eps_star = eps_star + (float(s) / 3.0) * I

        # Variable eigenstrains
        fn = material.get("_eigenstrain_func")
        if fn is not None:
            eps_star = eps_star + fn(T, material, model=self, dim=dim)

        return eps_star

    def sigma_th(self, T, material):
        """
        Eigenstress −ℂ : ε* moved to the right-hand side of the momentum
        balance, where ε* is the total inelastic eigenstrain returned by
        :meth:`eigenstrain`. For a purely thermal eigenstrain this reduces to
        the classical σ_th = −(3λ + 2G) α (T − T_ref) I (and −E α (T − T_ref)
        in regime "1d"); fuel swelling/creep enter through the same channel
        once the material supplies them.

        Damage coupling: when damage is active the eigenstress is degraded by
        g(D) = (1−D)² + K, mirroring σ_mech. Without this a fully damaged cell
        (g(D) ~ K ~ 1e-6) loses its stiffness on the LHS but still feels the
        full eigenstress on the RHS, driving the displacement to ~1/K and a
        numerical explosion; degrading both recovers the traction-free
        crack-face limit. (The damage *driver* uses its own "pancake"
        eth_zz = 0 convention in damage_model._thermal_eigenstrain, which
        affects only the ψ⁺ decision, not the displacement field here.)
        """
        eps_star = self.eigenstrain(T, material)
        dim = eps_star.ufl_shape[0]

        # Eigenstress ℂ : ε*. For the engineering 1D bar the stiffness is E
        # (uniaxial); otherwise the isotropic Lamé contraction λ tr(ε*) I + 2G ε*
        # — which for ε* = αΔT I recovers exactly the (3λ + 2G) αΔT thermal form.
        if self.regime == "1d":
            sigma_eig = material["E"] * eps_star
        else:
            sigma_eig = (
                material["lmbda"] * ufl.tr(eps_star) * ufl.Identity(dim)
                + 2.0 * material["G"] * eps_star
            )

        sigma_th = -sigma_eig

        if self.on.get("damage", False):
            sigma_th = self.degradation_function(self.D) * sigma_th

        return sigma_th

    def elastic_energy_density(self, u, material, T=None):
        """
        Elastic strain energy density.

        For linear-elastic (Lame) constitutive::

            psi_el = 0.5 * [lambda * tr(eps_el)^2 + 2G * eps_el : eps_el]

        where eps_el = eps(u) - alpha*(T - T_ref)*I is the elastic strain
        (total strain minus thermal eigenstrain). When T is None or the
        material has no thermal-expansion properties, eps_el = eps(u).

        For other constitutive modes (hyperelastic / plasticity /
        custom), fall back to 0.5 * sigma_mech(u, material) : eps(u) on
        total strain. Those modes don't have thermal coupling wired up in
        z3st, so the simpler formula is used.

        The energy released by cracking is the elastic strain energy: uniform
        thermal expansion in an unconstrained body produces total strain but
        zero elastic strain and zero releasable energy, so the bulk's psi_el
        is ~0.
        """
        mode = material.get("constitutive_mode", "lame")

        if mode == "lame":
            eps = self.epsilon(u)
            if T is not None and "alpha" in material and "T_ref" in material:
                # Thermal eigenstrain. In 2D Cartesian (plane strain or plane
                # stress), the z-component is suppressed because eps_zz is
                # geometrically constrained / sigma_zz is zero -- including
                # eth_zz would produce a spurious bulk psi_el. See the
                # docstring of DamageModel._thermal_eigenstrain, which carries
                # the same logic.
                factor = material["alpha"] * (T - material["T_ref"])
                dim = eps.ufl_shape[0]
                regime = str(getattr(self, "regime", "3d")).lower()
                if regime == "2d" and dim == 3:
                    I_inplane = ufl.as_tensor([
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 0.0],
                    ])
                    eth = factor * I_inplane
                else:
                    eth = factor * ufl.Identity(dim)
                eps_el = eps - eth
            else:
                eps_el = eps
            if str(getattr(self, "regime", "3d")).lower() == "1d":
                # Engineering 1D bar: psi_el = 0.5 E eps_xx^2, consistent with
                # the sigma = E eps law used for regime 1d. eps_el is the 3x3
                # padded strain with only the [0,0] component non-zero, so
                # inner(eps_el, eps_el) = eps_xx^2.
                psi_el = 0.5 * material["E"] * ufl.inner(eps_el, eps_el)
            else:
                psi_el = 0.5 * (
                    material["lmbda"] * ufl.tr(eps_el) ** 2
                    + 2.0 * material["G"] * ufl.inner(eps_el, eps_el)
                )
            # Damage degradation: mirrors sigma_mech / sigma_th. Simple
            # g(D)*psi, not g(D)*psi+ + psi-: a small overcount in compression
            # cells, which carry little D under the hybrid constraint.
            if self.on.get("damage", False):
                psi_el = self.degradation_function(self.D) * psi_el
            return psi_el

        # Non-linear / non-Lame fallback: total-strain energy density.
        # sigma_mech is already g(D)-degraded when damage is active.
        sigma = self.sigma_mech(u, material)
        eps = self.epsilon(u)
        return 0.5 * ufl.inner(sigma, eps)

    # --.. ..- .-.. .-.. --- staggered step --.. ..- .-.. .-.. ---
    # solver.py owns the staggered loop and the services this step calls
    # (get_solver_options, _stagger_residual, _adapt_relax, _bc_objects,
    # _value_at_step, _build_measures, aitken_omega); Spine inherits both.
    def _mechanical_step(self, u_new, u_old, bcs_m, rtol_mech, stag_tol_mech, prev_res_u, T_current):

        u_old.x.array[:] = u_new.x.array

        w = self.weight

        # Creep, or a plasticity / hyperelastic constitutive mode, makes σ(u)
        # nonlinear in u, so the step must go through the SNES path regardless
        # of the configured solver.
        creep_present = any(self.creep_active(m) for m in self.materials.values())
        nonlinear_constitutive = any(
            m.get("constitutive_mode", "lame") in ("plasticity", "hyperelastic")
            for m in self.materials.values()
        )
        linear = (
            self.mech_cfg["solver"] == "linear"
            and not creep_present
            and not nonlinear_constitutive
        )

        # Creep predictor at the current iterate, before assembling. Its
        # change feeds the convergence test.
        creep_pred_change = 0.0
        if creep_present:
            creep_pred_change = self.update_creep_predictor(u_new, T_current)

        # Forms are step-invariant: only Functions (u_new, T_current, creep
        # predictor/state, burnup) and Constants (contact pressure, BC values)
        # change between iterations, used by reference. Build once per step.
        cache = getattr(self, "_mech_cache", None)
        rebuild = (
            cache is None
            or cache["step"] != self.current_step
            or cache["u_new"] is not u_new
            or cache["T"] is not T_current
        )

        bcs_mech = self._bc_objects(self.dirichlet_mechanical)

        if rebuild:
            print("\n[INFO] Assembling mechanical problem...")
            if self.mech_cfg["solver"] == "linear" and creep_present:
                print("  [INFO] creep active → mechanical step promoted to the nonlinear (SNES) path")

            # --- update step-dependent displacement ---
            for _, bc_list in self.dirichlet_mechanical.items():
                for bc in bc_list:
                    # Skip BCs that are Clamp, Slip, etc. (not yet step-dependent)
                    if not isinstance(bc, dict):
                        continue

                    raw = bc.get("raw", None)
                    if isinstance(raw, list):
                        val = self._value_at_step(raw)
                        bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)
                        print(f"  [INFO] Updating Displacement Dirichlet on region {bc['id']} → {val}")

            # --- update step-dependent tractions ---
            for _, bc_list in self.traction.items():
                for bc in bc_list:
                    raw = bc.get("raw", None)

                    if isinstance(raw, list):
                        val = self._value_at_step(raw)
                    elif isinstance(raw, (int, float)):
                        val = raw
                    else:
                        raise RuntimeError(
                            f"Invalid traction 'raw' format (got {type(raw).__name__}: {raw!r}); "
                            f"expected a scalar or a list of length n_steps"
                        )

                    bc["const"].value = np.array(val, dtype=dolfinx.default_scalar_type)
                    print(f"  [INFO] Updating traction on region {bc['id']} → {val} Pa")

                    n_vec = self._regime_normal()

                    bc["value"] = bc["const"] * n_vec

            u_m, v_m = ufl.TrialFunction(self.V_m), ufl.TestFunction(self.V_m)
            a_m, L_m = 0, 0
            F_m = 0

            for label, material in self.materials.items():
                tag = self.label_map[label]
                dx = self.dx_tags[tag]
                print(f"  Building weak form, volume integrals (dx) for {label}, tag = {tag}")

                rho = dolfinx.default_scalar_type(material["rho"])
                g = dolfinx.default_scalar_type(self.g)

                regime = self.regime
                if self.mgr.tdim == 1:
                    body_force = dolfinx.fem.Constant(self.mesh, (-rho * g,))
                elif regime in ["axisymmetric", "2d"]:
                    # 2D: (F_r, F_z) or (F_x, F_y)
                    body_force = dolfinx.fem.Constant(self.mesh, (0.0, -rho * g))
                else:
                    # 3D: (F_x, F_y, F_z)
                    body_force = dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -rho * g))

                if linear:
                    sigma = self.sigma_mech(u_m, material)
                    a_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx
                    L_m += w * ufl.dot(body_force, v_m) * dx
                    # Eigenstress -C:ε* (thermal + material eigenstrains, assembled when the material requires it.
                    if self.applies_eigenstress(material):
                        L_m -= w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx
                else:
                    mode = material.get("constitutive_mode", "lame")
                    if self.creep_active(material):
                        # Condensed implicit creep stress (creep_model.py). The
                        # eigenstrain ε* is inside σ(u) — no separate eigenstress.
                        sigma = self.creep_stress(u_new, material, T_current, self.dt)
                        F_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx
                        F_m -= w * ufl.dot(body_force, v_m) * dx
                    elif mode == "hyperelastic":
                        F_m += self.hyperelastic_residual(u_new, v_m, material, dx, w)
                        F_m -= w * ufl.dot(body_force, v_m) * dx
                        if self.applies_eigenstress(material):
                            F_m += w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx
                    else:
                        sigma = self.sigma_mech(u_new, material)
                        F_m += w * ufl.inner(sigma, self.epsilon(v_m)) * dx - w * ufl.dot(body_force, v_m) * dx
                        # Eigenstress on the residual — mirrors the linear path.
                        if self.applies_eigenstress(material):
                            F_m += w * ufl.inner(self.sigma_th(T_current, material), self.epsilon(v_m)) * dx

            # Traction BCs
            for label in self.materials:
                for bc_info in self.traction[label]:
                    print(f"  Applying mechanical traction on subdomain id = {bc_info['id']}")
                    ds = self.ds_tags[bc_info["id"]]
                    if linear:
                        L_m += w * ufl.dot(bc_info["value"], v_m) * ds
                    else:
                        F_m -= w * ufl.dot(bc_info["value"], v_m) * ds

            # Contact traction (persistent pressure Constant, updated above)
            if self.on.get("contact", False):
                contact_form = self.contact_traction(v_m)
                if linear:
                    L_m += contact_form
                else:
                    F_m -= contact_form

            if linear:
                print("  Linear solver")
                petsc_opts_mech = self.get_solver_options(
                    solver_type=self.mech_cfg["linear_solver"],
                    physics="mechanical",
                    rtol=rtol_mech,
                )
                problem_m = dolfinx.fem.petsc.LinearProblem(
                    a_m,
                    L_m,
                    bcs=bcs_mech,
                    u=u_new,
                    petsc_options=petsc_opts_mech,
                    petsc_options_prefix="mechanical_",
                )
                # Elasticity AMG needs the rigid-body kernel to scale; attach
                # it to the operator (GAMG use it, LU/Hypre ignore it).
                if self.mech_cfg["linear_solver"].startswith("iterative"):
                    problem_m.A.setNearNullSpace(
                        build_rigid_body_nullspace(self.V_m, regime=self.regime)
                    )

                # Opt-in: project the floating rigid-body modes out of the solve
                # for a body the BCs leave rigid-singular. KSP removes the
                # kernel from RHS and solution -> unique minimal-norm
                # displacement. Default off; the standard BC-pinned case has no
                # nullspace and must not get one.
                if as_bool(self.mech_cfg.get("remove_rigid_nullspace", False)):
                    ns = build_constrained_rigid_nullspace(
                        self.V_m, bcs_mech, regime=self.regime
                    )
                    if ns is not None:
                        problem_m.A.setNullSpace(ns)
                        print("  [INFO] rigid-body nullspace removed from mechanical solve")
            else:
                print("  Non-linear solver (SNES Newton)")
                linear_solver = self.mech_cfg.get("linear_solver", "direct_mumps")

                # SNES Newton options + inner linear solver
                if linear_solver == "direct_mumps":
                    petsc_opts_mech = {
                        "snes_type": "newtonls",
                        "snes_linesearch_type": "basic",
                        "snes_atol": rtol_mech,
                        "snes_rtol": rtol_mech,
                        "snes_max_it": int(self.mech_cfg.get("snes_max_it", 50)),
                        "ksp_type": "preonly",
                        "pc_type": "lu",
                        "pc_factor_mat_solver_type": "mumps",
                    }
                else:
                    # Iterative inner solver (AMG / HYPRE)
                    ksp_opts = self.get_solver_options(
                        solver_type=linear_solver,
                        physics="mechanical",
                        rtol=rtol_mech,
                    )
                    petsc_opts_mech = {
                        "snes_type": "newtonls",
                        "snes_linesearch_type": "bt",
                        "snes_atol": rtol_mech,
                        "snes_rtol": rtol_mech,
                        "snes_max_it": 100,
                        "snes_divergence_tolerance": 1e10,
                        **ksp_opts,
                    }

                problem_m = NonlinearProblem(
                    F_m,
                    u_new,
                    bcs=bcs_mech,
                    petsc_options=petsc_opts_mech,
                    petsc_options_prefix="elasticity_",
                )

            self._mech_cache = {
                "step": self.current_step,
                "u_new": u_new,
                "T": T_current,
                "problem": problem_m,
            }

        problem_m = self._mech_cache["problem"]
        dolfinx.fem.set_bc(u_new.x.array, bcs_mech)
        problem_m.solve()

        # Penalty contact: measure the gap from the raw solve and set the
        # pressure used by the next solve.
        if self.on.get("contact", False):
            self.update_contact_pressure(u_new)

        # Relax. With Aitken Δ² enabled the relaxation factor is recomputed
        # each iteration from the last two raw residuals R_k = ũ_k − u_old_k:
        #   ω_{k+1} = −ω_k · (R_{k−1} · ΔR)/|ΔR|²,  ΔR = R_k − R_{k−1},
        # clamped to [relax_min, relax_max]. Dot products are global: restricted
        # to owned dofs (ghosts would be double-counted) and allreduce'd, so omega
        # is rank-independent under MPI. In serial this reduces to the local dot.
        if getattr(self, "relax_aitken", False):
            R = u_new.x.array - u_old.x.array
            R_prev = getattr(self, "_aitken_R_prev", None)
            no = self.V_m.dofmap.index_map.size_local * self.V_m.dofmap.index_map_bs
            omega = aitken_omega(
                R, R_prev, float(getattr(self, "_aitken_omega", self.relax_u)),
                self.mesh.comm, no, self.relax_min, self.relax_max)
            self._aitken_R_prev = R.copy()
            self._aitken_omega = omega
            self.relax_u = omega
            print(f"  [aitken] relax_u={omega:.3f}")

        u_new.x.array[:] = self.relax_u * u_new.x.array + (1 - self.relax_u) * u_old.x.array
        dolfinx.fem.set_bc(u_new.x.array, bcs_mech)

        conv_mech, norm_du, rel_norm_du, res_curr = self._stagger_residual(
            u_new, u_old, self.mech_cfg, stag_tol_mech, "u")

        # The creep predictor must be consistent with u as well — |Δu| alone
        # can pass on the first iteration of a step while Δγ₀ is still moving.
        if creep_present:
            print(f"  [creep] predictor rel change = {creep_pred_change:.3e}")
            pred_tol = max(stag_tol_mech, 1e-8)
            conv_mech = conv_mech and creep_pred_change < pred_tol

        # Heuristic grow/shrink controller — superseded by Aitken when enabled
        if self.relax_adaptive and not getattr(self, "relax_aitken", False):
            prev_res_u = self._adapt_relax("u", res_curr, prev_res_u)

        return conv_mech, norm_du, rel_norm_du, prev_res_u
