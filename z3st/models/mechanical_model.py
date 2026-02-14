# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import sys

import dolfinx
import numpy as np
import ufl
from petsc4py import PETSc


class MechanicalModel:
    def __init__(self):
        print("[MechanicalModel] initializer")
        self.traction = {}
        self.dirichlet_mechanical = {}

        # --. Mechanical model options --..
        self.mech_cfg = self.input_file.get("mechanical", {})

        print("[MechanicalModel] options loaded from input.yaml:")
        for key, value in self.mech_cfg.items():
            print(f"  {key:<20}: {value}")

    def set_mechanical_boundary_conditions(self, V_u_sub, V_u_map=None):
        """
        Apply mechanical boundary conditions for both staggered and mixed cases.

        If V_u_map is None → staggered (non-mixed spaces).
        If V_u_map is provided → mixed (collapsed subspace).

        To do:
        - ERROR if a region is assigned more than once
        """
        mechanical_bcs_defs = self.boundary_conditions.get("mechanical", {})

        seen_regions = {}

        for mat_type, bc_list in mechanical_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                key = (region_name, bc_type)
                if key in seen_regions:
                    print(
                        f"[WARNING] Duplicate mechanical BC of type '{bc_type}' defined for region '{region_name}' (previously in '{seen_regions[key]}', now in '{mat_type}')."
                    )
                else:
                    seen_regions[key] = mat_type

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
                facet_dim = self.fdim

                # --. Dirichlet --..
                if bc_type == "Dirichlet":
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
                    disp_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(initial_disp))

                    dofs = dolfinx.fem.locate_dofs_topological(V_u_sub, self.fdim, facets)
                    bc = dolfinx.fem.dirichletbc(disp_const, dofs, V_u_sub)

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

                # --. Neumann --..
                elif bc_type == "Neumann":
                    traction_value = bc_info.get("traction")
                    
                    # Handling list
                    if isinstance(traction_value, list):
                        raw_value = traction_value
                        initial_val = float(traction_value[0]) # starting from step 0
                    # Scalar
                    else:
                        raw_value = [traction_value] * self.n_steps
                        initial_val = float(traction_value)

                    # scalar constant (Pa)
                    traction_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(initial_val))

                    # normal vector according to the mechanical regime
                    regime = self.regime
                    if regime in ["axisymmetric", "2d", "plane_stress"]:
                        n_vec = ufl.as_vector([self.normal[0], self.normal[1]])
                    else:
                        n_vec = self.normal

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
                    
                # --. Slip (double component-wise blocking) --..
                elif bc_type == "Slip_x":
                    for i in [1, 2]:  # constrain u_y, u_z
                        V_m_sub = V_u_sub.sub(i)
                        boundary_dofs = dolfinx.fem.locate_dofs_topological(
                            (V_u_sub, V_m_sub), facet_dim, facets
                        )
                        boundary_dofs = np.concatenate(boundary_dofs).astype(np.int32)
                        bc_i = dolfinx.fem.dirichletbc(
                            dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(0.0)),
                            boundary_dofs,
                            V_m_sub,
                        )
                        self.dirichlet_mechanical[mat_type].append(bc_i)

                    print(
                        f"  [INFO] Slip_x mechanical BC on '{mat_type}' → u_y, u_z = 0.0 at region '{region_name}'"
                    )

                elif bc_type == "Slip_y":
                    for i in [0, 2]:  # constrain u_x, u_z
                        V_m_sub = V_u_sub.sub(i)
                        boundary_dofs = dolfinx.fem.locate_dofs_topological(
                            (V_u_sub, V_m_sub), facet_dim, facets
                        )
                        boundary_dofs = np.concatenate(boundary_dofs).astype(np.int32)
                        bc_i = dolfinx.fem.dirichletbc(
                            dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(0.0)),
                            boundary_dofs,
                            V_m_sub,
                        )
                        self.dirichlet_mechanical[mat_type].append(bc_i)

                    print(
                        f"  [INFO] Slip_y mechanical BC on '{mat_type}' → u_x, u_z = 0.0 at region '{region_name}'"
                    )

                elif bc_type == "Slip_z":
                    for i in [0, 1]:  # constrain u_x, u_y
                        V_m_sub = V_u_sub.sub(i)
                        boundary_dofs = dolfinx.fem.locate_dofs_topological(
                            (V_u_sub, V_m_sub), facet_dim, facets
                        )
                        boundary_dofs = np.concatenate(boundary_dofs).astype(np.int32)
                        bc_i = dolfinx.fem.dirichletbc(
                            dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(0.0)),
                            boundary_dofs,
                            V_m_sub,
                        )
                        self.dirichlet_mechanical[mat_type].append(bc_i)

                    print(
                        f"  [INFO] Slip_z mechanical BC on '{mat_type}' →  u_x, u_y = 0.0 at region '{region_name}'"
                    )

                # --. Clamp (single component-wise blocking) --..
                elif bc_type == "Clamp_x":
                    boundary_dofs_x = dolfinx.fem.locate_dofs_topological(
                        V_u_sub.sub(0), self.fdim, self.facet_tags.find(region_id)
                    )
                    bcx = dolfinx.fem.dirichletbc(
                        dolfinx.default_scalar_type(0), boundary_dofs_x, V_u_sub.sub(0)
                    )

                    self.dirichlet_mechanical[mat_type].append(bcx)

                    print(
                        f"  [INFO] Clamp_x mechanical BC on '{mat_type}' → Clamp_x at region '{region_name}'"
                    )

                elif bc_type == "Clamp_y":

                    val = bc_info.get("value", 0.0)

                    boundary_dofs_y = dolfinx.fem.locate_dofs_topological(
                        V_u_sub.sub(1), self.fdim, self.facet_tags.find(region_id)
                    )
                    bcy = dolfinx.fem.dirichletbc(
                        dolfinx.default_scalar_type(val), boundary_dofs_y, V_u_sub.sub(1)
                    )

                    self.dirichlet_mechanical[mat_type].append(bcy)

                    print(
                        f"  [INFO] Clamp_y mechanical BC on '{mat_type}' → Clamp_y at region '{region_name}'"
                    )

                elif bc_type == "Clamp_z":

                    regime = self.regime
                    if regime == "2d" or regime == "axisymmetric":
                        raise ValueError(
                            f"\n[ERROR] Boundary condition 'Clamp_z' is not allowed in 2D mode.\n"
                            f"        In 2D axisymmetric regime, the axial/vertical component is Y.\n"
                            f"        Please use 'Clamp_y' in your boundary_conditions.yaml for region '{region_name}'."
                        )

                    val = bc_info.get("value", 0.0)

                    boundary_dofs_z = dolfinx.fem.locate_dofs_topological(
                        V_u_sub.sub(2), self.fdim, self.facet_tags.find(region_id)
                    )
                    bcz = dolfinx.fem.dirichletbc(
                        dolfinx.default_scalar_type(val), boundary_dofs_z, V_u_sub.sub(2)
                    )
                    self.dirichlet_mechanical[mat_type].append(bcz)

                    print(
                        f"  [INFO] Clamp_z mechanical BC on '{mat_type}' → Clamp_z at region '{region_name}'"
                    )

                else:
                    print(f"  [ERROR] Unknown mechanical BC type '{bc_type}' for '{mat_type}'.")
                    print(f"  Available are: Dirichlet, Neumann, Clamp_x/y/z, Slip_x/y/z.")
                    sys.exit(1)

    def create_dirichlet_bc(self, mesh, V, dofs, displacement, tdim):
        """
        DirichletBC on a vectorial space, using Constant or Functiont.

        Parameters:
            mesh: dolfinx.mesh.Mesh
            V: dolfinx.fem.FunctionSpace (collapsed)
            dofs: np.ndarray
            displacement: list, tuple, np.ndarray
            tdim: int

        Returns:
            dolfinx.fem.DirichletBC
        """
        displacement = np.array(displacement, dtype=PETSc.ScalarType)

        if displacement.ndim == 1 and displacement.size == tdim:
            # print(f"  create_dirichlet_bc with Constant")
            return dolfinx.fem.dirichletbc(dolfinx.fem.Constant(mesh, displacement), dofs, V)
        elif displacement.ndim == 2 and displacement.shape == (tdim, 1):
            # print(f"  create_dirichlet_bc with Function")
            u_d = dolfinx.fem.Function(V)
            u_d.interpolate(lambda x: np.tile(displacement, (1, x.shape[1])))
            return dolfinx.fem.dirichletbc(u_d, dofs)
        else:
            raise ValueError(f"Formato displacement non valido: shape {displacement.shape}")

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

        else:
            # Default symmetric gradient: 0.5 * (grad(u) + grad(u).T)
            return ufl.sym(ufl.grad(u))

    def sigma_mech(self, u, material):
        """
        Mechanical Cauchy stress σ(u) with two selectable constitutive routes.

        Mode selection (by material dict)
        --.--.--.--.-

        - If ``material["constitutive"] == "voigt"``:

        σ = C_voigt · ε_voigt  (6×6 · 6×1) → mapped back to 3×3.

        * If ``material["C_matrix"]`` (6×6) is provided, use it (anisotropy allowed).
        * Else, build isotropic C from (λ, G) with shear blocks = 2G.

        Voigt order used here: [xx, yy, zz, yz, xz, xy] with **no factor 2** on shear strains.


        - Else (default = ``"lame"``):

        σ = λ tr(ε) I + 2 G ε


        Plane-stress handling
        --.--.-----

        If ``self.mech_regime == "plane_stress"``, the returned tensor is modified
        to enforce σ_zz = 0 in an x–y plane-stress sense (same behaviour as before).


        Notes
        -----
        * Assumes 3D (tdim == 3). For 2D, adapt Voigt packing/unpacking.
        * ε = sym(∇u) is small strain.

        """

        # mode = material["constitutive_mode"]
        mode = material.get("constitutive_mode", "lame")
        regime = self.regime

        if mode == "voigt":
            # small strain
            eps = self.epsilon(u)

            # ε in Voigt (6x1), order: [xx, yy, zz, yz, xz, xy]
            eps_voigt = ufl.as_vector(
                [
                    eps[0, 0],
                    eps[1, 1],
                    eps[2, 2],
                    eps[1, 2],  # yz
                    eps[0, 2],  # xz
                    eps[0, 1],  # xy
                ]
            )

            # Elasticity matrix C (6x6)
            if "C_matrix" in material and material["C_matrix"] is not None:
                C_user = material.get("C_matrix", None)
                if C_user is not None:
                    C_np = np.asarray(C_user, dtype=dolfinx.default_scalar_type)
                    if C_np.shape != (6, 6):
                        raise ValueError(f"C_matrix must be 6x6, got {C_np.shape}")
                    C = ufl.as_matrix(C_np.tolist())

            else:
                lmbda = material["lmbda"]
                G = material["G"]

                # isotropic, homogeneous
                C = ufl.as_matrix(
                    [
                        [lmbda + 2 * G, lmbda, lmbda, 0, 0, 0],
                        [lmbda, lmbda + 2 * G, lmbda, 0, 0, 0],
                        [lmbda, lmbda, lmbda + 2 * G, 0, 0, 0],
                        [0, 0, 0, 2 * G, 0, 0],  # yz
                        [0, 0, 0, 0, 2 * G, 0],  # xz
                        [0, 0, 0, 0, 0, 2 * G],  # xy
                    ]
                )

            # σ in Voigt, then map back to 3x3
            sigma_voigt = C * eps_voigt
            sigma = ufl.as_tensor(
                [
                    [sigma_voigt[0], sigma_voigt[5], sigma_voigt[4]],
                    [sigma_voigt[5], sigma_voigt[1], sigma_voigt[3]],
                    [sigma_voigt[4], sigma_voigt[3], sigma_voigt[2]],
                ]
            )

        else:
            # Plane-stress reduction (x–y plane)
            if regime == "plane_stress":

                # Modified Lame parameter for plane stress
                lmbda_ps = (
                    2 * material["G"] * material["lmbda"] / (material["lmbda"] + 2 * material["G"])
                )

                eps = self.epsilon(u)
                sigma = (
                    lmbda_ps * ufl.tr(eps) * ufl.Identity(self.tdim) + 2.0 * material["G"] * eps
                )

                s_xx = sigma[0, 0]
                s_xy = sigma[0, 1]
                s_yy = sigma[1, 1]
                return ufl.as_tensor(
                    [
                        [s_xx, s_xy, 0.0],
                        [s_xy, s_yy, 0.0],
                        [0.0, 0.0, 0.0],
                    ]
                )

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

    def sigma_th(self, T, material):
        """
        Returns the thermal stress tensor σ_th(T) = - (3λ + 2G) α (T - T_ref) I
        for a given temperature field T and material properties.

        Negative sign convention:
            compressive thermal stress for T > T_ref.

        Parameters:
            T (dolfinx.fem.Function): The temperature field.
            material (dict): Material properties.


        Returns:
            dolfinx.fem.Tensor: The thermal stress tensor.

        """
        regime = self.regime
        dim = 3 if regime in ["axisymmetric", "3d", "2d"] else self.tdim

        return (
            -(3 * material["lmbda"] + 2 * material["G"])
            * material["alpha"]
            * (T - material["T_ref"])
            * ufl.Identity(dim)
        )

    def elastic_energy_density(self, u, material):
        """
        Elastic strain energy density = 1/2 * sigma_mech(u) : epsilon(u).
        """

        sigma = self.sigma_mech(u, material)
        eps = self.epsilon(u)

        return 0.5 * ufl.inner(sigma, eps)

    def zero_displacement(self, mesh, dofs, V):
        """
        Returns a DirichletBC with zero displacement vector.

        Parameters:
            mesh (dolfinx.Mesh): The mesh object.
            dofs (array-like): Degrees of freedom where BC is applied.
            V (FunctionSpace): The function space.

        Returns:
            dolfinx.fem.DirichletBC: The zero displacement boundary condition.
        """
        return dolfinx.fem.dirichletbc(dolfinx.fem.Constant(mesh, (0.0, 0.0, 0.0)), dofs, V)

    def fixed_displacement(self, mesh, dofs, V, displacement):
        """
        Returns a DirichletBC with a specified displacement vector.

        Parameters:
            mesh (dolfinx.Mesh): The mesh object.
            dofs (array-like): Degrees of freedom where BC is applied.
            V (FunctionSpace): The function space (vector-valued).
            displacement (tuple or list): The displacement vector (e.g. (0.0, 0.0, 0.0)).

        Returns:
            dolfinx.fem.DirichletBC: The displacement boundary condition.
        """
        return dolfinx.fem.dirichletbc(dolfinx.fem.Constant(mesh, displacement), dofs, V)
