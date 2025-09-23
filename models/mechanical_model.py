import sys, numpy as np
import dolfinx, ufl
from petsc4py import PETSc

class MechanicalModel:
    def __init__(self):
        print("[MechanicalModel] initializer")
        self.traction = {}
        self.dirichlet_mechanical = {}

        data = self.input_file

        self.mech_options = data.get("mechanical", {}) # parsing
        if not self.mech_options:
            raise ValueError("[MechanicalModel] 'mechanical' missing in input.yaml.")

        print("[MechanicalModel] options loaded from input.yaml:")
        for key, value in self.mech_options.items():
            print(f"  {key:<20}: {value}")
            
    def set_mechanical_boundary_conditions(self, V_u_sub, V_u_map=None):
        """
        Apply mechanical boundary conditions for both staggered and mixed cases.

        If V_u_map is None → staggered (non-mixed spaces).
        If V_u_map is provided → mixed (collapsed subspace).

        To do:
        - ERROR if a region is assigned more than once
        """
        print("\nSetting mechanical boundary conditions...")
        mechanical_bcs_defs = self.boundary_conditions.get("mechanical_bcs", {})
        
        seen_regions = {}

        for mat_type, bc_list in mechanical_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                key = (region_name, bc_type)
                if key in seen_regions:
                    print(f"[WARNING] Duplicate mechanical BC of type '{bc_type}' defined for region '{region_name}' (previously in '{seen_regions[key]}', now in '{mat_type}').")
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
                    print(f"  [ERROR] Region '{region_name}' not found in label_map for mechanical BC.")
                    sys.exit(1)

                facets = self.ft.find(region_id)
                facet_dim = self.fdim

                # -------- Dirichlet (vector full BC) --------
                if bc_type == "Dirichlet":
                    displacement = bc_info.get("displacement")
                    if displacement is None:
                        print(f"[ERROR] Dirichlet BC on '{mat_type}' for region '{region_name}' has no displacement.")
                        sys.exit(1)

                    # locate dofs
                    dofs = dolfinx.fem.locate_dofs_topological(V_u_sub, self.fdim, facets)
                    bc = self.create_dirichlet_bc(self.mesh, V_u_sub, dofs, displacement, self.tdim)
                    self.dirichlet_mechanical[mat_type].append(bc)
                    print(f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {displacement} at region '{region_name}'")

                # -------- Neumann --------
                elif bc_type == "Neumann":
                    traction_value = bc_info.get("traction")
                    if traction_value is None:
                        print(f"  [ERROR] Neumann BC on '{mat_type}' for region '{region_name}' has no traction.")
                        sys.exit(1)

                    # Create a mutable Constant that stores the scalar traction magnitude
                    traction_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(traction_value))

                    # Build the traction expression as (Constant * normal)
                    traction_expr = traction_const * self.normal

                    # Store both: the expression (used in forms) and the Constant (for updates)
                    self.traction[mat_type].append({
                        'id': region_id,
                        'value': traction_expr,   # UFL expression used in weak form
                        'const': traction_const   # PETSc Constant (mutable in incremental loop)
                    })

                    print(f"  [INFO] Neumann mechanical BC on '{mat_type}' → traction {traction_value} Pa at region '{region_name}'")

                # -------- Generic Slip (component-wise blocking) --------
                elif bc_type == "Slip_x":
                    for i in [1, 2]:  # constrain u_y, u_z
                        V_m_sub = V_u_sub.sub(i)
                        boundary_dofs = dolfinx.fem.locate_dofs_topological(
                            (V_u_sub, V_m_sub),
                            facet_dim,
                            facets
                        )
                        boundary_dofs = np.concatenate(boundary_dofs).astype(np.int32)
                        bc_i = dolfinx.fem.dirichletbc(
                            dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(0.0)),
                            boundary_dofs,
                            V_m_sub
                        )
                        self.dirichlet_mechanical[mat_type].append(bc_i)

                    print(f"  [INFO] Slip_x mechanical BC on '{mat_type}' → u_y, u_z = 0.0 at region '{region_name}'")

                elif bc_type == "Slip_z":
                    for i in [0, 1]:  # constrain u_x, u_y
                        V_m_sub = V_u_sub.sub(i)
                        boundary_dofs = dolfinx.fem.locate_dofs_topological(
                            (V_u_sub, V_m_sub),
                            facet_dim,
                            facets
                        )
                        boundary_dofs = np.concatenate(boundary_dofs).astype(np.int32)
                        bc_i = dolfinx.fem.dirichletbc(
                            dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(0.0)),
                            boundary_dofs,
                            V_m_sub
                        )
                        self.dirichlet_mechanical[mat_type].append(bc_i)

                    print(f"  [INFO] Slip_z mechanical BC on '{mat_type}' →  u_x, u_y = 0.0 at region '{region_name}'")

                elif bc_type == "Clamp_x":
                    boundary_dofs_x = dolfinx.fem.locate_dofs_topological(V_u_sub.sub(0), self.fdim, self.ft.find(region_id))
                    bcx = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), boundary_dofs_x, V_u_sub.sub(0))

                    self.dirichlet_mechanical[mat_type].append(bcx)
                    
                    print(f"  [INFO] Clamp_x mechanical BC on '{mat_type}' → Clamp_x at region '{region_name}'")

                elif bc_type == "Clamp_y":
                    boundary_dofs_y = dolfinx.fem.locate_dofs_topological(V_u_sub.sub(1), self.fdim, self.ft.find(region_id))
                    bcy = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), boundary_dofs_y, V_u_sub.sub(1))

                    self.dirichlet_mechanical[mat_type].append(bcy)

                    print(f"  [INFO] Clamp_y mechanical BC on '{mat_type}' → Clamp_y at region '{region_name}'")

                elif bc_type == "Clamp_z":
                    val = bc_info.get("value", 0.0)

                    boundary_dofs_z = dolfinx.fem.locate_dofs_topological(V_u_sub.sub(2), self.fdim, self.ft.find(region_id))
                    bcz = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(val), boundary_dofs_z, V_u_sub.sub(2))
                    self.dirichlet_mechanical[mat_type].append(bcz)

                    print(f"  [INFO] Clamp_z mechanical BC on '{mat_type}' → u_z = {val} at region '{region_name}'")

                else:
                    print(f"  [ERROR] Unknown mechanical BC type '{bc_type}' for '{mat_type}'.")
                    print(f"  Available are: Dirichlet, Neumann, Clamp_x/y/z, Slip_z.")
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

    def constrain_point(self, mesh: dolfinx.mesh.Mesh, 
                        V: dolfinx.fem.FunctionSpace, 
                        point: np.ndarray) -> dolfinx.fem.DirichletBC:
        """
        Finds the mesh node closest to a given point and applies a zero
        Dirichlet boundary condition to it.

        Args:
            mesh: The dolfinx mesh object.
            V: The FunctionSpace on which the boundary condition is defined.
            point: A numpy array of shape (dim,) specifying the target coordinates
                (e.g., np.array([0.1, 0.2, 0.0])).

        Returns:
            A dolfinx.fem.DirichletBC object that constrains the node.
        """
        # 1. Get all node coordinates from the mesh geometry.
        # The shape is (num_nodes, 3) for a 3D mesh.
        all_coords = mesh.geometry.x

        # 2. Ensure the input point has the correct dimension.
        point = np.asarray(point, dtype=mesh.geometry.x.dtype)
        if point.shape[0] != mesh.geometry.dim:
            raise ValueError(f"Point dimension {point.shape[0]} does not match mesh dimension {mesh.geometry.dim}")

        # 3. Calculate the squared Euclidean distance from every node to the target point.
        # Using squared distance avoids a costly square root and gives the same minimum.
        distances_sq = np.sum((all_coords - point)**2, axis=1)

        # 4. Find the index (and the coordinates) of the node with the minimum distance.
        closest_node_idx = np.argmin(distances_sq)
        closest_coord = all_coords[closest_node_idx]

        # Optional: Print a message to inform the user which node was selected.
        # This is very useful for debugging.
        print(f"\nRequested point: {point}. Constraining closest node at: {closest_coord}\n")

        # 5. Find the degrees of freedom (dofs) geometrically located at that specific node.
        # The lambda function identifies the point(s) that are extremely close to our
        # `closest_coord`. np.isclose is used for safe floating-point comparison.
        # The `.all(axis=1)` ensures that all components (x, y, z) of a coordinate match.
        dofs = dolfinx.fem.locate_dofs_geometrical(
            V,
            lambda x: np.isclose(x.T, closest_coord).all(axis=1)
        )

        # 6. Create the Dirichlet Boundary Condition.
        # We constrain the located dofs to a zero vector of the appropriate dimension.
        # The value is wrapped in a dolfinx.fem.Constant.
        u_bc = dolfinx.fem.Constant(mesh, np.zeros(mesh.geometry.dim, dtype=dolfinx.default_scalar_type))
        bc = dolfinx.fem.dirichletbc(u_bc, dofs, V)

        return bc

    def epsilon(self, u):
        # return 0.5 * (ufl.grad(u) + ufl.transpose(ufl.grad(u)))
        return ufl.sym(ufl.grad(u))

    def sigma_mech(self, u, material):
        """
        Mechanical Cauchy stress σ(u) with two selectable constitutive routes:

        Mode selection (by material dict):
        - material["constitutive"] == "voigt":
            σ = C_voigt · ε_voigt   (6x6 · 6x1)  → mapped back to 3x3.
            * If material["C_matrix"] (6x6) is provided, use it (anisotropy allowed).
            * Else build isotropic C from (λ, G) with shear blocks = 2G.
            Voigt order used here: [xx, yy, zz, yz, xz, xy] with NO factor-2 on shear strains.
        - else (default "lame"):
            σ = λ tr(ε) I + 2 G ε

        Plane-stress handling:
        If self.regime == "plane_stress", the returned tensor is modified to enforce σ_zz = 0
        in an x–y plane-stress sense (same behaviour you had).

        Notes
        -----
        * Assumes 3D (tdim == 3). For 2D, adapt Voigt packing/unpacking.
        * ε = sym(∇u) is small strain.
        """
        
        # mode = material["constitutive_mode"]
        mode = material.get("constitutive_mode","lame")

        if mode == "voigt":
            # small strain
            eps = self.epsilon(u)

            # ε in Voigt (6x1), order: [xx, yy, zz, yz, xz, xy]
            eps_voigt = ufl.as_vector([
                eps[0, 0],
                eps[1, 1],
                eps[2, 2],
                eps[1, 2],  # yz
                eps[0, 2],  # xz
                eps[0, 1],  # xy
            ])

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
                G     = material["G"]

                # isotropic, homogeneous
                C = ufl.as_matrix([
                    [lmbda+2*G, lmbda,     lmbda,     0,   0,   0],
                    [lmbda,     lmbda+2*G, lmbda,     0,   0,   0],
                    [lmbda,     lmbda,     lmbda+2*G, 0,   0,   0],
                    [0,         0,         0,       2*G,  0,   0],  # yz
                    [0,         0,         0,         0, 2*G,   0], # xz
                    [0,         0,         0,         0,   0, 2*G]  # xy
                ])

            # σ in Voigt, then map back to 3x3
            sigma_voigt = C * eps_voigt
            sigma = ufl.as_tensor([
                [sigma_voigt[0], sigma_voigt[5], sigma_voigt[4]],
                [sigma_voigt[5], sigma_voigt[1], sigma_voigt[3]],
                [sigma_voigt[4], sigma_voigt[3], sigma_voigt[2]]
            ])

        else:
            # Default isotropic Lamé
            eps = self.epsilon(u)
            sigma = material["lmbda"] * ufl.tr(eps) * ufl.Identity(self.tdim) + 2.0 * material["G"] * eps

        # Optional plane-stress reduction (x–y plane)
        if getattr(self, "regime", "3D") == "plane_stress":
            s_xx = sigma[0, 0]
            s_xy = sigma[0, 1]
            s_yy = sigma[1, 1]
            return ufl.as_tensor([
                [s_xx, s_xy, 0.0],
                [s_xy, s_yy, 0.0],
                [0.0,  0.0,  0.0],
            ])

        return sigma

    def sigma_th(self, T, material):
        return (3 * material["lmbda"] + 2 * material["G"]) * material["alpha"] * (T - material["T_ref"]) * ufl.Identity(self.tdim)

    def elastic_energy_density(self, u, material):
        """
        Elastic strain energy density (small strain linear elasticity).
        Uses w = 1/2 * sigma_mech(u) : epsilon(u).
        Works also with self.regime == 'plane_stress' because sigma_mech() already handles it.
        """

        sigma = self.sigma_mech(u, material)
        eps   = self.epsilon(u)
        
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