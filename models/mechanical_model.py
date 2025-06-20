import sys, numpy as np
import dolfinx, ufl
from petsc4py import PETSc

class MechanicalModel:
    def __init__(self):
        print("__MechanicalModel initializer__")
        self.a_mechanical = {}
        self.L_mechanical = {}
        self.traction = {}
        self.dirichlet_mechanical = {}

    def set_mechanical_boundary_conditions(self, V_u_sub, V_u_map=None):
        """
        Apply mechanical boundary conditions for both staggered and mixed cases.

        If V_u_map is None → staggered (non-mixed spaces).
        If V_u_map is provided → mixed (collapsed subspace).

        To do:
        - ERROR if a region is assigned more than once
        """
        print("Setting mechanical boundary conditions...")
        mechanical_bcs_defs = self.boundary_conditions.get("mechanical_bcs", {})
        
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

                    if V_u_map is None:  # staggered mode
                        dofs = dolfinx.fem.locate_dofs_topological(V_u_sub, self.fdim, facets)
                    else:  # mixed mode
                        dofs_mixed = dolfinx.fem.locate_dofs_topological((self.W.sub(0), V_u_sub), self.fdim, facets)
                        dofs_mixed = np.array(dofs_mixed, dtype=np.int32)
                        dofs = V_u_map[dofs_mixed[:, 1]]

                    bc = self.create_dirichlet_bc(self.mesh, V_u_sub, dofs, displacement, self.tdim)
                    self.dirichlet_mechanical[mat_type].append(bc)
                    print(f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {displacement} at region '{region_name}'")

                # -------- Neumann --------
                elif bc_type == "Neumann":
                    traction_value = bc_info.get("traction")
                    if traction_value is None:
                        print(f"  [ERROR] Neumann BC on '{mat_type}' for region '{region_name}' has no traction.")
                        sys.exit(1)

                    traction_constant = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(traction_value)) * self.normal
                    self.traction[mat_type].append({'id': region_id, 'value': traction_constant})
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

                    print(f"  [INFO] Slip_x mechanical BC on '{mat_type}' → u_y, u_z fixed at region '{region_name}'")

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

                    print(f"  [INFO] Slip_z mechanical BC on '{mat_type}' →  u_x, u_y fixed at region '{region_name}'")

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
        return ufl.sym(ufl.grad(u))

    def sigma_mech(self, u, material):

        stress_tensor = material["lmbda"] * ufl.tr(self.epsilon(u)) * ufl.Identity(self.tdim) + 2 * material["G"] * self.epsilon(u)

        if self.regime == "plane_stress":
            """
            Modified constitutive model for 3D mesh under plane stress assumption in x-y plane.
            The component σ_zz is set to zero by modifying the strain tensor appropriately.
            """

            s_xx = stress_tensor[0,0]
            s_xy = stress_tensor[0,1]
            s_yy = stress_tensor[1,1]

            stress_tensor_plane_stress = ufl.as_tensor([
                [s_xx, s_xy, 0.0],
                [s_xy, s_yy, 0.0],
                [0.0,  0.0,  0.0]
            ])

            return stress_tensor_plane_stress
        
        else:
            return stress_tensor

    def sigma_th(self, T, material):
        return (3 * material["lmbda"] + 2 * material["G"]) * material["alpha"] * (T - material["T_ref"]) * ufl.Identity(self.tdim)

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