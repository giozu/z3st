import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np

class MechanicalModel:
    def __init__(self):
        print("__MechanicalModel initializer__")
        self.a_mechanical = {}
        self.L_mechanical = {}
        self.bcs_mechanical = {}
        self.traction = {}
        self.traction_label = {}
        self.mechanical_residuals = {}
        self.i_mech = {}
        self.u = {}
        self.u_old = {}

    def set_mechanical_boundary_conditions(self, V_u_sub, V_u_map=None):
        """
        Apply mechanical boundary conditions for both staggered and mixed cases.

        If V_u_map is None → staggered (non-mixed spaces).
        If V_u_map is provided → mixed (collapsed subspace).
        """
        print("Setting mechanical boundary conditions...")
        mechanical_bcs_defs = self.boundary_conditions.get("mechanical_bcs", {})

        dirichlet_bcs_list = []
        self.bcs_info.setdefault('traction', [])
        self.bcs_info['traction'].clear()

        for mat_type, bc_list in mechanical_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                if region_name is None or bc_type is None:
                    print(f"  [WARNING] Incomplete mechanical BC definition for '{mat_type}'. Skipping.")
                    continue

                region_id = self.label_map.get(region_name)
                if region_id is None:
                    print(f"  [WARNING] Region '{region_name}' not found in label_map for mechanical BC. Skipping.")
                    continue

                facets = self.ft.find(region_id)
                facet_dim = self.fdim
                # -------- Dirichlet (vector full BC) --------
                if bc_type == "Dirichlet":
                    displacement = bc_info.get("displacement")
                    if displacement is None:
                        print(f"[WARNING] Dirichlet BC on '{mat_type}' for region '{region_name}' has no displacement. Skipping.")
                        continue

                    if V_u_map is None:  # staggered mode
                        dofs = dolfinx.fem.locate_dofs_topological(V_u_sub, self.fdim, facets)
                    else:  # mixed mode
                        dofs_mixed = dolfinx.fem.locate_dofs_topological((self.W.sub(0), V_u_sub), self.fdim, facets)
                        dofs_mixed = np.array(dofs_mixed, dtype=np.int32)
                        dofs = V_u_map[dofs_mixed[:, 1]]

                    bc = self.create_dirichlet_bc(self.mesh, V_u_sub, dofs, displacement, self.tdim)
                    dirichlet_bcs_list.append(bc)
                    print(f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {displacement} at region '{region_name}'")

                # -------- Neumann --------
                elif bc_type == "Neumann":
                    traction_val = bc_info.get("traction")
                    if traction_val is None:
                        print(f"  [WARNING] Neumann BC on '{mat_type}' for region '{region_name}' has no traction. Skipping.")
                        continue

                    traction_constant = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(traction_val)) * self.normal
                    self.bcs_info['traction'].append({'id': region_id, 'value': traction_constant})
                    print(f"  [INFO] Neumann mechanical BC on '{mat_type}' → traction {traction_val} Pa at region '{region_name}'")

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
                        dirichlet_bcs_list.append(bc_i)

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
                        dirichlet_bcs_list.append(bc_i)

                    print(f"  [INFO] Slip_z mechanical BC on '{mat_type}' →  u_x, u_y fixed at region '{region_name}'")

                elif bc_type == "Slip_xy":
                    for i in [2]:  # constrain only u_z
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
                        dirichlet_bcs_list.append(bc_i)

                    print(f"  [INFO] Slip_xy mechanical BC on '{mat_type}' → u_z fixed at region '{region_name}'")

                elif bc_type == "Clamp_x":
                    boundary_dofs_x = dolfinx.fem.locate_dofs_topological(V_u_sub.sub(0), self.fdim, self.ft.find(region_id))
                    bcx = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), boundary_dofs_x, V_u_sub.sub(0))

                    dirichlet_bcs_list.append(bcx)
                    print(f"  [INFO] Clamp_x mechanical BC on '{mat_type}' → Clamp_x at region '{region_name}'")

                elif bc_type == "Clamp_y":
                    boundary_dofs_y = dolfinx.fem.locate_dofs_topological(V_u_sub.sub(1), self.fdim, self.ft.find(region_id))
                    bcx = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), boundary_dofs_y, V_u_sub.sub(1))

                    dirichlet_bcs_list.append(bcx)
                    print(f"  [INFO] Clamp_y mechanical BC on '{mat_type}' → Clamp_y at region '{region_name}'")

                elif bc_type == "Clamp_z":
                    boundary_dofs_z = dolfinx.fem.locate_dofs_topological(V_u_sub.sub(2), self.fdim, self.ft.find(region_id))
                    bcz = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), boundary_dofs_z, V_u_sub.sub(2))

                    dirichlet_bcs_list.append(bcz)
                    print(f"  [INFO] Clamp_z mechanical BC on '{mat_type}' → Clamp_z at region '{region_name}'")
                    
                else:
                    print(f"  [WARNING] Unknown mechanical BC type '{bc_type}' for '{mat_type}'.")
                    print(f"  Available are: Dirichlet, Neumann, Clamp_x/y/z, Slip_xy/z. Skipping.")

        return dirichlet_bcs_list

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
        return material["lmbda"] * ufl.tr(self.epsilon(u)) * ufl.Identity(self.tdim) + 2 * material["mu"] * self.epsilon(u)

    def sigma_th(self, T, material):
        return (3 * material["lmbda"] + 2 * material["mu"]) * material["alpha"] * (T - material["T_ref"]) * ufl.Identity(self.tdim)

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

    def assemble_mechanical_solver(self, material_label):

        material = self.materials[material_label]

        T = self.T[material_label]

        self.dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=self.label_map[material_label])

        self.strain = self.epsilon(self.v_m)
        self.stress = self.sigma_mech(self.u_m, material)

        body_force = dolfinx.fem.Constant(self.mesh, (0, 0, -material["rho"] * self.g))

        self.a_mechanical[material_label] = ufl.inner(self.stress, self.strain) * self.dx_local

        self.L_mechanical[material_label] = (
            ufl.inner(body_force, self.v_m) * self.dx_local
            + ufl.inner(self.sigma_th(T, material), ufl.sym(ufl.grad(self.v_m))) * self.dx_local
            + ufl.inner(self.traction[material_label], self.v_m) * self.ds(self.traction_label[material_label])
        )
 
    def solve_mechanical_problem(self, material_label, max_iter=1000, tol=1e-8):

        a_form = dolfinx.fem.form(self.a_mechanical[material_label])
        L_form = dolfinx.fem.form(self.L_mechanical[material_label])

        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)

        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=self.bcs_mechanical[material_label])
        A.assemble()

        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [self.bcs_mechanical[material_label]])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, self.bcs_mechanical[material_label])

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        # solver.setType(PETSc.KSP.Type.GMRES)
        # solver.getPC().setType("jacobi")
        
        # solver.setType("lu")
        
        solver.setType("gmres")
        solver.getPC().setType("gamg")
        # solver.getPC().setType("hypre")

        solver.setTolerances(rtol=tol, atol=tol, max_it=max_iter)

        residuals = []

        def monitor(ksp, its, rnorm):
            residuals.append(rnorm)

        solver.setMonitor(monitor)

        if self.debug_mode is True:
            print(f"[DEBUG] Solving mechanical problem for {material_label}")

        solver.solve(b, self.u[material_label].x.petsc_vec)
        self.u[material_label].x.scatter_forward()
        self.mechanical_residuals[material_label] = residuals

        return solver.getIterationNumber()


    