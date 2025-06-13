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

    def set_mechanical_boundary_conditions(self):

        print(f"Setting mechanical boundary conditions...")
        mechanical_bcs_defs = self.boundary_conditions.get("mechanical_bcs", {})

        # self.V_u_sub: independent (collapsed) FunctionSpace
        # self.V_u_map: NumPy array that maps the DoFs from the collapsed space back to the mixed space.
        self.V_u_sub, self.V_u_map = self.W.sub(0).collapse() # need to .collapse() to define BCs / define Functions from it / interpolate
        self.V_u_map = np.array(self.V_u_map, dtype=np.int32)

        dirichlet_bcs_list = []
        
        self.bcs_info.setdefault('neumann_thermal', [])
        self.bcs_info['neumann_thermal'].clear()
    
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

                if bc_type == "Dirichlet":
                    displacement = bc_info.get("displacement")
                    if displacement is None:
                        print(f"  [WARNING] Dirichlet BC on '{mat_type}' for region '{region_name}' has no displacement. Skipping.")
                        continue

                    facets = self.ft.find(region_id)
                    dofs = dolfinx.fem.locate_dofs_topological(self.V_u_sub, self.fdim, facets)
                    
                    bc = self.create_dirichlet_bc(self.mesh, self.V_u_sub, dofs, displacement, self.tdim)

                    dirichlet_bcs_list.append(bc)
                    print(f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {displacement} at region '{region_name}'")

                    print("DirichletBC successfully created")

                    dirichlet_bcs_list.append(bc)
                    print(f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {displacement} at region '{region_name}'")

                elif bc_type == "Neumann":
                    traction_val = bc_info.get("traction")
                    if traction_val is None:
                        print(f"  [WARNING] Neumann BC on '{mat_type}' for region '{region_name}' has no traction. Skipping.")
                        continue
                    
                    traction_constant = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(traction_val))
                    self.bcs_info['traction'].append({'id': region_id, 'value': traction_constant})
                    print(f"  [INFO] Neumann mechanical BC on '{mat_type}' → traction {traction_val} Pa at region '{region_name}'")
                
                else:
                    print(f"  [WARNING] Unknown mechanical BC type '{bc_type}' for '{mat_type}'. Skipping.")
        
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
            print(f"create_dirichlet_bc with Constant")
            return dolfinx.fem.dirichletbc(dolfinx.fem.Constant(mesh, displacement), dofs, V)
        elif displacement.ndim == 2 and displacement.shape == (tdim, 1):
            print(f"  create_dirichlet_bc with Function")
            u_d = dolfinx.fem.Function(V)
            u_d.interpolate(lambda x: np.tile(displacement, (1, x.shape[1])))
            return dolfinx.fem.dirichletbc(u_d, dofs)
        else:
            raise ValueError(f"Formato displacement non valido: shape {displacement.shape}")

    # def set_mechanical_boundary_conditions(self):
    #     print(f"Mechanical boundary conditions")
    #     self.mechanical_bcs = self.boundary_conditions.get("mechanical_bcs", {})

    #     for mat_type, _ in self.materials.items():
    #         bc_list = self.mechanical_bcs.get(mat_type, [])

    #         self.bcs_mechanical[mat_type] = []
    #         self.traction[mat_type] = dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0))
    #         self.traction_label[mat_type] = None

    #         for bc_info in bc_list:
    #             bc_type = bc_info.get("type")
    #             region_name = bc_info.get("region")

    #             if region_name is None or bc_type is None:
    #                 print(f"  [WARNING] Incomplete mechanical BC definition for '{mat_type}'")
    #                 continue

    #             region_id = self.label_map.get(region_name)
    #             if region_id is None:
    #                 print(f"  [WARNING] Region '{region_name}' not found in label_map for mechanical BC.")
    #                 continue

    #             dofs = self.locateMechanicalDofs(region_id)

    #             if bc_type == "Dirichlet":
    #                 displacement = bc_info.get("displacement")
    #                 if displacement is None:
    #                     print(f"  [WARNING] Dirichlet BC on '{mat_type}' but no displacement specified.")
    #                     continue

    #                 if isinstance(displacement, str):
    #                     displacement = eval(displacement)
    #                 elif isinstance(displacement, list):
    #                     displacement = tuple(displacement)

    #                 self.bcs_mechanical[mat_type].append(
    #                     self.fixed_displacement(self.mesh, dofs, self.V_m, displacement)
    #                 )

    #                 print(f"  [INFO] Dirichlet mechanical BC on '{mat_type}' → {displacement} at region '{region_name}'")

    #             elif bc_type == "Neumann":
    #                 traction_scalar = bc_info.get("traction")
    #                 if traction_scalar is None:
    #                     print(f"  [WARNING] Neumann BC on '{mat_type}' but no traction specified.")
    #                     continue

    #                 self.traction[mat_type] = dolfinx.fem.Constant(self.mesh, float(traction_scalar)) * self.normal
    #                 self.traction_label[mat_type] = region_id

    #                 print(f"  [INFO] Neumann mechanical BC on '{mat_type}' → traction {traction_scalar} Pa at region '{region_name}'")

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


    