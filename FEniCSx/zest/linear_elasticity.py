import numpy as np
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc

from FEniCSx.zest.finiteElementSettings import Mechanics

class LinearElasticity(Mechanics):
    def volume(self):
        """
        Compute the volume of the fuel and cladding.
        """
        
        # Integration measures for each subdomain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
        dx_fuel = dx_subdomain(1)      # Marker ID for fuel region
        dx_cladding = dx_subdomain(2)  # Marker ID for cladding region

        # Compute volumes
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))
        volume_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_fuel))
        volume_cladding = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_cladding))

        # Print results
        print("------------")
        print("Mesh volumes")
        print("------------")
        print(f"Total volume, from mesh = {volume} m3")
        print(f"Fuel volume, from mesh = {volume_fuel} m3")
        print(f"Cladding volume, from mesh = {volume_cladding} m3")

        return volume, volume_fuel, volume_cladding

    def parameters(self, rho, g, lmbda, mu):
        """
        Set material and physical parameters.

        Parameters:
        rho (float): Density.
        g (float): Gravitational acceleration.
        lmbda (float): First Lamé parameter.
        mu (float): Second Lamé parameter.
        """
        self.rho = rho
        self.g = g
        self.lmbda = lmbda
        self.mu = mu

    def compute_material_properties(material):
        G = material["E"] / (2 * (1 + material["nu"]))
        lame1 = material["nu"] * material["E"] / ((1 + material["nu"]) * (1 - 2 * material["nu"]))
        return G, lame1

    def validate_material(material):
        required_keys = ["E", "nu", "alpha", "rho", "geometry"]
        missing = [key for key in required_keys if key not in material]
        if missing:
            raise KeyError(f"Missing keys in material: {missing}")

    def boundary_conditions(self, radius_fi, radius_fo, height_fuel, radius_ci, radius_co, height_cladding, gap_pressure, coolant_pressure):
        """
        Define and set the boundary conditions.
        """

        # ---------
        # Fuel dofs
        # ---------

        self.dofs_D_bottom_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(1))

        self.dofs_D_outer_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(3))

        self.dofs_D_top_fuel = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(4))

        # -------------
        # Cladding dofs
        # -------------

        self.dofs_D_bottom_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(2))

        self.dofs_D_inner_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(5))

        self.dofs_D_outer_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(6))

        self.dofs_D_top_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(7))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Clamped boundary conditions
        u_clamped = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, 0.0)
        )

        self.bc_top_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_top_cladding,
            self.V
        )

        self.bc_top_fuel  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_top_fuel,
            self.V
        )
        
        self.bc_bottom_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_bottom_cladding,
            self.V
        )

        self.bc_bottom_fuel  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_bottom_fuel,
            self.V
        )
        
        self.bc = [self.bc_bottom_fuel, self.bc_bottom_cladding, self.bc_top_fuel, self.bc_top_cladding]
        # self.bc = [self.bc_bottom_fuel, self.bc_bottom_cladding]
        # self.bc = [self.bc_top_fuel, self.bc_top_cladding]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.gap_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(gap_pressure)
        )

        self.coolant_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(coolant_pressure)
        )
        
        self.traction_gap = self.gap_pressure * self.normal
        self.traction_coolant = self.coolant_pressure * self.normal

    def boundary_conditions_thick_cylinder(self, radius_ci, radius_co, height_cladding, gap_pressure, coolant_pressure):
        """
        Define and set the boundary conditions.
        """

        self.bottom_cladding_label = 1
        self.inner_cladding_label = 2
        self.outer_cladding_label = 3
        self.top_cladding_label = 4

        # -------------
        # Cladding dofs
        # -------------

        self.dofs_D_bottom_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.bottom_cladding_label))

        self.dofs_D_inner_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.inner_cladding_label))

        self.dofs_D_outer_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.outer_cladding_label))

        self.dofs_D_top_cladding = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.top_cladding_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Clamped boundary conditions
        u_clamped = dolfinx.fem.Constant(
            self.mesh, (0.0, 0.0, 0.0)
        )

        self.bc_top_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_top_cladding,
            self.V
        )

        self.bc_bottom_cladding  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_D_bottom_cladding,
            self.V
        )

        self.bc = [self.bc_bottom_cladding, self.bc_top_cladding]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.gap_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(gap_pressure)
        )

        self.coolant_pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(coolant_pressure)
        )
        
        self.traction_gap = self.gap_pressure * self.normal
        self.traction_coolant = self.coolant_pressure * self.normal

    def boundary_conditions_full_cylinder(self, pressure):
        """
        Define and set the boundary conditions.
        """

        self.bottom_label = 1
        self.lateral_label = 2
        self.top_label= 3

        # ----
        # Dofs
        # ----

        self.dofs_bottom  = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.bottom_label))

        self.dofs_lateral = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.lateral_label))

        self.dofs_top     = dolfinx.fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(self.top_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Clamped boundary conditions
        u_clamped = dolfinx.fem.Constant(
            self.mesh, (0.0, 0.0, 0.0)
        )

        self.bc_top  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_top,
            self.V
        )

        self.bc_bottom  = dolfinx.fem.dirichletbc(
            u_clamped,
            self.dofs_bottom,
            self.V
        )

        self.bc = [self.bc_bottom, self.bc_top]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.pressure = dolfinx.fem.Constant(
            self.mesh,
            dolfinx.default_scalar_type(pressure)
        )
        
        self.traction = self.pressure * self.normal

    def assemble_system(self):
        """
        Assemble the system for the linear elasticity problem.
        """

        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        def epsilon(u):
            # return ufl.sym(ufl.grad(u))
            return 0.5 * (ufl.grad(u) + ufl.grad(u).T)

        # -----------
        # Body forces
        # -----------
        self.f = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, - self.rho * self.g) # or ScalarType((0.0, 0.0, - self.rho * self.g))
        )

        # Mesure "ds"
        self.ds = ufl.Measure(
            "ds",
            domain=self.mesh,
            subdomain_data=self.ft
        )

        # ds_fo = self.ds(self.outer_fuel)        # or ds(10) with subdomain_data=self.ft
        # ds_ci = self.ds(self.inner_cladding)    # or ds(9) with subdomain_data=self.ft
        # ds_co = self.ds(self.outer_cladding)    # or ds(9) with subdomain_data=self.ft

        ds_fo = self.ds(3) # fuel outer
        ds_ci = self.ds(5) # cladding inner
        ds_co = self.ds(6) # cladding outer

        ds_fo = self.ds(3) # fuel outer
        ds_ci = self.ds(self.inner_cladding_label) # cladding inner
        ds_co = self.ds(self.outer_cladding_label) # cladding outer

        # -------------------
        # Material properties
        # -------------------
        
        # lame1 = dolfinx.fem.Function(self.V_scalar)
        # cells_0 = dolfinx.mesh.locate_entities(self.mesh, self.tdim, 1)
        # cells_1 = dolfinx.mesh.locate_entities(self.mesh, self.tdim, 2)

        self.a = ufl.inner(sigma(self.u, self.lmbda, self.mu), epsilon(self.v)) * ufl.dx
        self.L_body_forces = ufl.inner(self.f, self.v) * ufl.dx
        self.L_tractions = ufl.inner(self.traction_coolant, self.v) * ds_co + ufl.inner(self.traction_gap, self.v) * ds_ci # + ufl.inner(self.traction_gap, self.v) * ds_fo

        self.L = self.L_body_forces + self.L_tractions

        # Extraction of the linear and bi-linear form
        F = self.a - self.L # full variational form
        a, L = ufl.system(F)

        # We generate and compile the C code for these expressions using dolfinx.fem.form
        a_compiled = dolfinx.fem.form(a)
        L_compiled = dolfinx.fem.form(L)

        # Assemble the stiffness matrix
        A = dolfinx.fem.petsc.assemble_matrix(a_compiled)
        A.assemble()

        # Save the matrix to a text file
        viewer = PETSc.Viewer().createASCII("stiffness_matrix.txt")
        A.view(viewer)

        print(f"Matrix A is symmetric after assembly: {A.isSymmetric(1e-5)}")
        # If the matrix A isn't symmetric we need to exclude a whole class of iterative solvers (Conjugate Gradient).

    def assemble_system_full_cylinder(self):
        """
        Assemble the system for the linear elasticity problem.
        """

        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        def epsilon(u):
            return 0.5 * (ufl.grad(u) + ufl.grad(u).T)

        # -----------
        # Body forces
        # -----------
        self.f = dolfinx.fem.Constant(
            self.mesh,
            (0.0, 0.0, - self.rho * self.g)
        )

        # Mesure "ds"
        self.ds = ufl.Measure(
            "ds",
            domain=self.mesh,
            subdomain_data=self.ft
        )

        self.a = ufl.inner(sigma(self.u, self.lmbda, self.mu), epsilon(self.v)) * ufl.dx

        self.L_body_forces = ufl.inner(self.f, self.v) * ufl.dx
        self.L_tractions = ufl.inner(self.traction, self.v) * self.ds(self.lateral_label)
        
        self.L = self.L_body_forces + self.L_tractions

        # Extraction of the linear and bi-linear form
        F = self.a - self.L # full variational form
        a, L = ufl.system(F)

        # We generate and compile the C code for these expressions using dolfinx.fem.form
        a_compiled = dolfinx.fem.form(a)
        L_compiled = dolfinx.fem.form(L)

        # Assemble the stiffness matrix
        A = dolfinx.fem.petsc.assemble_matrix(a_compiled)
        A.assemble()

        # Save the matrix to a text file
        viewer = PETSc.Viewer().createASCII("stiffness_matrix.txt")
        A.view(viewer)

        print(f"Matrix A is symmetric after assembly: {A.isSymmetric(1e-5)}")
        # If the matrix A isn't symmetric we need to exclude a whole class of iterative solvers (Conjugate Gradient).

    def solve(self):
        """
        Solve the linear elasticity problem.
        """

        # PETSc options
        # LU factorization, a direct solver method (similar to Gaussian elimination).
        # For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES is better.

        # petsc_options = {
        #     "ksp_type": "preonly",
        #     "pc_type": "lu",
        #     "pc_factor_mat_solver_type": "mumps",
        # }

        petsc_options = {
            "ksp_type": "gmres",
            "ksp_rtol": 1e-6,
            "ksp_atol": 1e-10,
            "ksp_max_it": 1000,
            "pc_type": "none",
        }

        jit_options = {
            "timeout": 20,
        }

        # Linear problem definition
        self.problem = dolfinx.fem.petsc.LinearProblem(
            self.a,
            self.L,
            bcs=self.bc,
            petsc_options=petsc_options,
            jit_options=jit_options,
        )

        # Solve the problem
        self.displacement = self.problem.solve()
        # self.displacement.x.scatter_forward() # scatter_forward is necessary for dolfinx ?

        # Compute stress and strain fields
        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

        def epsilon(u):
            # return ufl.sym(ufl.grad(u))
            return 0.5*(ufl.grad(u) + ufl.grad(u).T)

        self.stress_field = sigma(self.displacement, self.lmbda, self.mu)
        self.strain_field = epsilon(self.displacement)

        # Solver options
        solver = self.problem.solver
        viewer = PETSc.Viewer().createASCII("log.solver")
        solver.view(viewer)
        # solver_output = open("lu_output.txt", "r")
        # for line in solver_output.readlines():
        #     print(line)

    def volumetric_stress(self):
        """
        Compute the volumetric component of the stress field.
        The volumetric stress is the trace of the stress tensor divided by the geometric dimension.
        It is related to the change in volume of the material, and the hydrostatic pressure.
        """

        I = ufl.Identity(self.gdim)
        volumetric_stress = ufl.tr(self.stress_field) * I / self.gdim

        return volumetric_stress

    def deviatoric_stress(self):
        """
        Compute the deviatoric component of the stress field.
        The deviatoric stress is the stress tensor minus the volumetric stress.
        It is related to the shape change of the material, and the shear stress.
        """
        
        I = ufl.Identity(self.gdim)
        deviatoric_stress = self.stress_field - ufl.tr(self.stress_field) * I / self.gdim

        return deviatoric_stress

    def get_results(self):
        """
        Get the computed displacement, stress, and strain fields.

        Returns:
        displacement (Function): The displacement field.
        stress_field (ufl.Expr): The stress field.
        strain_field (ufl.Expr): The strain field.
        """
        return self.displacement, self.stress_field, self.strain_field
    
    def test_problem_3d_rect(self):

        Lx = 0.1
        Ly = 1.0
        Lz = 2.0

        domain = dolfinx.mesh.create_box(MPI.COMM_WORLD, [np.array([0, 0, 0]), np.array([Lx, Ly, Lz])], [10, 10, 10], cell_type=dolfinx.mesh.CellType.hexahedron)
        fdim = domain.topology.dim - 1

        V = dolfinx.fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

        # ----------------------        
        # Locate boundary facets
        # ----------------------

        left_facets = dolfinx.mesh.locate_entities_boundary(
            domain, fdim, lambda x: np.isclose(x[0], Lx)
        )
        
        right_facets = dolfinx.mesh.locate_entities_boundary(
            domain, fdim, lambda x: np.isclose(x[0], 0.0)
        )
        
        # ---------------------
        # Mesh marker, meshtags
        # ---------------------

        num_facets = domain.topology.index_map(2).size_local
        markers = np.zeros(num_facets, dtype=np.int32)

        left = 1
        right = 2

        markers[left_facets] = left
        markers[right_facets] = right

        facet_marker = dolfinx.mesh.meshtags(domain, fdim, np.arange(num_facets, dtype=np.int32), markers)

        # ----------------------------------------------------------------
        # Locate dofs for Dirichlet boundary conditions from facet_markers
        # ----------------------------------------------------------------

        # dofs_D_left = dolfinx.fem.locate_dofs_topological(
        #     V, facet_marker.dim, facet_marker.find(left)
        # )

        # dofs_D_right = dolfinx.fem.locate_dofs_topological(
        #     V, facet_marker.dim, facet_marker.find(right)
        # )

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        def boundary_bottom(x):
            return np.isclose(x[2], 0)
        
        def boundary_top(x):
            return np.isclose(x[2], Lz)
        
        boundary_facets_bottom = dolfinx.mesh.locate_entities_boundary(domain, fdim, boundary_bottom)

        boundary_facets_top = dolfinx.mesh.locate_entities_boundary(domain, fdim, boundary_top)

        u_D = np.array([0, 0, 0], dtype=dolfinx.default_scalar_type)

        bc_bottom = dolfinx.fem.dirichletbc(u_D, dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets_bottom), V)

        bc_top = dolfinx.fem.dirichletbc(u_D, dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets_top), V)

        bc = [bc_bottom, bc_top]

        normal = ufl.FacetNormal(domain)

        pressure = dolfinx.fem.Constant(
            domain,
            dolfinx.default_scalar_type(0.01)
        )

        p0 = dolfinx.fem.Constant(
            domain,
            dolfinx.default_scalar_type(0.0)
        )
        
        T = pressure * normal

        ds = ufl.Measure(
            "ds",
            domain=domain,
            subdomain_data=facet_marker
        )

        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u):
            return 100 * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * 100 * epsilon(u)

        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        f = dolfinx.fem.Constant(domain, dolfinx.default_scalar_type((0, 0, 0)))
        
        a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
        L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds(left) + ufl.dot(p0*normal, v) * ds(right)

        problem = dolfinx.fem.petsc.LinearProblem(a, L, bcs=bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        uh = problem.solve()

        pyvista.start_xvfb()

        p = pyvista.Plotter()
        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(V)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        grid["u"] = uh.x.array.reshape((geometry.shape[0], 3))
        actor_0 = p.add_mesh(grid, style="wireframe", color="k")
        warped = grid.warp_by_vector("u", factor=1.5)
        actor_1 = p.add_mesh(warped, show_edges=True)
        p.show_axes()
        if not pyvista.OFF_SCREEN:
            p.show()
        else:
            figure_as_array = p.screenshot("displacement.png")

        with dolfinx.io.XDMFFile(domain.comm, "test_problem_3d_rect_u.xdmf", "w")as xdmf:
            xdmf.write_mesh(domain)
            uh.name = "u"
            xdmf.write_function(uh)

        gdim = domain.geometry.dim
        mesh_nodes = len(domain.geometry.x)

        stress_space = dolfinx.fem.functionspace(domain, ("Lagrange", 1, (gdim, gdim)))
        stress_expression = dolfinx.fem.Expression(sigma(uh), stress_space.element.interpolation_points())
        stress_values = dolfinx.fem.Function(stress_space)
        stress_values.interpolate(stress_expression)

        output_file = "test_problem_3d_rect_sigma.xdmf"
        with dolfinx.io.XDMFFile(domain.comm, output_file, "w") as xdmf:
            xdmf.write_mesh(domain)
            stress_values.name = "stress"
            xdmf.write_function(stress_values)
        print(f"Stress tensor written to {output_file}")


