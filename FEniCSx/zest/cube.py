import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from finiteElementSettings import FiniteElementSetting

class Cube(FiniteElementSetting):
    def volume(self, edge):
        """
        Compute the volume of the domain.
        """
        self.edge = edge
        
        # Integration measures the domain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)

        # Compute volume
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))

        # Print results
        print("\n")
        print("-----------")
        print("Mesh volume")
        print("-----------")
        print(f"Total volume, from mesh = {volume:.3f} m3")

        analytical_volume = self.edge**3

        print("\n")
        print("-----------------")
        print("Analytical volume")
        print("-----------------")
        print(f"Total volume, analytical = {analytical_volume:.3f} m3")
        print("\n")
        print(f"Relative error (%) = {abs(volume - analytical_volume) / analytical_volume * 100:.3e}")
        print("\n")

        return volume

    def parameters(self, rho, g, E, poisson, k, lhr, alpha):
        """
        Set material and physical parameters.

        Parameters:
        rho (float): Density.
        g (float): Gravitational acceleration.
        E (float): Young/Elastic Modulus.
        poisson (float): Poisson's ratio.
        k (float): Thermal conductivity.
        lhr (float): Linear heat rate.
        q_third (float): Volumetric power.
        alpha (float): Thermal expansion coefficient.
        """

        self.rho = rho
        self.g = g
        self.E = E
        self.poisson = poisson
        self.lmbda =  E * poisson / ((1 + poisson) * (1 - 2 * poisson))
        self.mu = E / (2 * (1 + poisson))
        self.bulk_modulus = E / (3 * (1 - 2 * poisson))
        self.k = k
        self.lhr = lhr
        self.q_third = lhr / (self.edge**2)
        self.alpha = alpha
        self.T_ref = 300.0

        print("----------")
        print("Parameters")
        print("----------")
        print(f"Edge: {self.edge} m")
        print(f"Material density: {self.rho} kg/m3")
        print(f"Gravitational acceleration: {self.g} m/s2")
        print(f"Elastic Modulus: {self.E/1e9} GPa")
        print(f"Poisson's ratio: {self.poisson}")
        print(f"Bulk modulus: {self.bulk_modulus/1e9:.2f} GPa")
        print(f"First Lamé parameter: {self.lmbda/1e9:.2f} GPa")
        print(f"Second Lamé parameter: {self.mu/1e9:.2f} GPa")
        print(f"Thermal conductivity: {self.k} W/mK")
        print(f"Linear heat rate: {self.lhr/1000} kW/m")
        print(f"Volumetric power: {self.q_third/1000:.2f} kW/m3")
        print(f"Thermal expansion coefficient: {self.alpha} 1/K")
        print(f"Reference temperature: {self.T_ref:.2f} K")
        print("\n")

    def boundary_conditions(self, pressure):
        """
        Define and set the boundary conditions.
        """

        print("-------------------")
        print("Boundary conditions")
        print("-------------------")
        print(f"Pressure = {pressure*1e-6:.2f} MPa")
        print("\n")

        self.bottom_label = 1
        self.lateral_1_label = 2
        self.lateral_2_label = 3
        self.lateral_3_label = 4
        self.lateral_4_label = 5
        self.top_label= 6

        # Verify face tags
        # for label in [self.bottom_label, self.lateral_1_label, self.lateral_2_label, self.lateral_3_label, self.lateral_4_label, self.top_label]:
        #     entities = self.ft.find(label)
        #     print(f"Tag {label}: {entities}")

        # ----
        # Dofs
        # ----

        self.dofs_bottom    = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.bottom_label))
        self.dofs_lateral_1 = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.lateral_1_label))
        self.dofs_lateral_2 = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.lateral_2_label))
        self.dofs_lateral_3 = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.lateral_3_label))
        self.dofs_lateral_4 = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.lateral_4_label))
        self.dofs_top       = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.top_label))
        
        def lateral_boundary_x_0(x):
            return np.logical_and(np.isclose(x[0], 0.0), x[1] < self.edge, x[2] < self.edge)

        def lateral_boundary_x_L(x):
            return np.logical_and(np.isclose(x[0], self.edge), x[1] < self.edge, x[2] < self.edge)
        
        def lateral_boundary_3(x):
            return np.logical_and(np.isclose(x[1], 0.0), x[0] < self.edge, x[2] < self.edge)

        def lateral_boundary_4(x):
            return np.logical_and(np.isclose(x[1], self.edge), x[0] < self.edge, x[2] < self.edge)

        def bottom_boundary(x):
            return np.logical_and(np.isclose(x[2], 0.0), x[0] < self.edge, x[1] < self.edge)

        def top_boundary(x):
            return np.logical_and(np.isclose(x[2], self.edge), x[0] < self.edge, x[1] < self.edge)
        
        # Boundary facets
        boundary_facets_lateral_x_0 = dolfinx.mesh.locate_entities_boundary(self.mesh, self.fdim, lateral_boundary_x_0)
        boundary_facets_lateral_x_L = dolfinx.mesh.locate_entities_boundary(self.mesh, self.fdim, lateral_boundary_x_L)
        boundary_facets_lateral_3 = dolfinx.mesh.locate_entities_boundary(self.mesh, self.fdim, lateral_boundary_3)
        boundary_facets_lateral_4 = dolfinx.mesh.locate_entities_boundary(self.mesh, self.fdim, lateral_boundary_4)
        boundary_facets_bottom    = dolfinx.mesh.locate_entities_boundary(self.mesh, self.fdim, bottom_boundary)
        boundary_facets_top       = dolfinx.mesh.locate_entities_boundary(self.mesh, self.fdim, top_boundary)

        # Boundary DoFs
        self.dofs_lateral_x_0_boundary = dolfinx.fem.locate_dofs_topological(self.V_m, self.fdim, boundary_facets_lateral_x_0)
        self.dofs_lateral_x_L_boundary = dolfinx.fem.locate_dofs_topological(self.V_m, self.fdim, boundary_facets_lateral_x_L)
        self.dofs_bottom_boundary = dolfinx.fem.locate_dofs_topological(self.V_m, self.fdim, boundary_facets_bottom)

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Mechanical boundary conditions
        # self.bc_mech_top = dolfinx.fem.dirichletbc(dolfinx.default_scalar_type(0), self.dofs_top, self.V_m.sub(2))
        # self.bc_mech_top = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_top, self.V_m)
        # self.bc_mech_lateral = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m)
        self.bc_mech_bottom = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom_boundary, self.V_m)

        self.bc_mech = []

        # Thermal boundary conditions
        bc_thermal_x_0 = dolfinx.fem.dirichletbc(self.T_ref, self.dofs_lateral_x_0_boundary, self.V_t)
        bc_thermal_x_L = dolfinx.fem.dirichletbc(self.T_ref, self.dofs_lateral_x_L_boundary, self.V_t)
        self.bc_thermal = [bc_thermal_x_0, bc_thermal_x_L]
        # self.bc_thermal_lateral_1 = dolfinx.fem.dirichletbc(400.0, self.dofs_lateral_1, self.V_t)
        # self.bc_thermal_lateral_3 = dolfinx.fem.dirichletbc(400.0, self.dofs_lateral_3, self.V_t)
        # self.bc_thermal = [self.bc_thermal_lateral_1, self.bc_thermal_lateral_3]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)
        self.traction = dolfinx.fem.Constant(self.mesh, pressure) * self.normal
        # self.traction_bottom = dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type((0.0, 0.0, pressure)))
        # self.traction_up = dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type((0.0, 0.0, - pressure)))
        # self.traction_y_0 = dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type((pressure, 0.0, 0.0)))

        print("Traction")
        print(self.traction)

    def assemble_and_solve(self):
        """
        Assemble and solve the system for the linear thermo-elastic problem.
        """

        # Mesure "ds"
        self.ds = ufl.Measure("ds", domain=self.mesh)
        ds_1 = self.ds(self.lateral_1_label)
        ds_2 = self.ds(self.lateral_2_label)
        ds_3 = self.ds(self.lateral_3_label)
        ds_4 = self.ds(self.lateral_4_label)

        print("---------------")
        print("Thermal problem")
        print("---------------")

        self.a_thermal = ufl.inner(self.k * ufl.grad(self.u_t), ufl.grad(self.v_t)) * ufl.dx
        self.L_thermal = dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type(self.q_third)) * self.v_t * ufl.dx
        
        petsc_options_sym = {
            "ksp_type": "preonly",
            "pc_type": "lu", # LU factorization, a direct solver method (similar to Gaussian elimination)
            "pc_factor_mat_solver_type": "mumps",
        }

        # For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES, is better.
        petsc_options_non_sym = {
            "ksp_type": "gmres",
            "ksp_rtol": 1e-6,
            "ksp_atol": 1e-10,
            "ksp_max_it": 1000,
            "pc_type": "none",
        }

        jit_options = {
            "timeout": 20,
        }

        self.assemble_and_print_matrix(a = self.a_thermal, L = self.L_thermal, filename = 'thermal_stiffness_matrix.txt')

        self.thermal_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_thermal,
            self.L_thermal,
            bcs = self.bc_thermal,
            petsc_options=petsc_options_sym,
            jit_options=jit_options,
        )

        self.temperature = self.thermal_problem.solve()
        self.temperature.x.scatter_forward() # to propagate the data among the processors
        self.temperature.name = "Temperature"

        print(f"Analytical thermal jump: {self.q_third * (self.edge / 2)**2 / (2 *self.k):.2f} K")
        print("\n")
        
        print("-------------------")
        print("Mechcanical problem")
        print("-------------------")

        def epsilon(u):
            """
            Strain-displacement compatibility relations
            """
            return ufl.sym(ufl.grad(u))
        
        def sigma(u):
            """
            Stress-strain constitutive relation.
            """
            return self.lmbda * ufl.tr(epsilon(u)) * ufl.Identity(len(u)) + 2.0 * self.mu * epsilon(u)

        # Variational forms
        self.a_mech = ufl.inner(sigma(self.u_m), epsilon(self.v_m)) * ufl.dx
        self.L_body_forces = ufl.dot(dolfinx.fem.Constant(self.mesh, dolfinx.default_scalar_type((0.0, 0.0, - self.rho * self.g))), self.v_m) * ufl.dx
        # self.L_tractions = ufl.dot(self.traction, self.v_m) * ds_1 + ufl.dot(self.traction, self.v_m) * ds_2 + ufl.dot(self.traction, self.v_m) * ds_3 + ufl.dot(self.traction, self.v_m) * ds_4
        self.L_tractions = ufl.dot(self.traction, self.v_m) * self.ds
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        # dx = ufl.Measure("dx", domain=self.mesh)
        # T_avg = dolfinx.fem.assemble_scalar(dolfinx.fem.form(self.temperature * dx)) / dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx))
        self.L_thermal = (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.div(self.v_m) * ufl.dx
        self.L_mech = self.L_body_forces + self.L_thermal + self.L_tractions

        self.assemble_and_print_matrix(a = self.a_mech, L = self.L_mech, filename = 'mech_stiffness_matrix.txt')

        # Linear problem definition
        self.mech_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_mech,
            self.L_mech,
            bcs = self.bc_mech,
            petsc_options = petsc_options_sym,
            jit_options = jit_options,
        )

        # Non-linear problem definition
        # F_mech = self.a_mech - self.L_body_forces
        # mech_problem = dolfinx.fem.petsc.NonlinearProblem(F_mech, self.u_m, bcs=self.bc_mech)
        # solver = dolfinx.nls.petsc.NewtonSolver(MPI.COMM_WORLD, mech_problem)
        # solver.convergence_criterion = "incremental"
        # solver.rtol = 1e-6
        # solver.report = True
        # solver.solve(self.u_m)

        # self.displacement = self.u_m
        self.displacement = self.mech_problem.solve()
        self.displacement.name = "Displacement"
        self.displacement.x.scatter_forward() # to propagate the data among the processors
        self.stress = sigma(self.displacement)
        self.strain = epsilon(self.displacement)

        # Solver options
        solver = self.thermal_problem.solver
        viewer = PETSc.Viewer().createASCII("log.thermal_solver")
        solver.view(viewer)

        solver = self.mech_problem.solver
        viewer = PETSc.Viewer().createASCII("log.mech_solver")
        solver.view(viewer)

    def volumetric_stress(self):
        """
        Compute the volumetric component of the stress field.
        The volumetric stress is the trace of the stress tensor divided by the geometric dimension.
        It is related to the change in volume of the material, and the hydrostatic pressure.
        """

        I = ufl.Identity(self.gdim)
        volumetric_stress = ufl.tr(self.stress) * I / self.gdim

        return volumetric_stress

    def deviatoric_stress(self):
        """
        Compute the deviatoric component of the stress field.
        The deviatoric stress is the stress tensor minus the volumetric stress.
        It is related to the shape change of the material, and the shear stress.
        """
        
        I = ufl.Identity(self.gdim)
        deviatoric_stress = self.stress - ufl.tr(self.stress) * I / self.gdim

        return deviatoric_stress

    def get_results(self):
        """
        Get the computed displacement, stress, and strain fields.

        Returns:
        displacement (Function): The displacement field.
        stress (ufl.Expr): The stress field.
        strain (ufl.Expr): The strain field.
        """

        return self.displacement, self.stress, self.strain, self.temperature
