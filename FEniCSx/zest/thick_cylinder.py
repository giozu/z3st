import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from finiteElementSettings import FiniteElementSetting

class ThickCylinder(FiniteElementSetting):
    def volume(self, outer_radius, inner_radius, height):
        """
        Compute the volume of the domain.
        """

        self.outer_radius = outer_radius
        self.inner_radius = inner_radius
        self.height = height
        
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

        analytical_volume = np.pi * (self.outer_radius**2 - self.inner_radius**2) * self.height

        print("\n")
        print("-----------------")
        print("Analytical volume")
        print("-----------------")
        print(f"Total volume, analytical = {analytical_volume:.3f} m3")
        print("\n")
        print("Relative error (%) = ", abs(volume - analytical_volume) / analytical_volume * 100)
        print("\n")

        return volume

    def parameters(self, rho, g, lmbda, mu, k, lhr, alpha):
        """
        Set material and physical parameters.

        Parameters:
        rho (float): Density.
        g (float): Gravitational acceleration.
        lmbda (float): First Lamé parameter.
        mu (float): Second Lamé parameter.
        k (float): Thermal conductivity.
        lhr (float): Linear heat rate.
        q_third (float): Volumetric power.
        alpha (float): Thermal expansion coefficient.
        """
        self.rho = rho
        self.g = g
        self.lmbda = lmbda
        self.mu = mu
        self.k = k
        self.lhr = lhr
        self.q_third = lhr / (np.pi * (self.outer_radius**2 - self.inner_radius**2))
        self.alpha = alpha

        print("\n")
        print("----------")
        print("Parameters")
        print("----------")
        print(f"Inner radius: {self.inner_radius} m")
        print(f"Outer radius: {self.outer_radius} m")
        print(f"Height: {self.height} m")
        print(f"Material density: {self.rho} kg/m3")
        print(f"Gravitational acceleration: {self.g} m/s2")
        print(f"First Lamé parameter: {self.lmbda/1e9:.2f} GPa")
        print(f"Second Lamé parameter: {self.mu/1e9:.2f} GPa")
        print(f"Thermal conductivity: {self.k} W/mK")
        print(f"Linear heat rate: {self.lhr/1000} kW/m")
        print(f"Volumetric power: {self.q_third/1000:.2f} kW/m3")
        print(f"Thermal expansion coefficient: {self.alpha} 1/K")
        print("\n")

    def boundary_conditions(self, inner_pressure, outer_pressure):
        """
        Define and set the boundary conditions.
        """

        print("\n")
        print("-------------------")
        print("Boundary conditions")
        print("-------------------")
        print(f"Inner pressure = {inner_pressure*1e-6} MPa")
        print(f"Outer pressure = {outer_pressure*1e-6} MPa")
        print("\n")

        self.bottom_label = 1
        self.inner_lateral_label = 2
        self.outer_lateral_label = 3
        self.top_label= 4

        # ----
        # Dofs
        # ----

        self.dofs_bottom        = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.bottom_label))
        self.dofs_inner_lateral = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.inner_lateral_label))
        self.dofs_top           = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.top_label))
        self.dofs_outer_lateral = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.outer_lateral_label))

        self.dofs_outer_lateral_t = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.outer_lateral_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Mechanical boundary conditions
        self.bc_mech_top = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_top, self.V_m)
        self.bc_mech_bottom = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m)
        self.bc_mech = [self.bc_mech_bottom, self.bc_mech_top]

        # Thermal boundary conditions
        self.bc_thermal_outer = dolfinx.fem.dirichletbc(400.0, self.dofs_outer_lateral_t, self.V_t)
        self.bc_thermal = [self.bc_thermal_outer]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)

        self.inner_pressure = dolfinx.fem.Constant(self.mesh, inner_pressure)
        self.inner_traction = self.inner_pressure * self.normal

        self.outer_pressure = dolfinx.fem.Constant(self.mesh, outer_pressure)
        self.outer_traction = self.outer_pressure * self.normal

    def assemble_and_solve(self):
        """
        Assemble and solve the system for the linear thermo-elastic problem.
        """

        # ---------------
        # Thermal problem
        # ---------------
        self.a_thermal = ufl.inner(self.k * ufl.grad(self.u_t), ufl.grad(self.v_t)) * ufl.dx
        self.L_thermal = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.q_third)) * self.v_t * ufl.dx

        petsc_options = {
            "ksp_type": "preonly",
            "pc_type": "lu", # LU factorization, a direct solver method (similar to Gaussian elimination)
            "pc_factor_mat_solver_type": "mumps",
        }

        # For non-symmetric problems, a Krylov solver for non-symmetric systems, such as GMRES, is better.
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

        self.assemble_and_print_matrix(a = self.a_thermal, L = self.L_thermal, filename = 'thermal_stiffness_matrix.txt')

        self.thermal_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_thermal,
            self.L_thermal,
            bcs = self.bc_thermal,
            petsc_options=petsc_options,
            jit_options=jit_options,
        )

        self.temperature = self.thermal_problem.solve()

        print(f"Analytical thermal jump: {self.lhr / ( 4.0 * np.pi * self.k) * (1. - self.inner_radius**2/(self.outer_radius**2 - self.inner_radius**2) * 2. * np.log(self.outer_radius/self.inner_radius)):.2f} K")
        print("\n")

        # -----------------
        # Mechnical problem
        # -----------------

        # Mesure "ds"
        self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft)
        ds_i = self.ds(self.inner_lateral_label)
        ds_o = self.ds(self.outer_lateral_label)

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
        self.a_mech = ufl.inner(sigma(self.u_m), epsilon(self.v_m))  * ufl.dx
        self.L_body_forces = ufl.inner(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, - self.rho * self.g)), self.v_m) * ufl.dx
        self.L_tractions = ufl.inner(self.inner_traction, self.v_m) * ds_i + ufl.inner(self.outer_traction, self.v_m) * ds_o
        # self.L_thermal = (3 * self.lmbda + 2 * self.mu) * self.alpha * self.temperature * ufl.div(self.v_m) * ufl.dx
        self.L_mech = self.L_body_forces + self.L_tractions # + self.L_thermal

        self.assemble_and_print_matrix(a = self.a_mech, L = self.L_mech, filename = 'mech_stiffness_matrix.txt')

        # Linear problem definition
        self.mech_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_mech,
            self.L_mech,
            bcs = self.bc_mech,
            petsc_options = petsc_options,
            jit_options = jit_options,
        )

        self.displacement = self.mech_problem.solve()
        # self.displacement.x.scatter_forward() # scatter_forward is necessary for dolfinx ?
        self.stress_field = sigma(self.displacement)
        self.strain_field = epsilon(self.displacement)

        # Solver options
        solver = self.thermal_problem.solver
        viewer = PETSc.Viewer().createASCII("log.thermal_solver")
        solver.view(viewer)
        # solver_output = open("log.thermal_solver", "r")
        # for line in solver_output.readlines():
        #     print(line)

        solver = self.mech_problem.solver
        viewer = PETSc.Viewer().createASCII("log.mech_solver")
        solver.view(viewer)
        # solver_output = open("log.mech_solver", "r")
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

        # Save the displacement field to a VTK file
        vtk = dolfinx.io.VTKFile(self.mesh.comm, "linear_elasticity.pvd", "w")
        vtk.write_function(self.displacement)
        vtk.close()

        return self.displacement, self.stress_field, self.strain_field, self.temperature
