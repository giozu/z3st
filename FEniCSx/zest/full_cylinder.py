import dolfinx.fem.petsc
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from finiteElementSettings import FiniteElementSetting

class FullCylinder(FiniteElementSetting):
    def volume(self, outer_radius, inner_radius, height):
        """
        Compute the volume of the domain.
        """
        self.outer_radius = outer_radius
        self.inner_radius = inner_radius
        self.height = height
        self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
        
        # Integration measures the domain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)

        # Compute volume
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))

        # Print results
        print("\n")
        print("-----------")
        print("Mesh volume")
        print("-----------")
        print(f"Total volume, from mesh = {volume*1e9:.3e} mm3")

        analytical_volume = self.area * self.height

        print("\n")
        print("-----------------")
        print("Analytical volume")
        print("-----------------")
        print(f"Total volume, analytical = {analytical_volume*1e9:.3e} mm3")
        print("\n")
        print(f"Relative error (%) = {abs(volume - analytical_volume) / analytical_volume * 100:.3e}")
        print("\n")

    def parameters(self, rho, g, E, poisson, k, lhr, alpha, T_ref, cp):
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
        T_ref (float): Reference temperature.
        cp (float): Specific heat capacity.
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
        self.q_third = lhr / self.area
        self.alpha = alpha
        self.T_ref = T_ref
        self.cp = cp
       
        print("----------")
        print("Parameters")
        print("----------")
        print(f"Inner radius: {self.inner_radius:.3e} m")
        print(f"Outer radius: {self.outer_radius:.3e} m")
        print(f"Height: {self.height:.3e} m")
        print(f"Material density: {self.rho:.2f} kg/m3")
        print(f"Gravitational acceleration: {self.g:.2f} m/s2")
        print(f"Elastic modulus: {self.E/1e9:.2f} GPa")
        print(f"Poisson's ratio: {self.poisson:.2f}")
        print(f"Bulk modulus: {self.bulk_modulus/1e9:.2f} GPa")
        print(f"First Lamé parameter: {self.lmbda/1e9:.2f} GPa")
        print(f"Second Lamé parameter: {self.mu/1e9:.2f} GPa")
        print(f"Thermal conductivity: {self.k:.2f} W/mK")
        print(f"Specific heat capacity: {self.cp:.2f} J/kg-K")
        print(f"Linear heat rate: {self.lhr/1000:.2f} kW/m")
        print(f"Volumetric power: {self.q_third/1000000:.2f} MW/m3")
        print(f"Thermal expansion coefficient: {self.alpha:.2e} 1/K")
        print(f"Reference temperature: {self.T_ref:.2f} K")
        print("\n")

    def boundary_conditions(self, pressure, temperature):
        """
        Define and set the boundary conditions.
        """

        print("-------------------")
        print("Boundary conditions")
        print("-------------------")
        print(f"Pressure = {pressure*1e-6:.2f} MPa")
        print(f"Temperature = {temperature:.2f} K")
        print("\n")

        self.bottom_label = 1
        self.lateral_label = 2
        self.top_label= 3

        # ----
        # Dofs
        # ----

        self.dofs_bottom  = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.bottom_label))
        self.dofs_lateral = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.lateral_label))
        self.dofs_top     = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.top_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Mechanical boundary conditions
        # self.bc_mech_top = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_top, self.V_m)
        self.bc_mech_bottom = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m)
        # self.bc_mech_lateral = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m)

        self.bc_mech = []

        # Thermal boundary conditions
        self.T_bc = temperature
        self.bc_thermal_lateral = dolfinx.fem.dirichletbc(temperature, self.dofs_lateral, self.V_t)
        self.bc_thermal = [self.bc_thermal_lateral]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)
        self.pressure = pressure
        self.traction = dolfinx.fem.Constant(self.mesh, pressure) * self.normal

    def assemble_and_solve(self):
        """
        Assemble and solve the system for the linear thermo-elastic problem.
        """

        # Direct solver + MUMPS
        petsc_options_lu = {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "ksp_rtol": 1e-8,
            "ksp_max_it": 1000,
            "ksp_monitor": None
        }

        # CG-solver + jacobi
        petsc_options_cg = {
            "ksp_type": "cg",
            "pc_type": "jacobi",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-10,
            "ksp_max_it": 1000,
            "ksp_monitor": None
        }

        # GMRES solver + pre-conditioner AMG (Algebraic Multigrid)
        petsc_options_gmres = {
            "ksp_type": "gmres",
            "pc_type": "hypre",
            "pc_hypre_type": "boomeramg",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-10,
            "ksp_max_it": 1000,
            "pc_type": None,
        }

        # GMRES iterative solver + ILU pre-conditioner
        petsc_options_ilu = {
            "ksp_type": "gmres",
            "pc_type": "ilu",
            "ksp_rtol": 1e-8,
            "ksp_max_it": 500,
            "ksp_monitor": None
        }

        jit_options = {
            "timeout": 20,
        }

        # Integration mesures "ds", "dx"
        metadata = {"quadrature_degree": 4}

        self.ds = ufl.Measure("ds", domain=self.mesh, metadata=metadata)
        # ds_lateral = self.ds(self.lateral_label)
        ds = self.ds
        dx = ufl.Measure("dx", domain=self.mesh, metadata=metadata)


        print("---------------")
        print("Thermal problem")
        print("---------------")

        characteristic_time = self.rho * self.cp / self.k * self.area
        print(f"Characteristic time: {characteristic_time:.3e} s")

        self.a_thermal = ufl.dot(self.k * ufl.grad(self.u_t), ufl.grad(self.v_t)) * ufl.dx
        self.L_thermal = dolfinx.fem.Constant(self.mesh, self.q_third) * self.v_t * ufl.dx
        
        self.assemble_and_print_matrix(a = self.a_thermal, L = self.L_thermal, filename = 'thermal_stiffness_matrix.txt')

        self.thermal_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_thermal,
            self.L_thermal,
            bcs = self.bc_thermal,
            petsc_options=petsc_options_cg,
            jit_options=jit_options,
        )

        self.temperature = self.thermal_problem.solve()
        self.temperature.x.scatter_forward() # to propagate the data among the processors
        self.temperature.name = "temperature"

        thermal_resistance = 4.0 * np.pi * self.k
        print(f"Analytical thermal jump: {self.lhr / thermal_resistance:.2f} K")
        print("\n")

        print("------------------")
        print("Mechanical problem")
        print("------------------")

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
        self.L_body_forces = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, - self.rho * self.g)), self.v_m) * ufl.dx
        self.L_tractions = ufl.dot(self.traction, self.v_m) * ds
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        self.L_thermal = (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.div(self.v_m) * ufl.dx
        self.L_mech = self.L_body_forces + self.L_thermal + self.L_tractions

        self.assemble_and_print_matrix(a = self.a_mech, L = self.L_mech, filename = 'mech_stiffness_matrix.txt')

        # Linear problem definition
        self.mech_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_mech,
            self.L_mech,
            bcs = self.bc_mech,
            petsc_options = petsc_options_cg,
            jit_options = jit_options,
        )

        self.displacement = self.mech_problem.solve()
        self.displacement.name = "displacement"
        self.displacement.x.scatter_forward() # to propagate the data among the processors
        self.stress = self.lmbda * ufl.tr(epsilon(self.displacement)) * ufl.Identity(len(self.displacement)) + 2.0 * self.mu * epsilon(self.displacement) - (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))
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
