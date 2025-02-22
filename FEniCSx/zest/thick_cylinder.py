import dolfinx.fem.petsc
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
        print(f"Thermal expansion coefficient: {self.alpha:.2e} 1/K")
        print(f"Reference temperature: {self.T_ref:.2f} K")
        print("\n")

    def boundary_conditions(self, inner_pressure, outer_pressure, temperature):
        """
        Define and set the boundary conditions.
        """

        print("-------------------")
        print("Boundary conditions")
        print("-------------------")
        print(f"Inner pressure = {inner_pressure*1e-6:.2f} MPa")
        print(f"Outer pressure = {outer_pressure*1e-6:.2f} MPa")
        print(f"Temperature = {temperature:.2f} K")
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

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Mechanical boundary conditions
        self.bc_mech_bottom = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m)
        self.bc_mech_top = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_top, self.V_m)
        self.bc_mech_inner_lateral = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_inner_lateral, self.V_m)
        self.bc_mech_outer_lateral = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_outer_lateral, self.V_m)

        self.bc_mech = []

        # Thermal boundary conditions
        self.T_bc = temperature
        self.bc_thermal_outer = dolfinx.fem.dirichletbc(temperature, self.dofs_outer_lateral, self.V_t)
        self.bc_thermal = [self.bc_thermal_outer]

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)
        self.inner_traction = dolfinx.fem.Constant(self.mesh, inner_pressure) * self.normal
        self.outer_traction = dolfinx.fem.Constant(self.mesh, outer_pressure) * self.normal

        self.heat_flux = self.lhr / (2.0 * np.pi * self.inner_radius) # W/m2

    def assemble_and_solve(self):
        """
        Assemble and solve the system for the linear thermo-elastic problem.
        """
        # Integration mesures "ds", "dx"
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure("dx", domain=self.mesh, metadata=metadata)
        self.ds = ufl.Measure("ds", domain=self.mesh, metadata=metadata, subdomain_data=self.ft)
        self.ds_i = self.ds(self.inner_lateral_label)
        self.ds_o = self.ds(self.outer_lateral_label)

        self.solver_options()
       
        self.thermal_solver()
        self.mechanical_solver()
        
    def solver_options(self):
        # Direct solver + MUMPS
        self.petsc_options_lu = {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "ksp_rtol": 1e-8,
            "ksp_max_it": 1000,
            "ksp_monitor": None
        }

        # CG-solver + jacobi
        self.petsc_options_cg = {
            "ksp_type": "cg",
            "pc_type": "jacobi",
            "ksp_rtol": 1e-8,
            "ksp_atol": 1e-10,
            "ksp_max_it": 2000,
            "ksp_monitor": None
        }

        # GMRES solver + pre-conditioner AMG (Algebraic Multigrid)
        self.petsc_options_gmres = {
            "ksp_type": "gmres",
            "pc_type": "hypre",
            "pc_hypre_type": "boomeramg",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-10,
            "ksp_max_it": 1000,
            "pc_type": None,
        }

        # GMRES iterative solver + ILU pre-conditioner
        self.petsc_options_ilu = {
            "ksp_type": "gmres",
            "pc_type": "ilu",
            "ksp_rtol": 1e-8,
            "ksp_atol": 1e-10,
            "ksp_max_it": 2000,
            "ksp_monitor": None
        }

        self.jit_options = {
            "timeout": 20,
        }

    def thermal_solver(self):

        print("--------------")
        print("Thermal solver")
        print("--------------")

        characteristic_time = self.rho * self.cp / self.k * self.area
        print(f"Characteristic time: {characteristic_time:.3e} s")

        self.a_thermal = ufl.dot(self.k * ufl.grad(self.u_t), ufl.grad(self.v_t)) * self.dx
        self.L_thermal = self.heat_flux * self.v_t * self.ds_i

        integral_value = dolfinx.fem.assemble_scalar(dolfinx.fem.form(self.heat_flux * self.ds_i))
        print(f"Total power from the integrated value of the heat flux on the inner wall: {integral_value:.3f} W")
        print(f"Total power from the analytical value: {self.lhr * self.height:.3f} W")

        # self.assemble_and_print_matrix(a = self.a_thermal, L = self.L_thermal, filename = 'thermal_stiffness_matrix.txt')

        self.thermal_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_thermal,
            self.L_thermal,
            bcs = self.bc_thermal,
            petsc_options=self.petsc_options_lu,
            jit_options=self.jit_options,
        )

        self.temperature = self.thermal_problem.solve()
        self.temperature.x.scatter_forward() # to propagate the data among the processors
        self.temperature.name = "temperature"

        # print(f"Analytical thermal jump: {self.lhr / ( 4.0 * np.pi * self.k) * (1. - self.inner_radius**2/(self.outer_radius**2 - self.inner_radius**2) * 2. * np.log(self.outer_radius/self.inner_radius)):.2f} K")
        thermal_resistance = np.log(self.outer_radius/self.inner_radius) / (2.0*np.pi*self.k)
        print(f"Analytical thermal jump: {self.lhr * thermal_resistance:.2f} K")
        print(f"Min temperature from thermal solver: {np.min(self.temperature.x.array):.2f} K")
        print(f"Max temperature from thermal solver: {np.max(self.temperature.x.array):.2f} K")
        print("\n")

    def mechanical_solver(self):

        print("-----------------")
        print("Mechanical solver")
        print("-----------------")

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
        self.a_mech = ufl.inner(sigma(self.u_m), epsilon(self.v_m)) * self.dx
        self.L_body_forces = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, - self.rho * self.g)), self.v_m) * self.dx
        self.L_tractions = ufl.dot(self.inner_traction, self.v_m) * self.ds_i + ufl.inner(self.outer_traction, self.v_m) * self.ds_o
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        self.L_thermal_exp = (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.div(self.v_m) * self.dx
        self.L_mech = self.L_body_forces + self.L_thermal_exp + self.L_tractions

        # self.assemble_and_print_matrix(a = self.a_mech, L = self.L_mech, filename = 'mech_stiffness_matrix.txt')

        # Linear problem definition
        self.mech_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_mech,
            self.L_mech,
            bcs = self.bc_mech,
            petsc_options = self.petsc_options_ilu,
            jit_options = self.jit_options,
        )

        self.displacement = self.mech_problem.solve()
        self.displacement.name = "displacement"
        self.displacement.x.scatter_forward() # to propagate the data among the processors
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        self.stress = sigma(self.displacement) - (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))
        self.strain = epsilon(self.displacement)

        max_displacement = np.max(np.abs(self.displacement.x.array))
        print(f"\nMax displacement: {max_displacement:.3e} m")

    def get_results(self):
        """
        Get the computed displacement, stress, and strain fields.

        Returns:
        displacement (Function): The displacement field.
        stress (ufl.Expr): The stress field.
        strain (ufl.Expr): The strain field.
        """

        return self.displacement, self.stress, self.strain, self.temperature

    def vonMises_equivalent(self):
        """
        Compute and visualize the Von Mises stress field.
        """

        # Compute deviatoric stress
        deviatoric_stress = self.stress - (1.0 / 3.0) * ufl.tr(self.stress) * ufl.Identity(len(self.displacement))
        
        # Compute Von Mises stress
        von_Mises = ufl.sqrt(3.0 / 2.0 * ufl.inner(deviatoric_stress, deviatoric_stress))

        # Define function space for Von Mises stress
        V_von_mises = dolfinx.fem.functionspace(self.mesh, ("DG", 0))

        # Interpolate the computed Von Mises stress
        von_mises_expr = dolfinx.fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
        von_mises_func = dolfinx.fem.Function(V_von_mises)
        von_mises_func.interpolate(von_mises_expr)

        # Start PyVista visualization
        import pyvista
        p = pyvista.Plotter()

        # Create the plotter and grid for visualization
        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(self.V_m)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
        
        # Attach vector values to grid and warp grid by vector
        grid["u"] = self.displacement.x.array.reshape((geometry.shape[0], 3))
        
        p.add_mesh(grid, style="wireframe", color="k") 
        
        # Attach Von Mises stress data to the grid
        warped = grid.warp_by_vector("u", factor=1)
        warped.cell_data["VonMises"] = von_mises_func.x.array
        warped.set_active_scalars("VonMises")

        # Plot the results
        p.add_mesh(warped)
        p.show_axes()

        # Show the plot (if not in off-screen mode)
        if not pyvista.OFF_SCREEN:
            p.show()
        else:
            p.screenshot("output/von_mises_stress.png")

        import os
        output_dir = "output"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with dolfinx.io.XDMFFile(self.mesh.comm, "output/von_mises_stress.xdmf", "w") as xdmf:
            xdmf.write_mesh(self.mesh)
            von_mises_func.name = "VonMises"
            xdmf.write_function(von_mises_func)

        return von_mises_func
