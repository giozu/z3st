import dolfinx.fem.petsc
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from finiteElementSettings import FiniteElementSetting

class FuelRod(FiniteElementSetting):
    def volume(self, fuel_outer_radius, cladding_inner_radius, cladding_outer_radius, height):
        """
        Compute the volume of the domain.
        """

        self.fuel_outer_radius = fuel_outer_radius
        self.cladding_inner_radius = cladding_inner_radius
        self.cladding_outer_radius = cladding_outer_radius

        self.height = height
        
        self.area_fuel = np.pi * self.fuel_outer_radius**2
        self.area_cladding = np.pi * (self.cladding_outer_radius**2 - self.cladding_inner_radius**2)

        # Integration measures the domain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
        self.dx_fuel = dx_subdomain(1)      # Marker ID for fuel region
        self.dx_cladding = dx_subdomain(2)  # Marker ID for cladding region

        # Compute volumes
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))
        volume_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_fuel))
        volume_cladding = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_cladding))

        # Print results
        print("\n")
        print("-----------")
        print("Mesh volume")
        print("-----------")
        print(f"Fuel volume, from mesh = {volume_fuel*1e9:.3e} mm3")
        print(f"Cladding volume, from mesh = {volume_cladding*1e9:.3e} mm3")

        print("\n")
        print("-----------------")
        print("Analytical volume")
        print("-----------------")
        
        analytical_volume = (self.area_fuel + self.area_cladding)*self.height

        print(f"Fuel volume, analytical = {self.area_fuel*self.height*1e9:.3e} mm3")
        print(f"Cladding volume, analytical = {self.area_cladding*self.height*1e9:.3e} mm3")
        print("\n")
        print(f"Relative error (%) = {abs(volume - (analytical_volume)) / analytical_volume * 100:.3e}")
        print("\n")

    def parameters(self, fuel, cladding, filling_gas, g, lhr, T_ref):
        """
        Set material and physical parameters for both the fuel and the cladding.

        Parameters:
        fuel (dict): Material parameters for the fuel with keys:
            "rho"   : density [kg/m³],
            "E"     : Young's modulus [Pa],
            "nu"    : Poisson's ratio,
            "k"     : thermal conductivity [W/mK],
            "alpha" : thermal expansion coefficient [1/K].
        cladding (dict): Material parameters for the cladding with the same keys.
        g (float): Gravitational acceleration [m/s²].
        lhr (float): Linear heat rate [W/m].
        T_ref (float): Reference temperature [K].
        cp (float): Specific heat capacity [J/kgK].
        """
        # Set global parameters
        self.g = g
        self.lhr = lhr
        self.T_ref = T_ref

        # Set fuel parameters and compute elastic properties for fuel
        self.fuel = {}
        self.fuel["rho"] = fuel["rho"]
        self.fuel["E"] = fuel["E"]
        self.fuel["nu"] = fuel["nu"]
        self.fuel["k"] = fuel["k"]
        self.fuel["alpha"] = fuel["alpha"]
        self.fuel["lmbda"] = fuel["E"] * fuel["nu"] / ((1 + fuel["nu"]) * (1 - 2 * fuel["nu"]))
        self.fuel["mu"] = fuel["E"] / (2 * (1 + fuel["nu"]))
        self.fuel["bulk_modulus"] = fuel["E"] / (3 * (1 - 2 * fuel["nu"]))
        self.fuel["cp"] = fuel["cp"]

        # Set cladding parameters and compute elastic properties for cladding
        self.cladding = {}
        self.cladding["rho"] = cladding["rho"]
        self.cladding["E"] = cladding["E"]
        self.cladding["nu"] = cladding["nu"]
        self.cladding["k"] = cladding["k"]
        self.cladding["alpha"] = cladding["alpha"]
        self.cladding["lmbda"] = cladding["E"] * cladding["nu"] / ((1 + cladding["nu"]) * (1 - 2 * cladding["nu"]))
        self.cladding["mu"] = cladding["E"] / (2 * (1 + cladding["nu"]))
        self.cladding["bulk_modulus"] = cladding["E"] / (3 * (1 - 2 * cladding["nu"]))
        self.cladding["cp"] = cladding["cp"]

        # Compute the volumetric heat generation (q_third) using the fuel area.
        self.q_third = lhr / self.area_fuel

        # Compute the (radial) gap thickness
        self.gap_thickness = self.cladding_inner_radius - self.fuel_outer_radius

        # Print fuel parameters
        print("---------------")
        print("Fuel parameters")
        print("---------------")
        print(f"Fuel inner radius: {self.fuel_outer_radius:.3e} m")
        print(f"Fuel density: {self.fuel['rho']:.2f} kg/m³")
        print(f"Fuel Elastic modulus: {self.fuel['E']/1e9:.2f} GPa")
        print(f"Fuel Poisson's ratio: {self.fuel['nu']:.2f}")
        print(f"Fuel Bulk modulus: {self.fuel['bulk_modulus']/1e9:.2f} GPa")
        print(f"Fuel Lamé parameters: {self.fuel['lmbda']/1e9:.2f} GPa, {self.fuel['mu']/1e9:.2f} GPa")
        print(f"Fuel Thermal conductivity: {self.fuel['k']:.2f} W/mK")
        print(f"Fuel Thermal expansion coefficient: {self.fuel['alpha']:.2e} 1/K")
        print(f"Fuel specific heat capacity: {self.fuel['cp']:.2f} J/kgK")
        print("-------------------")
        
        # Print cladding parameters
        print("Cladding parameters")
        print("-------------------")
        print(f"Cladding inner radius: {self.cladding_inner_radius:.3e} m")
        print(f"Cladding outer radius: {self.cladding_outer_radius:.3e} m")
        print(f"Gap thickness: {(self.cladding_inner_radius - self.fuel_outer_radius)*1e6:.3f} micron")
        print(f"Cladding thickness: {(self.cladding_outer_radius - self.cladding_inner_radius)*1e6:.3f} micron")
        print(f"Cladding density: {self.cladding['rho']:.2f} kg/m³")
        print(f"Cladding Elastic modulus: {self.cladding['E']/1e9:.2f} GPa")
        print(f"Cladding Poisson's ratio: {self.cladding['nu']:.2f}")
        print(f"Cladding Bulk modulus: {self.cladding['bulk_modulus']/1e9:.2f} GPa")
        print(f"Cladding Lamé parameters: {self.cladding['lmbda']/1e9:.2f} GPa, {self.cladding['mu']/1e9:.2f} GPa")
        print(f"Cladding Thermal conductivity: {self.cladding['k']:.2f} W/mK")
        print(f"Cladding Thermal expansion coefficient: {self.cladding['alpha']:.2e} 1/K")
        print(f"Cladding specific heat capacity: {self.cladding['cp']:.2f} J/kgK")
        print("-----------------")

        # Print cladding parameters
        print(f"Filling gas: {filling_gas['name']}")
        self.filling_gas = filling_gas
        print("-----------------")
        
        # Print global parameters
        print("Global parameters")
        print("-----------------")
        print(f"Gravitational acceleration: {g:.2f} m/s²")
        print(f"Linear heat rate: {lhr/1000:.2f} kW/m")
        print(f"Reference temperature: {T_ref:.2f} K")
        print("\n")

    def boundary_conditions(self, gap_pressure, coolant_pressure, temperature):
        """
        Define and set the boundary conditions.
        """

        print("-------------------")
        print("Boundary conditions")
        print("-------------------")
        print(f"Gap pressure = {gap_pressure*1e-6:.2f} MPa")
        print(f"Coolant pressure = {coolant_pressure*1e-6:.2f} MPa")
        print(f"Temperature = {temperature:.2f} K")
        print("\n")

        self.fuel_bottom_label = 1
        self.cladding_bottom_label = 2
        self.fuel_outer_label = 3
        self.cladding_outer_label = 4
        self.cladding_inner_label = 5
        self.cladding_top_label = 6
        self.fuel_top_label = 7

        # ----
        # Dofs
        # ----

        self.dofs_fuel_bottom     = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.fuel_bottom_label))
        self.dofs_cladding_bottom = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.cladding_bottom_label))
        self.dofs_fuel_outer      = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.fuel_outer_label))
        self.dofs_cladding_outer  = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.cladding_outer_label))
        self.dofs_cladding_inner  = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.cladding_inner_label))
        self.dofs_cladding_top    = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.cladding_top_label))
        self.dofs_fuel_top        = dolfinx.fem.locate_dofs_topological(V=self.V_m, entity_dim=self.fdim, entities=self.ft.find(self.fuel_top_label))

        # -----------------------------
        # Dirichlet boundary conditions
        # -----------------------------

        # Mechanical boundary conditions
        self.bc_mech_fuel_bottom = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_fuel_bottom, self.V_m)
        self.bc_mech_cladding_bottom = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_cladding_bottom, self.V_m)

        self.bc_mech = [self.bc_mech_fuel_bottom, self.bc_mech_cladding_bottom]

        # Thermal boundary conditions
        self.T_bc = temperature
        self.bc_Tco = dolfinx.fem.dirichletbc(temperature, self.dofs_cladding_outer, self.V_t)
        self.bc_Tfo = dolfinx.fem.dirichletbc(temperature, self.dofs_fuel_outer, self.V_t)
    
        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)
        self.gap_traction = dolfinx.fem.Constant(self.mesh, gap_pressure) * self.normal
        self.coolant_traction = dolfinx.fem.Constant(self.mesh, coolant_pressure) * self.normal

        self.heat_flux = self.lhr / (2.0 * np.pi * self.cladding_inner_radius) # W/m2

    def assemble_and_solve(self):
        """
        Assemble and solve the system for the linear thermo-elastic problem.
        """
        # Integration mesures "ds", "dx"
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure("dx", domain=self.mesh, metadata=metadata)
        self.ds = ufl.Measure("ds", domain=self.mesh, metadata=metadata, subdomain_data=self.ft)
        self.ds_ci = self.ds(self.cladding_inner_label)
        self.ds_co = self.ds(self.cladding_outer_label)
        self.ds_fo = self.ds(self.fuel_outer_label)

        self.solver_options()
       
        self.thermal_solver()
        self.mechanical_solver()
        
    def solver_options(self):
        # Direct solver + MUMPS
        self.petsc_options_lu = {
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-7,
            "ksp_max_it": 2000,
            "ksp_monitor": None
        }

        # CG-solver + jacobi
        self.petsc_options_cg = {
            "ksp_type": "cg",
            "pc_type": "jacobi",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-7,
            "ksp_max_it": 2000,
            "ksp_monitor": None
        }

        # GMRES solver + pre-conditioner AMG (Algebraic Multigrid)
        self.petsc_options_gmres = {
            "ksp_type": "gmres",
            "pc_type": "hypre",
            "pc_hypre_type": "boomeramg",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-7,
            "ksp_max_it": 2000,
            "pc_type": None,
        }

        # GMRES iterative solver + ILU pre-conditioner
        self.petsc_options_ilu = {
            "ksp_type": "gmres",
            "pc_type": "ilu",
            "ksp_rtol": 1e-7,
            "ksp_atol": 1e-7,
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

        # Characteristic time, fuel thermal evolution
        characteristic_time = self.fuel["rho"] * self.fuel["cp"] / self.fuel["k"] * self.area_fuel
        print(f"Characteristic time (fuel): {characteristic_time:.3e} s")

        # Characteristic time, cladding thermal evolution
        characteristic_time = self.cladding["rho"] * self.cladding["cp"] / self.cladding["k"] * self.area_cladding
        print(f"Characteristic time (cladding): {characteristic_time:.3e} s")

        a_thermal_fuel = ufl.dot(self.fuel["k"] * ufl.grad(self.u_t), ufl.grad(self.v_t)) * self.dx_fuel
        a_thermal_cladding = ufl.dot(self.cladding["k"] * ufl.grad(self.u_t), ufl.grad(self.v_t)) * self.dx_cladding
        self.a_thermal = a_thermal_fuel + a_thermal_cladding

        L_thermal_fuel = dolfinx.fem.Constant(self.mesh, self.q_third) * self.v_t * self.dx_fuel
        L_thermal_cladding = self.heat_flux * self.v_t * self.ds_ci
        self.L_thermal = L_thermal_fuel + L_thermal_cladding

        integral_value = dolfinx.fem.assemble_scalar(dolfinx.fem.form(self.heat_flux * self.ds_ci))
        print(f"Total power from the integrated value of the heat flux on the inner wall: {integral_value:.3f} W")
        print(f"Total power from the analytical value: {self.lhr * self.height:.3f} W")

        self.thermal_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_thermal,
            self.L_thermal,
            bcs=[self.bc_Tco, self.bc_Tfo],
            petsc_options=self.petsc_options_lu,
            jit_options=self.jit_options,
        )

        self.temperature = self.thermal_problem.solve()
        self.temperature.x.scatter_forward()

        temperature_cladding_inner = self.temperature.x.array[self.dofs_cladding_inner]
        print(f"Average temperature on inner cladding: {np.mean(temperature_cladding_inner):.2f} K")

        # Gap conductance
        temperature_gap = np.mean(temperature_cladding_inner)
        for i in range(15):
            k_gas = self.filling_gas['A'] * 1e-4 * temperature_gap**0.79 # W/m-K
            h_gap = k_gas / self.gap_thickness

            thermal_resistance_gap = 1./(2.0 * np.pi * self.fuel_outer_radius * h_gap )

            d_temperature_gap = self.lhr * thermal_resistance_gap

            print(f"Analytical thermal jump - gap: {self.lhr * thermal_resistance_gap:.2f} K")

            temperature_gap = np.mean(temperature_cladding_inner) + d_temperature_gap

        self.bc_Tfo = dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, np.mean(temperature_cladding_inner) + d_temperature_gap), self.dofs_fuel_outer, self.V_t)

        self.thermal_problem = dolfinx.fem.petsc.LinearProblem(
            self.a_thermal,
            self.L_thermal,
            bcs=[self.bc_Tco, self.bc_Tfo],
            petsc_options=self.petsc_options_lu,
            jit_options=self.jit_options,
        )

        self.temperature = self.thermal_problem.solve()
        self.temperature.x.scatter_forward()
        
        self.temperature.name = "temperature"

        self.assemble_and_print_matrix(a = self.a_thermal, L = self.L_thermal, filename = 'thermal_stiffness_matrix.txt')

        # print(f"Analytical thermal jump: {self.lhr / ( 4.0 * np.pi * self.k) * (1. - self.inner_radius**2/(self.outer_radius**2 - self.inner_radius**2) * 2. * np.log(self.outer_radius/self.inner_radius)):.2f} K")
        thermal_resistance_cladding = np.log(self.cladding_outer_radius/self.cladding_inner_radius) / (2.0*np.pi*self.cladding["k"])
        thermal_resistance_fuel = 1./(4.*np.pi*self.fuel["k"])
        print(f"Analytical thermal jump - cladding: {self.lhr * thermal_resistance_cladding:.2f} K")
        print(f"Analytical thermal jump - gap: {self.lhr * thermal_resistance_gap:.2f} K")
        print(f"Analytical thermal jump - fuel: {self.lhr * thermal_resistance_fuel:.2f} K")
        print(f"Min temperature from thermal solver: {np.min(self.temperature.x.array):.2f} K")
        print(f"Max temperature from thermal solver: {np.max(self.temperature.x.array):.2f} K")

        print("\n")

    def mechanical_solver(self):

        print("-----------------")
        print("Mechanical solver")
        print("-----------------")

        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
        dx_fuel = dx_subdomain(1)
        dx_cladding = dx_subdomain(2)

        def epsilon(u):
            """
            Strain-displacement compatibility relation.
            """
            return ufl.sym(ufl.grad(u))

        def sigma(u, lmbda, mu):
            """
            Stress-strain constitutive relation.
            """
            return lmbda * ufl.tr(epsilon(u)) * ufl.Identity(len(u)) + 2.0 * mu * epsilon(u)

        # Mechanical stiffness (bilinear form) split into fuel and cladding contributions.
        a_mech_fuel = ufl.inner(sigma(self.u_m, self.fuel["lmbda"], self.fuel["mu"]), epsilon(self.v_m)) * dx_fuel
        a_mech_cladding = ufl.inner(sigma(self.u_m, self.cladding["lmbda"], self.cladding["mu"]), epsilon(self.v_m)) * dx_cladding
        self.a_mech = a_mech_fuel + a_mech_cladding

        # Body force contributions: different densities for fuel and cladding.
        L_body_forces_fuel = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -self.fuel["rho"] * self.g)), self.v_m) * dx_fuel
        L_body_forces_cladding = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -self.cladding["rho"] * self.g)), self.v_m) * dx_cladding
        self.L_body_forces = L_body_forces_fuel + L_body_forces_cladding

        # Traction forces remain defined on the boundary (using ds measures already set up).
        self.L_tractions = (ufl.dot(self.gap_traction, self.v_m) * self.ds_fo +
                            ufl.dot(self.gap_traction, self.v_m) * self.ds_ci +
                            ufl.dot(self.coolant_traction, self.v_m) * self.ds_co)

        # Thermal expansion contributions: split the term for fuel and cladding.
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        L_thermal_exp_fuel = (3 * self.fuel["lmbda"] + 2 * self.fuel["mu"]) * self.fuel["alpha"] * (self.temperature - T_ref) * ufl.div(self.v_m) * dx_fuel
        L_thermal_exp_cladding = (3 * self.cladding["lmbda"] + 2 * self.cladding["mu"]) * self.cladding["alpha"] * (self.temperature - T_ref) * ufl.div(self.v_m) * dx_cladding
        self.L_thermal_exp = L_thermal_exp_fuel + L_thermal_exp_cladding

        # Total mechanical load is the sum of the body forces, thermal expansion, and traction forces.
        self.L_mech = self.L_body_forces + self.L_thermal_exp + self.L_tractions

        self.assemble_and_print_matrix(a = self.a_mech, L = self.L_mech, filename = 'mech_stiffness_matrix.txt')

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
        self.displacement.x.scatter_forward()
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)

        self.stress_fuel = sigma(self.displacement, self.fuel["lmbda"], self.fuel["mu"]) - (3 * self.fuel["lmbda"] + 2 * self.fuel["mu"]) * self.fuel["alpha"] * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))

        self.stress_cladding = sigma(self.displacement, self.cladding["lmbda"], self.cladding["mu"]) - (3 * self.cladding["lmbda"] + 2 * self.cladding["mu"]) * self.cladding["alpha"] * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))

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

        return self.displacement, self.stress_fuel, self.stress_cladding, self.strain, self.temperature

    def vonMises_equivalent(self):
        """
        Compute and visualize the Von Mises stress field.
        """

        # Compute deviatoric stress
        deviatoric_stress = self.stress_cladding - (1.0 / 3.0) * ufl.tr(self.stress_cladding) * ufl.Identity(len(self.displacement))
        
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
