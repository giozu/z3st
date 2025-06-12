import dolfinx.fem.function
import dolfinx.fem.petsc
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from core.finite_element_setup import FiniteElementSetup

class FullCylinder_time(FiniteElementSetup):
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

    def temporal_solver(self, ti, tf, n_steps):
        """
        Time-solver settings
        """

        self.ti = ti
        self.tf = tf
        self.n_steps = n_steps
        self.dt = (tf-ti) / n_steps

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

        # Mesure "ds"
        self.ds = ufl.Measure("ds", domain=self.mesh)
        ds_lateral = self.ds(self.lateral_label)
        ds = self.ds

        print("---------------")
        print("Thermal problem")
        print("---------------")

        characteristic_time = self.rho * self.cp / self.k * self.area
        print(f"Characteristic time: {characteristic_time:.3e} s")

        # Define solution variable, and interpolate initial solution for visualization in Paraview
        xdmf_temperature = dolfinx.io.XDMFFile(self.mesh.comm, "temperature.xdmf", "w")
        xdmf_temperature.write_mesh(self.mesh)

        xdmf_displacement = dolfinx.io.XDMFFile(self.mesh.comm, "displacement.xdmf", "w")
        xdmf_displacement.write_mesh(self.mesh)

        t = self.ti

        self.temperature = dolfinx.fem.Function(self.V_t)
        self.temperature.name = "temperature"
        self.temperature.interpolate(lambda x: self.T_bc*np.ones(x.shape[1]))
        xdmf_temperature.write_function(self.temperature, t)

        T_0 = dolfinx.fem.Function(self.V_t)
        T_0.name = "temperature_0"
        T_0.interpolate(lambda x: self.T_bc*np.ones(x.shape[1]))

        self.a_thermal = (self.u_t * self.v_t * ufl.dx) + self.dt * ufl.dot((self.k / (self.rho * self.cp)) * ufl.grad(self.u_t), ufl.grad(self.v_t)) * ufl.dx
        self.L_thermal = (T_0 + self.dt / (self.rho * self.cp) * dolfinx.fem.Constant(self.mesh, self.q_third)) * self.v_t * ufl.dx
        
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
        self.displacement = dolfinx.fem.Function(self.V_m)

        self.a_mech = ufl.inner(sigma(self.displacement), epsilon(self.v_m)) * ufl.dx
        self.L_body_forces = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, - self.rho * self.g)), self.v_m) * ufl.dx
        self.L_tractions = ufl.dot(self.traction, self.v_m) * ds
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        self.L_thermal_exp = (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.div(self.v_m) * ufl.dx
        self.L_mech = self.L_body_forces + self.L_thermal_exp + self.L_tractions

        F = self.a_mech - self.L_mech
        problem_mech = dolfinx.fem.petsc.NonlinearProblem(F, self.displacement, self.bc_mech)

        # Newton solver
        from dolfinx.nls.petsc import NewtonSolver
        solver_mech = NewtonSolver(self.mesh.comm, problem_mech)
        solver_mech.atol = 1e-8
        solver_mech.rtol = 1e-8
        solver_mech.convergence_criterion = "incremental"

        # linear algebra structures
        bilinear_form_thermal = dolfinx.fem.form(self.a_thermal)
        linear_form_thermal = dolfinx.fem.form(self.L_thermal)

        A_th = dolfinx.fem.petsc.assemble_matrix(bilinear_form_thermal, bcs = self.bc_thermal)
        A_th.assemble()
        b_th = dolfinx.fem.petsc.create_vector(linear_form_thermal)

        solver_th = PETSc.KSP().create(self.mesh.comm)
        solver_th.setOperators(A_th)
        solver_th.setType(PETSc.KSP.Type.PREONLY)
        solver_th.getPC().setType(PETSc.PC.Type.LU)

        print(f"Time: {t} s")
        for i in range(self.n_steps):
            t += self.dt
            print(f"Time: {t} s")

            # Update the right hand side reusing the initial vector
            with b_th.localForm() as loc_b:
                loc_b.set(0)

            dolfinx.fem.petsc.assemble_vector(b_th, linear_form_thermal)

            # Apply Dirichlet boundary condition to the vector
            dolfinx.fem.petsc.apply_lifting(b_th, [bilinear_form_thermal], [self.bc_thermal])
            b_th.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            dolfinx.fem.petsc.set_bc(b_th, self.bc_thermal)

            # Solve linear problem
            solver_th.solve(b_th, self.temperature.x.petsc_vec)
            self.temperature.x.scatter_forward()

            # Update solution at previous time step (u_n)
            T_0.x.array[:] = self.temperature.x.array

            num_its, converged = solver_mech.solve(self.displacement)
            assert (converged)
            self.displacement.x.scatter_forward()

            # Write solution to file
            xdmf_temperature.write_function(self.temperature, t)
            xdmf_displacement.write_function(self.displacement, t)

        xdmf_temperature.close()
        xdmf_displacement.close()

        thermal_resistance = 4.0 * np.pi * self.k
        print("\n")
        print(f"Analytical thermal jump: {self.lhr / thermal_resistance:.2f} K")
        print("\n")

        self.stress = self.lmbda * ufl.tr(epsilon(self.displacement)) * ufl.Identity(len(self.displacement)) + 2.0 * self.mu * epsilon(self.displacement) - (3 * self.lmbda + 2 * self.mu) * self.alpha * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))
        self.strain = epsilon(self.displacement)

        # Solver options
        # solver = self.thermal_problem.solver
        # viewer = PETSc.Viewer().createASCII("log.thermal_solver")
        # solver.view(viewer)

        # solver = self.mech_problem.solver
        # viewer = PETSc.Viewer().createASCII("log.mech_solver")
        # solver.view(viewer)

    def get_results(self):
        """
        Get the computed displacement, stress, and strain fields.

        Returns:
        displacement (Function): The displacement field.
        stress (ufl.Expr): The stress field.
        strain (ufl.Expr): The strain field.
        """

        return self.displacement, self.stress, self.strain, self.temperature
