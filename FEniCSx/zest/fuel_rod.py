import dolfinx.fem.petsc
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from finiteElementSettings import FiniteElementSetting
from geometryManager import GeometryManager

from utils import get_ksp_reason_name

class FuelRod(FiniteElementSetting):
    def volume(self, fuel_outer_radius, cladding_inner_radius, cladding_outer_radius, height, log_path="output/simulation_log.txt"):
        """
        Compute the volume of the domain.
        """

        self.fuel_outer_radius = fuel_outer_radius
        self.cladding_inner_radius = cladding_inner_radius
        self.cladding_outer_radius = cladding_outer_radius
        self.height = height
        self.gap_thickness = self.cladding_inner_radius - self.fuel_outer_radius
        self.area_fuel = np.pi * self.fuel_outer_radius**2
        self.area_cladding = np.pi * (self.cladding_outer_radius**2 - self.cladding_inner_radius**2)

        # Integration measures the domain
        dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
        self.dx_fuel = dx_subdomain(1)      # Marker ID for fuel region
        self.dx_cladding = dx_subdomain(2)  # Marker ID for cladding region

        # Compute mesh volumes
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))
        volume_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_fuel))
        volume_cladding = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_cladding))
        
        analytical_volume = (self.area_fuel + self.area_cladding)*self.height
        rel_error = abs(volume - analytical_volume) / analytical_volume * 100

        with open(log_path, "a") as f:
            f.write("📦 Mesh Volumes\n")
            f.write("──────────────────────────────────────────────\n")
            f.write(f"{'Fuel volume (from mesh)':<32}: {volume_fuel*1e9:.3e} mm³\n")
            f.write(f"{'Cladding volume (from mesh)':<32}: {volume_cladding*1e9:.3e} mm³\n")
            f.write(f"{'Fuel volume (analytical)':<32}: {self.area_fuel * self.height * 1e9:.3e} mm³\n")
            f.write(f"{'Cladding volume (analytical)':<32}: {self.area_cladding * self.height * 1e9:.3e} mm³\n")
            f.write(f"{'Relative error':<32}: {rel_error:.3e} %\n\n")

    def initialize_geometry(self):
        self.geometry_manager = GeometryManager(self.mesh, self.fuel_outer_radius, self.cladding_inner_radius, self.cladding_outer_radius, self.height)
        self.geometry_manager.initialize_geometry()

    def update_state(self):
        """
        Update derived physical quantities after mesh deformation.
        Must be called after update_geometry().
        """

        self.fuel_outer_radius = self.geometry_manager.fuel_outer_radius
        self.cladding_inner_radius = self.geometry_manager.cladding_inner_radius
        self.cladding_outer_radius = self.geometry_manager.cladding_outer_radius
        self.height = self.geometry_manager.height

        self.area_fuel = np.pi * self.fuel_outer_radius**2
        self.q_third = self.lhr / self.area_fuel
        self.heat_flux = self.lhr / (2.0 * np.pi * self.cladding_inner_radius)
        self.gap_thickness = self.cladding_inner_radius - self.fuel_outer_radius

        print(f"   fuel outer radius = {self.fuel_outer_radius:.6e} m")
        print(f"   cladding inner radius = {self.cladding_inner_radius:.6e} m")
        print(f"   cladding outer radius = {self.cladding_outer_radius:.6e} m")
        print(f"   fuel area = {self.area_fuel:.4e} m²")
        print(f"   gap thickness = {self.gap_thickness*1e6:.3f} micron")
        print(f"   q''' = {self.q_third:.2f} W/m³")
        print(f"   heat flux at cladding inner surface = {self.heat_flux:.2f} W/m²")

    def parameters(self, fuel, cladding, filling_gas, g, lhr, T_ref, log_path="output/simulation_log.txt"):
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
        self.filling_gas = filling_gas

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
        # self.q_third = self.lhr / self.area_fuel
        # self.heat_flux = self.lhr / (2.0 * np.pi * self.cladding_inner_radius) # W/m2

        # characteristic_time = self.fuel["rho"] * self.fuel["cp"] / self.fuel["k"] * self.area_fuel
        # print(f"\nCharacteristic time (fuel): {characteristic_time:.3e} s")
        # characteristic_time = self.cladding["rho"] * self.cladding["cp"] / self.cladding["k"] * self.area_cladding
        # print(f"Characteristic time (cladding): {characteristic_time:.3e} s")

        with open(log_path, "a") as f:
            f.write("🔥 Fuel Parameters\n")
            f.write("──────────────────────────────────────────────\n")
            f.write(f"{'Inner radius':<32}: {self.fuel_outer_radius:.3e} m\n")
            f.write(f"{'Density':<32}: {self.fuel['rho']:.2f} kg/m³\n")
            f.write(f"{'Elastic modulus':<32}: {self.fuel['E']/1e9:.2f} GPa\n")
            f.write(f"{'Poisson\'s ratio':<32}: {self.fuel['nu']:.2f}\n")
            f.write(f"{'Bulk modulus':<32}: {self.fuel['bulk_modulus']/1e9:.2f} GPa\n")
            f.write(f"{'Lamé λ, μ':<32}: {self.fuel['lmbda']/1e9:.2f} GPa, {self.fuel['mu']/1e9:.2f} GPa\n")
            f.write(f"{'Thermal conductivity':<32}: {self.fuel['k']:.2f} W/mK\n")
            f.write(f"{'Thermal expansion coefficient':<32}: {self.fuel['alpha']:.2e} 1/K\n")
            f.write(f"{'Specific heat capacity':<32}: {self.fuel['cp']:.2f} J/kgK\n\n")

            f.write("🧊 Cladding Parameters\n")
            f.write("──────────────────────────────────────────────\n")
            f.write(f"{'Inner radius':<32}: {self.cladding_inner_radius:.3e} m\n")
            f.write(f"{'Outer radius':<32}: {self.cladding_outer_radius:.3e} m\n")
            f.write(f"{'Thickness':<32}: {(self.cladding_outer_radius - self.cladding_inner_radius)*1e6:.3f} micron\n")
            f.write(f"{'Density':<32}: {self.cladding['rho']:.2f} kg/m³\n")
            f.write(f"{'Elastic modulus':<32}: {self.cladding['E']/1e9:.2f} GPa\n")
            f.write(f"{'Poisson\'s ratio':<32}: {self.cladding['nu']:.2f}\n")
            f.write(f"{'Bulk modulus':<32}: {self.cladding['bulk_modulus']/1e9:.2f} GPa\n")
            f.write(f"{'Lamé λ, μ':<32}: {self.cladding['lmbda']/1e9:.2f} GPa, {self.cladding['mu']/1e9:.2f} GPa\n")
            f.write(f"{'Thermal conductivity':<32}: {self.cladding['k']:.2f} W/mK\n")
            f.write(f"{'Thermal expansion coefficient':<32}: {self.cladding['alpha']:.2e} 1/K\n")
            f.write(f"{'Specific heat capacity':<32}: {self.cladding['cp']:.2f} J/kgK\n\n")

            f.write("🧪 Gas and Global Conditions\n")
            f.write("──────────────────────────────────────────────\n")
            f.write(f"{'Gap thickness':<32}: {self.gap_thickness*1e6:.3f} micron\n")
            f.write(f"{'Filling gas':<32}: {self.filling_gas['name']}\n")
            f.write(f"{'Gravitational acceleration':<32}: {self.g:.2f} m/s²\n")
            f.write(f"{'Reference temperature':<32}: {self.T_ref:.2f} K\n\n")

    def boundary_conditions(self, gap_pressure, coolant_pressure, temperature, log_path="output/simulation_log.txt"):
        """
        Define and set the boundary conditions.
        """

        with open(log_path, "a") as f:
            f.write("🔧 Boundary Conditions\n")
            f.write("──────────────────────────────────────────────\n")
            f.write(f"{'Gap pressure':<32}: {gap_pressure*1e-6:.2f} MPa\n")
            f.write(f"{'Coolant pressure':<32}: {coolant_pressure*1e-6:.2f} MPa\n")
            f.write(f"{'Boundary temperature':<32}: {temperature:.2f} K\n\n")

        self.fuel_bottom_label = 1
        self.cladding_bottom_label = 2
        self.fuel_outer_label = 3
        self.cladding_outer_label = 4
        self.cladding_inner_label = 5
        self.cladding_top_label = 6
        self.fuel_top_label = 7

        self.dx = ufl.Measure("dx", domain=self.mesh, metadata={"quadrature_degree": 4})
        self.ds = ufl.Measure("ds", domain=self.mesh, metadata={"quadrature_degree": 4}, subdomain_data=self.ft)
        self.ds_ci = self.ds(self.cladding_inner_label)
        self.ds_co = self.ds(self.cladding_outer_label)
        self.ds_fo = self.ds(self.fuel_outer_label)

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
        self.T_bc_co = dolfinx.fem.Constant(self.mesh, temperature)
        self.bc_Tco = dolfinx.fem.dirichletbc(self.T_bc_co, self.dofs_cladding_outer, self.V_t)

        self.T_bc_fo = dolfinx.fem.Constant(self.mesh, temperature)
        self.bc_Tfo = dolfinx.fem.dirichletbc(self.T_bc_fo, self.dofs_fuel_outer, self.V_t)

        # ---------------------------
        # Neumann boundary conditions
        # ---------------------------

        self.normal = ufl.FacetNormal(self.mesh)
        self.gap_traction = dolfinx.fem.Constant(self.mesh, gap_pressure) * self.normal
        self.coolant_traction = dolfinx.fem.Constant(self.mesh, coolant_pressure) * self.normal

    def initialize_fields(self):
        """
        Initialize fields like temperature and displacement to physically meaningful values.
        This should be called once before the time loop.
        """
        print("🧪 Initializing temperature and displacement fields...")

        # temperature
        self.temperature = dolfinx.fem.Function(self.V_t)
        self.temperature.name = "temperature"
        self.temperature.x.array[:] = float(self.T_bc_fo.value)
        print(f"\t🌡️  Temperature initialized at {float(self.T_bc_fo.value):.2f} K")

        # displacement
        self.displacement = dolfinx.fem.Function(self.V_m)
        self.displacement.name = "displacement"
        self.displacement.x.array[:] = 0.0
        print(f"\t📦 Displacement initialized at 0.0 m")

        self.total_thermal_iterations = 0.0
        self.total_mech_iterations = 0.0


    def solve(self, max_iter=10, tol=1e-2):
        """
        Iterative thermo-mechanical solver (steady-state) with fixed mesh.
        Converges when the relative change in maximum radial displacement is below tolerance.
        """

        print("======================")
        print("Thermo-mechanical loop")
        print("======================")

        converged = False
        u_r_max_prev = None

        for it in range(1, max_iter + 1):
            print(f"\n🔁 Thermo-mechanical iteration {it}")

            it_therm = self.thermal_solver()
    
            # T = dolfinx.fem.assemble_scalar(dolfinx.fem.form(self.temperature * self.dx_fuel))
            # volume_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_fuel))
            # T_avg = T / volume_fuel
            # print(f"📊 T_avg fuel (after thermal solver) = {T_avg:.2f} K")

            # relaxation factor
            omega = 0.7
            if it > 1:
                self.temperature.x.array[:] = (1 - omega) * temperature_prev + omega * self.temperature.x.array
            temperature_prev = self.temperature.x.array.copy()

            it_mech = self.mechanical_solver()
            
            # Compute radial displacement
            u = self.displacement.x.array.reshape(self.mesh.geometry.x.shape)
            disp_r = np.sqrt(u[:, 0]**2 + u[:, 1]**2)
            u_r_max = np.max(disp_r)

            # Convergence check
            if u_r_max_prev is not None:

                den = max(abs(u_r_max_prev), 1e-12)
                delta_u = abs(u_r_max - u_r_max_prev) / den

                print(f"\tΔu_r_rel = {delta_u:.2e}")

                if delta_u < tol:
                    print(f" >> Thermo-mechanical convergence achieved.")
                    print(f" >> Iterations: thermal = {it_therm}, mechanical = {it_mech}")

                    self.total_thermal_iterations += it_therm
                    self.total_mech_iterations += it_mech

                    converged = True
                    break
            else:
                print(f"\tInitial u_r_max = {u_r_max:.3e} m")

            # Logging
            # r_geo = np.linalg.norm(self.mesh.geometry.x[:, :2], axis=1)
            # T_center = self.temperature.x.array[np.argmin(r_geo)]
            # print(f"   → Max u_r = {u_r_max:.3e} m")
            # print(f"   → T_center = {T_center:.1f} K")

            u_r_max_prev = u_r_max
        
            gap = self.cladding_inner_radius - self.fuel_outer_radius - 2 * u_r_max
            if gap < 0:
                print(f"⚠️ Warning: apparent gap closure (gap = {gap*1e6:.2f} μm)")

        if not converged:
            print("⚠️ Warning: convergence not reached in max_iter.")

    # def solver_options(self):
        # ksp_rtol = 1e-6
        # ksp_atol = 1e-6

        # # Direct solver + MUMPS
        # self.petsc_options_lu = {
        #     "ksp_type": "preonly",
        #     "pc_type": "lu",
        #     "pc_factor_mat_solver_type": "mumps",
        #     "ksp_rtol": ksp_rtol,
        #     "ksp_atol": ksp_atol,
        #     "ksp_max_it": 2000,
        #     "ksp_monitor": None
        # }

        # # CG-solver + jacobi
        # self.petsc_options_cg = {
        #     "ksp_type": "cg",
        #     "pc_type": "jacobi",
        #     "ksp_rtol": ksp_rtol,
        #     "ksp_atol": ksp_atol,
        #     "ksp_max_it": 2000,
        #     "ksp_monitor": None
        # }

        # # GMRES solver + pre-conditioner AMG (Algebraic Multigrid)
        # self.petsc_options_gmres = {
        #     "ksp_type": "gmres",
        #     "pc_type": "hypre",
        #     "pc_hypre_type": "boomeramg",
        #     "ksp_rtol": ksp_rtol,
        #     "ksp_atol": ksp_atol,
        #     "ksp_max_it": 2000,
        #     "pc_type": None,
        # }

        # # GMRES iterative solver + ILU pre-conditioner
        # self.petsc_options_ilu = {
        #     "ksp_type": "gmres",
        #     "pc_type": "ilu",
        #     "ksp_rtol": ksp_rtol,
        #     "ksp_atol": ksp_atol,
        #     "ksp_max_it": 2000,
        #     "ksp_monitor": None
        # }

        # self.jit_options = {
        #     "timeout": 20,
        # }

    def assemble_thermal_system(self):
        """
        Assemble bilinear and linear forms for the thermal problem.
        Includes contributions from both fuel and cladding regions.
        """
        T = self.temperature
        # k_fuel = self.fuel["k"]
        k_fuel = 1.0 / (0.0452 + 2.46e-4 * T)
        k_clad = self.cladding["k"]

        k_fuel_expr = dolfinx.fem.Expression(k_fuel, self.V_scalar.element.interpolation_points())
        k_func = dolfinx.fem.Function(self.V_scalar)
        k_func.interpolate(k_fuel_expr)

        k_avg_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(k_func * self.dx_fuel)) / dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_fuel))

        print(f"\tAverage fuel thermal conductivity = {k_avg_fuel:.3f} W/mK")

        # Bilinear form (stiffness matrix)
        a_fuel = ufl.dot(k_fuel * ufl.grad(self.u_t), ufl.grad(self.v_t)) * self.dx_fuel
        a_clad = ufl.dot(k_clad * ufl.grad(self.u_t), ufl.grad(self.v_t)) * self.dx_cladding
        self.a_thermal = a_fuel + a_clad

        # Linear form (source terms)
        L_fuel = dolfinx.fem.Constant(self.mesh, self.q_third) * self.v_t * self.dx_fuel
        L_clad = self.heat_flux * self.v_t * self.ds_ci
        self.L_thermal = L_fuel + L_clad

    def solve_temperature_field(self, bcs):
        """
        Solve the thermal problem using PETSc KSP solver.
        Takes updated boundary conditions as input.
        """
        # Variational problem
        F = self.a_thermal - self.L_thermal
        a_form = dolfinx.fem.form(ufl.lhs(F))
        L_form = dolfinx.fem.form(ufl.rhs(F))

        # Assemble matrix and RHS
        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)
        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=bcs)
        A.assemble()
        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [bcs])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, bcs)

        # Solve system
        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.GMRES)
        solver.getPC().setType(PETSc.PC.Type.HYPRE)
        solver.setTolerances(rtol=1e-5, atol=1e-5, max_it=500)
        # solver.setMonitor(lambda ksp, its, rnorm: print(f"[KSP] Iteration {its}, residual = {rnorm:.3e}"))

        solver.solve(b, self.temperature.x.petsc_vec)
        its = solver.getIterationNumber()

        reason = solver.getConvergedReason()
        its = solver.getIterationNumber()

        print(f"\t[KSP] Converged in {its} iterations, reason = {reason}: {get_ksp_reason_name(reason)}")

        self.temperature.x.scatter_forward()

        return its

    def compute_gap_conductance(self, T_gap_initial):
        """
        Iteratively computes the temperature drop across the gap (ΔT_gap),
        given an initial guess for the gap temperature.
        """
        if self.lhr < 1e-6:
            return 0.0

        T_gap = T_gap_initial
        delta_T_old = 1.0
        tol = 1e-3
        error = 1.0
        max_iter = 20
        iteration = 0

        while error > tol and iteration < max_iter:
            iteration += 1
            k_gas = self.filling_gas['A'] * T_gap**0.79
            h_gap = k_gas / self.gap_thickness
            R_gap = 1.0 / (2.0 * np.pi * self.fuel_outer_radius * h_gap)
            delta_T_new = self.lhr * R_gap
            error = abs(delta_T_new - delta_T_old) / max(abs(delta_T_old), 1e-6)

            # print(f"  Iteration {iteration:02d} | ΔT_gap: {delta_T_new:.2f} K | Error: {error:.2e} | h_gap: {h_gap:.2f} W/m²K")
            # print(f"     → r_fuel_outer = {self.fuel_outer_radius:.6e} m")
            # print(f"     → LHR = {self.lhr:.2f} W/m")
            # print(f"     → gap_thickness = {self.gap_thickness:.6e} m")

            delta_T_old = delta_T_new
            T_gap = T_gap_initial + delta_T_new  # update based on current ΔT

        print(f"\tIteration {iteration:02d} | ΔT_gap: {delta_T_new:.2f} K | Error: {error:.2e} | h_gap: {h_gap:.2f} W/m²K")
        return delta_T_new

    def thermal_solver(self):
        """
        Thermal solver:
        1. Solves temperature with initial BCs.
        2. Computes gap conductance using the cladding inner surface temperature.
        3. Updates fuel outer surface BC and resolves temperature field.
        """

        print("Thermal solver")
        print("--------------")

        if self.gap_thickness < 0:
            raise RuntimeError(f"❌ Gap thickness < 0: {self.gap_thickness*1e6:.2f} microns")

        # Assemble forms for a and L (dependent on q_third and heat_flux)
        self.assemble_thermal_system()

        # Initial solve with current boundary conditions
        bcs = [self.bc_Tco, self.bc_Tfo]
        its = self.solve_temperature_field(bcs)

        # Extract temperature at cladding inner surface
        T_clad_inner = np.mean(self.temperature.x.array[self.dofs_cladding_inner])
        T_fuel_outer = float(self.T_bc_fo.value)

        T_gap_guess = 0.5 * (T_clad_inner + T_fuel_outer)
        delta_T_gap = self.compute_gap_conductance(T_gap_guess)

        # Update fuel outer surface temperature boundary condition
        self.T_bc_fo = dolfinx.fem.Constant(self.mesh, T_clad_inner + delta_T_gap)
        self.bc_Tfo = dolfinx.fem.dirichletbc(self.T_bc_fo, self.dofs_fuel_outer, self.V_t)

        # Solve again with updated BC
        bcs = [self.bc_Tco, self.bc_Tfo]
        its += self.solve_temperature_field(bcs)

        # print(f"Cladding outer temperature = {float(self.T_bc_co.value)} K")
        # print(f"Cladding inner temperature = {np.mean(self.temperature.x.array[self.dofs_cladding_inner])} K")
        # print(f"Fuel outer temperature = {float(self.T_bc_fo.value)} K")

        # Postprocess (optional): print analytical reference values
        R_clad = np.log(self.cladding_outer_radius / self.cladding_inner_radius) / (2.0 * np.pi * self.cladding["k"])
        R_fuel = 1.0 / (4.0 * np.pi * self.fuel["k"])

        # print(f"Analytical temperature drop - Cladding: {self.lhr * R_clad:.2f} K")
        # print(f"Analytical temperature drop - Gap:      {delta_T_gap:.2f} K")
        # print(f"Analytical temperature drop - Fuel:     {self.lhr * R_fuel:.2f} K")

        return its

    def mechanical_solver(self):

        print("Mechanical solver")
        print("-----------------")

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
        a_mech_fuel = ufl.inner(sigma(self.u_m, self.fuel["lmbda"], self.fuel["mu"]), epsilon(self.v_m)) * self.dx_fuel
        a_mech_cladding = ufl.inner(sigma(self.u_m, self.cladding["lmbda"], self.cladding["mu"]), epsilon(self.v_m)) * self.dx_cladding
        self.a_mech = a_mech_fuel + a_mech_cladding

        # Body force contributions: different densities for fuel and cladding.
        L_body_forces_fuel = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -self.fuel["rho"] * self.g)), self.v_m) * self.dx_fuel
        L_body_forces_cladding = ufl.dot(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, -self.cladding["rho"] * self.g)), self.v_m) * self.dx_cladding
        self.L_body_forces = L_body_forces_fuel + L_body_forces_cladding

        # Traction forces remain defined on the boundary (using ds measures already set up).
        self.L_tractions = (ufl.dot(self.gap_traction, self.v_m) * self.ds_fo +
                            ufl.dot(self.gap_traction, self.v_m) * self.ds_ci +
                            ufl.dot(self.coolant_traction, self.v_m) * self.ds_co)

        # Thermal expansion contributions: split the term for fuel and cladding.
        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        L_thermal_exp_fuel = (3 * self.fuel["lmbda"] + 2 * self.fuel["mu"]) * self.fuel["alpha"] * (self.temperature - T_ref) * ufl.div(self.v_m) * self.dx_fuel
        L_thermal_exp_cladding = (3 * self.cladding["lmbda"] + 2 * self.cladding["mu"]) * self.cladding["alpha"] * (self.temperature - T_ref) * ufl.div(self.v_m) * self.dx_cladding
        self.L_thermal_exp = L_thermal_exp_fuel + L_thermal_exp_cladding

        # Total mechanical load is the sum of the body forces, thermal expansion, and traction forces.
        self.L_mech = self.L_body_forces + self.L_thermal_exp + self.L_tractions

        # self.assemble_and_print_matrix(a = self.a_mech, L = self.L_mech, filename = 'mech_stiffness_matrix.txt')

        # Mechanical solver - low-level PETSc setup
        a2 = dolfinx.fem.form(self.a_mech)
        L2 = dolfinx.fem.form(self.L_mech)

        A2 = dolfinx.fem.petsc.create_matrix(a2)
        b2 = dolfinx.fem.petsc.create_vector(L2)

        solver2 = PETSc.KSP().create(self.mesh.comm)
        solver2.setOperators(A2)
        solver2.setType(PETSc.KSP.Type.GMRES)  # or "cg" or "preonly"
        # solver2.getPC().setType(PETSc.PC.Type.ILU)        # Incomplete LU
        solver2.getPC().setType(PETSc.PC.Type.HYPRE)      # Algebraic multigrid
        # PETSc.Options().setValue("pc_hypre_type", "boomeramg")
        solver2.setTolerances(rtol=1e-5, atol=1e-5, max_it=1000)
        # solver2.setMonitor(lambda ksp, its, rnorm: print(f"[KSP] Iteration {its}, residual = {rnorm:.3e}"))

        # Apply BCs and assemble
        A2.zeroEntries()
        dolfinx.fem.petsc.assemble_matrix(A2, a2, bcs=self.bc_mech)
        A2.assemble()

        with b2.localForm() as loc:
            loc.set(0)

        dolfinx.fem.petsc.assemble_vector(b2, L2)
        dolfinx.fem.petsc.apply_lifting(b2, [a2], [self.bc_mech])
        b2.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b2, self.bc_mech)

        solver2.solve(b2, self.displacement.x.petsc_vec)
        its = solver2.getIterationNumber()

        reason = solver2.getConvergedReason()
        its = solver2.getIterationNumber()

        print(f"\t[KSP] Converged in {its} iterations, reason = {reason}: {get_ksp_reason_name(reason)}")

        # # Linear problem definition
        # self.mech_problem = dolfinx.fem.petsc.LinearProblem(
        #     self.a_mech,
        #     self.L_mech,
        #     bcs = self.bc_mech,
        #     petsc_options = self.petsc_options_ilu,
        #     jit_options = self.jit_options,
        # )
        # self.displacement = self.mech_problem.solve()

        self.displacement.x.scatter_forward()

        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)
        self.stress_fuel = sigma(self.displacement, self.fuel["lmbda"], self.fuel["mu"]) - (3 * self.fuel["lmbda"] + 2 * self.fuel["mu"]) * self.fuel["alpha"] * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))
        self.stress_cladding = sigma(self.displacement, self.cladding["lmbda"], self.cladding["mu"]) - (3 * self.cladding["lmbda"] + 2 * self.cladding["mu"]) * self.cladding["alpha"] * (self.temperature - T_ref) * ufl.Identity(len(self.displacement))
        self.strain = epsilon(self.displacement)

        # max_displacement = np.max(np.abs(self.displacement.x.array))
        # print(f"\nMax displacement: {max_displacement:.3e} m")

        # u = self.displacement.x.array.reshape(self.mesh.geometry.x.shape)
        # disp_r = np.sqrt(u[:, 0]**2 + u[:, 1]**2)
        # print(f"Max radial displacement: {np.max(disp_r):.3e} m")
        # print(f"Average radial displacement at fuel outer: {np.mean(disp_r[self.geometry_manager.idx_fuel_outer]):.3e} m")

        return its

    def compute_deformed_radii(self):
        """
        Compute and print the deformed fuel outer radius, cladding inner radius, and cladding outer radius
        after solving the mechanical problem.
        """
        # Extract displacement field
        displacement = self.displacement.x.array.reshape((-1, 3))  # 3D problem
        mesh_geometry = self.mesh.geometry.x  # Get the original mesh coordinates

        # Compute original radial positions
        r_initial = np.sqrt(mesh_geometry[:, 0]**2 + mesh_geometry[:, 1]**2)  # Radial coordinate
        r_displacement = np.sqrt(displacement[:, 0]**2 + displacement[:, 1]**2)  # Radial displacement

        # Compute deformed radii
        r_deformed = r_initial + r_displacement

        # Get max/min values for specific regions
        fuel_outer_radius_deformed = np.max(r_deformed[np.isclose(r_initial, self.fuel_outer_radius, atol=1e-6)])
        cladding_inner_radius_deformed = np.min(r_deformed[np.isclose(r_initial, self.cladding_inner_radius, atol=1e-6)])
        cladding_outer_radius_deformed = np.max(r_deformed[np.isclose(r_initial, self.cladding_outer_radius, atol=1e-6)])

        # Print results
        print("\nDeformed Geometry:")
        print(f"Fuel Outer Radius (Deformed): {fuel_outer_radius_deformed:.6e} m")
        print(f"Cladding Inner Radius (Deformed): {cladding_inner_radius_deformed:.6e} m")
        print(f"Cladding Outer Radius (Deformed): {cladding_outer_radius_deformed:.6e} m")

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
