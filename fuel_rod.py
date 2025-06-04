import dolfinx.fem.petsc
import ufl
import dolfinx
import dolfinx.fem.petsc
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np

from finiteElementSettings import FiniteElementSetting
from geometryManager import GeometryManager

from models.thermal_model import ThermalSolver
from models.mechanical_model import MechanicalSolver

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
        self.geometry_manager = GeometryManager(
            mesh=self.mesh,
            surface_tags=self.ft,
            fuel_outer_label=self.fuel_outer_label,
            cladding_inner_label=self.cladding_inner_label,
            cladding_outer_label=self.cladding_outer_label,
            cladding_bottom_label=self.cladding_bottom_label,
            cladding_top_label=self.cladding_top_label,
            fuel_outer_radius=self.fuel_outer_radius,
            cladding_inner_radius=self.cladding_inner_radius,
            cladding_outer_radius=self.cladding_outer_radius,
            height=self.height,
            dofs_fuel_outer=self.dofs_fuel_outer,
            dofs_cladding_inner=self.dofs_cladding_inner,
            dofs_cladding_outer=self.dofs_cladding_outer,
            dofs_cladding_bottom=self.dofs_cladding_bottom,
            dofs_cladding_top=self.dofs_cladding_top
        )

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

        self.thermal_solver.update(self.lhr, self.fuel_outer_radius)

        # self.area_fuel = np.pi * self.fuel_outer_radius**2
        # self.q_third = self.lhr / self.area_fuel
        # self.heat_flux = self.lhr / (2.0 * np.pi * self.cladding_inner_radius)
        # self.gap_thickness = self.cladding_inner_radius - self.fuel_outer_radius

        # print(f"   fuel outer radius = {self.fuel_outer_radius:.3f} μm")
        # print(f"   cladding inner radius = {self.cladding_inner_radius:.3f} μm")
        # print(f"   cladding outer radius = {self.cladding_outer_radius:.3f} μm")
        # print(f"   fuel area = {self.area_fuel:.4e} m²")
        # print(f"   gap thickness = {self.gap_thickness*1e6:.3f} micron")
        # print(f"   q''' = {self.q_third:.2f} W/m³")
        # print(f"   heat flux at cladding inner surface = {self.heat_flux:.2f} W/m²")


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
        self.bc_mech_fuel = [
            dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_fuel_bottom, self.V_m)
        ]

        self.bc_mech_cladding = [
            dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_cladding_bottom, self.V_m)
        ]

        print(f"Cladding BCs: {len(self.bc_mech_cladding)}")
        assert self.bc_mech_cladding and len(self.bc_mech_cladding) > 0, "❌ No Dirichlet BCs for cladding!"

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
        print("🧪 Initializing temperature and displacement fields...")

        # --- Temperature ---
        T_init_fuel = float(self.T_bc_fo.value)
        T_init_clad = float(self.T_bc_co.value)

        self.T_fuel = dolfinx.fem.Function(self.V_t)
        self.T_fuel.name = "T_fuel"
        self.T_fuel.x.array[:] = T_init_fuel

        self.T_cladding = dolfinx.fem.Function(self.V_t)
        self.T_cladding.name = "T_cladding"
        self.T_cladding.x.array[:] = T_init_clad

        # --- Displacement ---
        self.u_fuel = dolfinx.fem.Function(self.V_m)
        self.u_fuel.name = "u_fuel"
        self.u_fuel.x.array[:] = 0.0

        self.u_cladding = dolfinx.fem.Function(self.V_m)
        self.u_cladding.name = "u_cladding"
        self.u_cladding.x.array[:] = 0.0

        # --- Derived quantities ---
        self.q_third = self.lhr / self.area_fuel
        self.heat_flux = self.lhr / (2.0 * np.pi * self.cladding_inner_radius)

        # --- Solver objects ---
        self.thermal_solver = ThermalSolver(
            mesh=self.mesh,
            V_t=self.V_t,
            dx_fuel=self.dx_fuel,
            dx_cladding=self.dx_cladding,
            ds_ci=self.ds_ci,
            dofs_fuel_outer=self.dofs_fuel_outer,
            dofs_cladding_outer=self.dofs_cladding_outer,
            dofs_cladding_inner=self.dofs_cladding_inner,
            T_bc_fo=self.T_bc_fo,
            T_bc_co=self.T_bc_co,
            lhr=self.lhr,
            q_third=self.q_third,
            heat_flux=self.heat_flux,
            fuel_props=self.fuel,
            cladding_props=self.cladding,
            gap_thickness=self.gap_thickness,
            fuel_outer_radius=self.fuel_outer_radius,
            filling_gas=self.filling_gas,
            T_fuel=self.T_fuel,
            T_cladding=self.T_cladding,
            debug=getattr(self, "debug_mode", False)
        )

        self.mechanical_solver = MechanicalSolver(
            mesh=self.mesh,
            V_m=self.V_m,
            dx_fuel=self.dx_fuel,
            dx_cladding=self.dx_cladding,
            ds_fo=self.ds_fo,
            ds_ci=self.ds_ci,
            ds_co=self.ds_co,
            dofs_fuel_bottom=self.dofs_fuel_bottom,
            dofs_cladding_bottom=self.dofs_cladding_bottom,
            fuel_props=self.fuel,
            cladding_props=self.cladding,
            g=self.g,
            T_ref=self.T_ref,
            bc_mech_fuel=self.bc_mech_fuel,
            bc_mech_cladding=self.bc_mech_cladding,
            u_fuel=self.u_fuel,
            u_cladding=self.u_cladding,
            T_fuel=self.T_fuel,
            T_cladding=self.T_cladding,
            gap_traction=self.gap_traction,
            coolant_traction=self.coolant_traction,
            debug=getattr(self, "debug_mode", False)
        )

        self.total_thermal_iterations = 0
        self.total_mech_iterations = 0

    def solve(self, max_iter=10, tol=1e-2):
        """
        Iterative thermo-mechanical steady-state solver with fixed mesh.
        The solver alternates between thermal and mechanical sub-solvers, applying 
        under-relaxation to temperature and displacement fields at each iteration.

        Convergence is reached when the relative change in maximum radial displacement 
        (in the fuel) is below the specified tolerance.

        Parameters:
            max_iter (int): Maximum number of iterations.
            tol (float): Relative convergence tolerance on max radial displacement.

        Returns:
            bool: True if converged, False otherwise.
        """

        print("======================")
        print("Thermo-mechanical loop")
        print("======================")

        converged = False

        u_r_max_prev = None
        gap_mean_prev = None
        
        omega = 0.7  # Relaxation factor

        print(f"🔁 Using relaxation factor ω = {omega}")

        # Store initial values for relaxation
        T_fuel_prev     = np.copy(self.T_fuel.x.array)
        T_cladding_prev = np.copy(self.T_cladding.x.array)
        u_fuel_prev     = np.copy(self.u_fuel.x.array)
        u_cladding_prev = np.copy(self.u_cladding.x.array)

        for it in range(1, max_iter + 1):
            print(f"\n🔁 Thermo-mechanical iteration {it}")

            it_therm = self.thermal_solver.thermal_solver(self.u_t, self.v_t)
            
            T_vals = self.T_fuel.x.array
            
            if getattr(self, "debug_mode"):
                print(f"[DEBUG] T_fuel min/max = {T_vals.min():.2f} / {T_vals.max():.2f}")

            it_mech = self.mechanical_solver.mechanical_solver(self.u_m, self.v_m)

            # --- Relaxation step ---
            self.T_fuel.x.array[:]     = (1 - omega) * T_fuel_prev     + omega * self.T_fuel.x.array
            self.T_cladding.x.array[:] = (1 - omega) * T_cladding_prev + omega * self.T_cladding.x.array
            self.u_fuel.x.array[:]     = (1 - omega) * u_fuel_prev     + omega * self.u_fuel.x.array
            self.u_cladding.x.array[:] = (1 - omega) * u_cladding_prev + omega * self.u_cladding.x.array

            print(f"\n🔁 Relaxation applied with ω = {omega}")

            # Store new values for next iteration
            T_fuel_prev     = np.copy(self.T_fuel.x.array)
            T_cladding_prev = np.copy(self.T_cladding.x.array)
            u_fuel_prev     = np.copy(self.u_fuel.x.array)
            u_cladding_prev = np.copy(self.u_cladding.x.array)

            # Compute max radial displacement
            u_fuel_array = self.u_fuel.x.array.reshape(self.mesh.geometry.x.shape)
            disp_r = np.sqrt(u_fuel_array[:, 0]**2 + u_fuel_array[:, 1]**2)
            u_r_max = np.max(disp_r)

            (gap_mean, gap_min, gap_max,
            fuel_outer_mean, fuel_outer_min, fuel_outer_max,
            clad_inner_mean, clad_inner_min, clad_inner_max,
            clad_outer_mean, clad_outer_min, clad_outer_max) = self.geometry_manager.compute_deformed_geometry(self.u_fuel, self.u_cladding)

            print(f"   Gap mean = {gap_mean*1e6:.2f} μm, min = {gap_min*1e6:.2f} μm, max = {gap_max*1e6:.2f} μm")
            print(f"   Fuel outer radius (mean/min/max) = {fuel_outer_mean*1e6:.3f} / {fuel_outer_min*1e6:.3f} / {fuel_outer_max*1e6:.3f} μm")
            print(f"   Clad inner radius (mean/min/max) = {clad_inner_mean*1e6:.3f} / {clad_inner_min*1e6:.3f} / {clad_inner_max*1e6:.3f} μm")
            print(f"   Clad outer radius (mean/min/max) = {clad_outer_mean*1e6:.3f} / {clad_outer_min*1e6:.3f} / {clad_outer_max*1e6:.3f} μm")

            # --- Convergence check ---
            if u_r_max_prev is not None and gap_mean_prev is not None:
                delta_u = abs(u_r_max - u_r_max_prev) / max(abs(u_r_max_prev), 1e-12)
                delta_gap = abs(gap_mean - gap_mean_prev) / max(abs(gap_mean_prev), 1e-12)

                print(f"\tΔu_r_rel = {delta_u:.2e}")
                print(f"\tΔgap_rel = {delta_gap:.2e}")

                if delta_u < tol and delta_gap < tol:
                    print(f">> Thermo-mechanical convergence achieved.")
                    print(f">> Iterations: thermal = {it_therm}, mechanical = {it_mech}")

                    self.geometry_manager.geometry_summary(self.u_fuel, self.u_cladding)
                    self.geometry_manager.update_geometry(self.u_fuel, self.u_cladding)

                    self.total_thermal_iterations += it_therm
                    self.total_mech_iterations += it_mech

                    converged = True
                    break
            else:
                print(f"\tInitial u_r_max = {u_r_max*1e6:.2f} μm")
                print(f"\tInitial gap_mean = {gap_mean*1e6:.2f} μm")

            u_r_max_prev = u_r_max
            gap_mean_prev = gap_mean

        if not converged:
            print("⚠️ Warning: convergence not reached in max_iter.")

        return converged
    
    def assemble_thermal_system(self):
        """
        Solve the thermal problem separately in cladding and fuel domains.
        1. Solve cladding temperature field.
        2. Compute gap conductance and update T_bc_fo.
        3. Solve fuel temperature field.
        """
        # --- STEP 1: Cladding system ---
        print("\n⚙️  Assembling and solving cladding thermal system")
        a_clad, L_clad = self.assemble_cladding_system()
        its_clad = self.solve_thermal_problem(a_clad, L_clad, self.T_cladding, bcs=[self.bc_Tco])

        # --- STEP 2: Compute ΔT_gap ---
        T_clad_inner_values = self.T_cladding.x.array[self.dofs_cladding_inner]
        T_clad_inner_avg = np.mean(T_clad_inner_values)
        T_clad_inner_min = np.min(T_clad_inner_values)
        T_clad_inner_max = np.max(T_clad_inner_values)
        delta_T_clad_inner = T_clad_inner_max - T_clad_inner_min
        max_dev_from_avg = np.max(np.abs(T_clad_inner_values - T_clad_inner_avg))

        print(f"📊 T_clad @ inner surface → avg = {T_clad_inner_avg:.2f} K, "
              f"min = {T_clad_inner_min:.2f} K, max = {T_clad_inner_max:.2f} K, "
              f"ΔT = {delta_T_clad_inner:.2f} K, max|ΔT - avg| = {max_dev_from_avg:.2f} K")

        T_fuel_outer_bc = float(self.T_bc_fo.value)
        T_gap_guess = 0.5 * (T_clad_inner_avg + T_fuel_outer_bc)

        print(f"\n🔁 Estimating gap conductance from T_cladding_inner = {T_clad_inner_avg:.2f} K")
        delta_T_gap = self.compute_gap_conductance(T_gap_guess)

        # Update boundary condition on fuel outer surface
        self.T_bc_fo = dolfinx.fem.Constant(self.mesh, T_clad_inner_avg + delta_T_gap)
        self.bc_Tfo = dolfinx.fem.dirichletbc(self.T_bc_fo, self.dofs_fuel_outer, self.V_t)

        print(f"🔧 Updated T_bc_fo = {float(self.T_bc_fo.value):.2f} K (ΔT_gap = {delta_T_gap:.2f} K)")

        # --- STEP 3: Fuel system ---
        print("\n⚙️  Assembling and solving fuel thermal system")
        a_fuel, L_fuel = self.assemble_fuel_system()
        its_fuel = self.solve_thermal_problem(a_fuel, L_fuel, self.T_fuel, bcs=[self.bc_Tfo])

        return its_clad, its_fuel

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
        print(f"Fuel Outer Radius (Deformed): {fuel_outer_radius_deformed:.3f} μm")
        print(f"Cladding Inner Radius (Deformed): {cladding_inner_radius_deformed:.3f} μm")
        print(f"Cladding Outer Radius (Deformed): {cladding_outer_radius_deformed:.3f} μm")

    def get_results(self):
        """
        Compute and return displacement, strain, and stress fields for fuel and cladding.

        Returns:
            u_fuel (Function): Displacement field in the fuel.
            u_cladding (Function): Displacement field in the cladding.
            strain_fuel (ufl.Expr): Strain tensor in the fuel.
            strain_cladding (ufl.Expr): Strain tensor in the cladding.
            stress_fuel (ufl.Expr): Stress tensor in the fuel (with thermal expansion).
            stress_cladding (ufl.Expr): Stress tensor in the cladding (with thermal expansion).
            T_fuel (Function): Temperature field in the fuel.
            T_cladding (Function): Temperature field in the cladding.
        """

        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u, lmbda, mu):
            return lmbda * ufl.tr(epsilon(u)) * ufl.Identity(3) + 2 * mu * epsilon(u)

        T_ref = dolfinx.fem.Constant(self.mesh, self.T_ref)

        # Strain
        strain_fuel = epsilon(self.u_fuel)
        strain_cladding = epsilon(self.u_cladding)

        # Stress (contributo termico incluso)
        stress_fuel = (
            sigma(self.u_fuel, self.fuel["lmbda"], self.fuel["mu"]) -
            (3 * self.fuel["lmbda"] + 2 * self.fuel["mu"]) * self.fuel["alpha"] *
            (self.T_fuel - T_ref) * ufl.Identity(3)
        )

        stress_cladding = (
            sigma(self.u_cladding, self.cladding["lmbda"], self.cladding["mu"]) -
            (3 * self.cladding["lmbda"] + 2 * self.cladding["mu"]) * self.cladding["alpha"] *
            (self.T_cladding - T_ref) * ufl.Identity(3)
        )

        return (
            self.u_fuel,
            self.u_cladding,
            strain_fuel,
            strain_cladding,
            stress_fuel,
            stress_cladding,
            self.T_fuel,
            self.T_cladding
        )

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
