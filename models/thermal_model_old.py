import numpy as np
import dolfinx
import ufl
from petsc4py import PETSc


class ThermalSolver:
    def __init__(self, mesh, V_t, dx_fuel, dx_cladding, ds_ci, dofs_fuel_outer, dofs_cladding_outer, dofs_cladding_inner,
                 T_bc_fo, T_bc_co, lhr, q_third, heat_flux, fuel_props, cladding_props,
                 gap_thickness, fuel_outer_radius, filling_gas, T_fuel, T_cladding, debug=False):
        self.mesh = mesh
        self.V_t = V_t
        self.dx_fuel = dx_fuel
        self.dx_cladding = dx_cladding
        self.ds_ci = ds_ci
        self.dofs_fuel_outer = dofs_fuel_outer
        self.dofs_cladding_outer = dofs_cladding_outer
        self.dofs_cladding_inner = dofs_cladding_inner
        self.T_bc_fo = T_bc_fo
        self.T_bc_co = T_bc_co
        self.lhr = lhr
        self.q_third = q_third
        self.heat_flux = heat_flux
        self.fuel = fuel_props
        self.cladding = cladding_props
        self.gap_thickness = gap_thickness
        self.fuel_outer_radius = fuel_outer_radius
        self.filling_gas = filling_gas
        self.T_fuel = T_fuel
        self.T_cladding = T_cladding
        self.debug = debug

        if self.debug:
            print("[DEBUG] ThermalSolver.__init__")
            print(f"  Linear Heat Rate (lhr): {self.lhr} W/m")
            print(f"  Volumetric Heat Source (q_third): {self.q_third}")
            print(f"  Heat Flux: {self.heat_flux}")
            print(f"  Fuel Properties: {self.fuel}")
            print(f"  Cladding Properties: {self.cladding}")
            print(f"  Gap Thickness: {self.gap_thickness} m")
            print(f"  Fuel Outer Radius: {self.fuel_outer_radius} m")
            print(f"  Filling Gas Properties: {self.filling_gas}")
            print(f"  Dirichlet DOFs - Fuel Outer: {len(self.dofs_fuel_outer)}")
            print(f"  Dirichlet DOFs - Cladding Outer: {len(self.dofs_cladding_outer)}")
            print(f"  Dirichlet DOFs - Cladding Inner: {len(self.dofs_cladding_inner)}")
            print(f"  Initial BC T_bc_fo: {float(self.T_bc_fo.value)}")
            print(f"  Initial BC T_bc_co: {float(self.T_bc_co.value)}")

    def update(self, lhr: float, fuel_outer_radius: float):
        self.lhr = lhr
        self.fuel_outer_radius = fuel_outer_radius
        self.area_fuel = np.pi * fuel_outer_radius**2
        self.q_third = lhr / self.area_fuel
        self.heat_flux = lhr / (2.0 * np.pi * fuel_outer_radius)

        if self.debug:
            print("[DEBUG] ThermalSolver.update:")
            print(f"  LHR = {self.lhr*1e-3:.2f} kW/m")
            print(f"  Fuel outer radius = {fuel_outer_radius*1e3:.2f} mm")
            print(f"  Area = {self.area_fuel:.2e} m²")
            print(f"  q''' = {self.q_third:.2f} W/m³")
            print(f"  Heat flux = {self.heat_flux:.2f} W/m²")

    def assemble_fuel_system(self, u_t, v_t, T_fuel):
        k_fuel_expr = 1.0 / (0.0452 + 2.46e-4 * T_fuel)

        a_fuel = ufl.inner(k_fuel_expr * ufl.grad(u_t), ufl.grad(v_t)) * self.dx_fuel
        q_vol = dolfinx.fem.Constant(self.mesh, self.q_third)
        L_fuel = q_vol * v_t * self.dx_fuel

        if self.debug:
            print("[DEBUG] ThermalSolver.assemble_fuel_system")
    
            T_vals = T_fuel.x.array
            k_vals = 1.0 / (0.0452 + 2.46e-4 * T_vals)
    
            print("\n[DEBUG] Fuel conductivity:")
            print(f"  k_min = {np.min(k_vals):.4f}, k_max = {np.max(k_vals):.4f}, k_mean = {np.mean(k_vals):.4f}")

            L_fuel_form = dolfinx.fem.form(L_fuel)
            rhs_norm = dolfinx.fem.assemble_scalar(L_fuel_form)
            print("[DEBUG] q_third:", self.q_third)
            print("[DEBUG] Fuel RHS norm:", rhs_norm)

        return a_fuel, L_fuel

    def assemble_cladding_system(self, u_t, v_t):
        k_clad = self.cladding["k"]
        a_clad = ufl.inner(k_clad * ufl.grad(u_t), ufl.grad(v_t)) * self.dx_cladding
        heat_flux = dolfinx.fem.Constant(self.mesh, self.heat_flux)
        L_clad = heat_flux * v_t * self.ds_ci
        return a_clad, L_clad

    def solve_thermal_problem(self, a, L, T_field, bcs):
        a_form = dolfinx.fem.form(a)
        L_form = dolfinx.fem.form(L)

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

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        solver.setType(PETSc.KSP.Type.GMRES)
        solver.getPC().setType("jacobi")
        solver.setTolerances(rtol=1e-6, atol=1e-6, max_it=1000)

        solver.solve(b, T_field.x.petsc_vec)
        T_field.x.scatter_forward()

        return solver.getIterationNumber()

    def compute_gap_conductance(self, T_gap_initial):
        if self.lhr < 1e-6:
            if self.debug:
                print("[DEBUG] Linear Heat Rate too low:", self.lhr)
            return 0.0
        
        if self.debug:
            print("[DEBUG] Starting gap conductance calculation")
            print(f"[DEBUG] Inputs - T_gap_initial: {T_gap_initial}, gap_thickness: {self.gap_thickness}, "
                f"fuel_outer_radius: {self.fuel_outer_radius}, gas A: {self.filling_gas['A']}")

        T_gap = T_gap_initial
        delta_T_old = 1.0
        tol = 1e-3
        error = 1.0
        max_iter = 20
        iteration = 0

        while error > tol and iteration < max_iter:
            k_gas = self.filling_gas['A'] * T_gap**0.79
            h_gap = k_gas / self.gap_thickness
            R_gap = 1.0 / (2.0 * np.pi * self.fuel_outer_radius * h_gap)
            delta_T_new = self.lhr * R_gap
            error = abs(delta_T_new - delta_T_old) / max(abs(delta_T_old), 1e-6)

            if self.debug:
                print(f"[DEBUG] Iter {iteration}: T_gap = {T_gap:.4f}, k_gas = {k_gas:.4f}, "
                f"h_gap = {h_gap:.4f}, R_gap = {R_gap:.6f}, delta_T = {delta_T_new:.4f}, error = {error:.4e}")

            delta_T_old = delta_T_new
            T_gap = T_gap_initial + delta_T_new
            iteration += 1

        if iteration == max_iter:
            print("[WARNING] Gap conductance iteration reached max iterations")

        return delta_T_new

    def thermal_solver(self, u_t, v_t):
        # STEP 1: Solve cladding
        a_clad, L_clad = self.assemble_cladding_system(u_t, v_t)
        its_clad = self.solve_thermal_problem(a_clad, L_clad, self.T_cladding, [
            dolfinx.fem.dirichletbc(self.T_bc_co, self.dofs_cladding_outer, self.V_t)
        ])

        # STEP 2: Estimate gap BC
        T_clad_inner_values = self.T_cladding.x.array[self.dofs_cladding_inner]
        T_clad_inner_avg = np.mean(T_clad_inner_values)
        T_gap_guess = 0.5 * (T_clad_inner_avg + float(self.T_bc_fo.value))
        delta_T_gap = self.compute_gap_conductance(T_gap_guess)
        new_T_bc_fo = T_clad_inner_avg + delta_T_gap

        # Update BC and STEP 3: Solve fuel
        self.T_bc_fo = dolfinx.fem.Constant(self.mesh, new_T_bc_fo)
        a_fuel, L_fuel = self.assemble_fuel_system(u_t, v_t, self.T_fuel)

        if self.debug:
            # print("[DEBUG] Dirichlet DOFs fuel outer:", len(self.dofs_fuel_outer))
            print(f"[DEBUG] Temperature jump across the gap = {delta_T_gap} K")
            print("[DEBUG] dofs_cladding_inner length:", len(self.dofs_cladding_inner))
            print("[DEBUG] T_clad_inner_values:", T_clad_inner_values)
            print("[DEBUG] Initial T_bc_fo:", float(self.T_bc_fo.value))

        its_fuel = self.solve_thermal_problem(a_fuel, L_fuel, self.T_fuel, [
            dolfinx.fem.dirichletbc(self.T_bc_fo, self.dofs_fuel_outer, self.V_t)
        ])

        return its_clad + its_fuel
