import dolfinx
import ufl
from petsc4py import PETSc

class ThermalModel:
    def __init__(self):
        print("__ThermalModel initializer__")
        self.a_thermal = {}
        self.L_thermal = {}
        self.thermal_residuals = {}
        self.q_third = {}
        self.bcs_thermal = {}
        self.neumann = {}
        self.neumann_label = {}
        self.i_th = {}

    def set_thermal_boundary_conditions(self):
        print(f"Thermal boundary conditions")
        self.thermal_bcs = self.boundary_conditions.get("thermal_bcs", {})

        for mat_type, _ in self.materials.items():
            bc_list = self.thermal_bcs.get(mat_type, [])

            self.bcs_thermal[mat_type] = []

            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                if region_name is None or bc_type is None:
                    print(f"  [WARNING] Incomplete BC definition for '{mat_type}'")
                    continue

                if bc_type == "Dirichlet":
                    temperature = bc_info.get("temperature")
                    if temperature is None:
                        print(f"  [WARNING] Dirichlet BC on '{mat_type}' but no temperature specified.")
                        continue

                    region_id = self.label_map.get(region_name)
                    if region_id is None:
                        print(f"  [WARNING] Region '{region_name}' not found in label_map for Dirichlet.")
                        continue

                    dofs = self.locateThermalDofs(region_id)
                    self.bcs_thermal[mat_type].append(
                        self.fixed_temperature(self.mesh, dofs, self.V_t, dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(temperature)))
                    )
                    print(f"  [INFO] Dirichlet BC on '{mat_type}' → {temperature} K at region '{region_name}'")

                elif bc_type == "Neumann":
                    flux = bc_info.get("flux")
                    if flux is None:
                        print(f"[WARNING] Neumann BC on '{mat_type}' but no flux specified.")
                        continue

                    region_id = self.label_map.get(region_name)
                    if region_id is None:
                        print(f"[WARNING] Region '{region_name}' not found in label_map for Neumann.")
                        continue

                    ds_local = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=region_id)

                    flux_constant = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(flux))

                    if mat_type not in self.neumann:
                        self.neumann[mat_type] = flux_constant * self.v_t * ds_local
                    else:
                        self.neumann[mat_type] += flux_constant * self.v_t * ds_local

                    print(f"  [INFO] Neumann BC on '{mat_type}' → {flux} W/m² at region '{region_name}'")

                else:
                    print(f"  [WARNING] Unknown BC type '{bc_type}' for '{mat_type}'")

    def fixed_temperature(self, mesh, dofs, V, T):
        """
        Create a Dirichlet BC with fixed temperature.

        Parameters:
            mesh (dolfinx.Mesh): The mesh object.
            dofs (array-like): Degrees of freedom where BC is applied.
            V (FunctionSpace): The function space.
            T (float): Temperature value.

        Returns:
            dolfinx.fem.DirichletBC: Fixed temperature boundary condition.
        """
        return dolfinx.fem.dirichletbc(
            dolfinx.fem.Constant(mesh, PETSc.ScalarType(T)), dofs, V
        )

    def assemble_thermal_problem(self, material_label):
        
        k = self.materials[material_label]["k"]
        h_gap = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(0.5))  # gap conductance
        q_third = PETSc.ScalarType(self.q_third[material_label])

        if self.debug_mode is True:
            print(f"[DEBUG] Thermal conductivity k = {k} W/m/K")
            print(f"[DEBUG] Gap conductance h_gap = {h_gap.value} W/m²/K")
            print(f"[DEBUG] Volumetric heat source q_third = {q_third} W/m³")

        dx_local = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags, subdomain_id=self.label_map[material_label])

        a = ufl.dot(k * ufl.grad(self.u_t), ufl.grad(self.v_t)) * dx_local
        L = dolfinx.fem.Constant(self.mesh, q_third) * self.v_t * dx_local

        # Neumann BCs if present
        if material_label in self.neumann:
            L += self.neumann[material_label]

        # --- Coupling interface term ---
        # --- Robin interface weak coupling ---
        if material_label == "fuel":
            T_other = self.T["cladding"]
            subdomain_id = self.label_map.get("face_xx_1_1")
            ds_interface = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=subdomain_id)
            a += h_gap * ufl.inner(self.u_t, self.v_t) * ds_interface
            L += h_gap * T_other * self.v_t * ds_interface
            # L -= h_gap * (self.T["fuel"] - T_other) * self.v_t * ds_interface
        elif material_label == "cladding":
            T_other = self.T["fuel"]
            subdomain_id = self.label_map.get("face_xx_0_2")
            ds_interface = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=subdomain_id)
            a -= h_gap * ufl.inner(self.u_t, self.v_t) * ds_interface
            L -= h_gap * T_other * self.v_t * ds_interface
            # L -= h_gap * (T_other - self.T["cladding"]) * self.v_t * ds_interface

        # Store assembled forms
        self.a_thermal[material_label] = a
        self.L_thermal[material_label] = L

    def solve_thermal_problem(self, material_label, max_iter=1000, tol=1e-8):

        if self.debug_mode is True:
            print(f"[DEBUG] self.a_thermal type: {type(self.a_thermal[material_label])}")
            print(f"[DEBUG] integrands: {self.a_thermal[material_label].integrals()}")

        a_form = dolfinx.fem.form(self.a_thermal[material_label])
        L_form = dolfinx.fem.form(self.L_thermal[material_label])

        A = dolfinx.fem.petsc.create_matrix(a_form)
        b = dolfinx.fem.petsc.create_vector(L_form)

        dolfinx.fem.petsc.assemble_matrix(A, a_form, bcs=self.bcs_thermal[material_label])
        A.assemble()

        with b.localForm() as loc:
            loc.set(0.0)
        dolfinx.fem.petsc.assemble_vector(b, L_form)
        dolfinx.fem.petsc.apply_lifting(b, [a_form], [self.bcs_thermal[material_label]])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        dolfinx.fem.petsc.set_bc(b, self.bcs_thermal[material_label])

        solver = PETSc.KSP().create(self.mesh.comm)
        solver.setOperators(A)
        # solver.setType(PETSc.KSP.Type.GMRES)
        # solver.getPC().setType("jacobi")
        
        # solver.setType("lu")
        
        solver.setType("gmres")
        solver.getPC().setType("gamg")
        # solver.getPC().setType("hypre")

        solver.setTolerances(rtol=tol, atol=tol, max_it=max_iter)

        residuals = []

        def monitor(ksp, its, rnorm):
            residuals.append(rnorm)

        solver.setMonitor(monitor)

        if self.debug_mode is True:
            print(f"[DEBUG] Solving thermal problem for {material_label}")
            print(f"[DEBUG] Initial temperature: min = {self.T[material_label].x.array.min():.2f} K, max = {self.T[material_label].x.array.max():.2f} K")
            print(f"[DEBUG] Source term: q_third = {self.q_third[material_label]:.2f} W/m³")

        solver.solve(b, self.T[material_label].x.petsc_vec)
        self.T[material_label].x.scatter_forward()
        self.thermal_residuals[material_label] = residuals

        return solver.getIterationNumber()
