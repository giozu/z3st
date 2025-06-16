import dolfinx
import ufl
from petsc4py import PETSc
import numpy as np

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
        self.T = {}
        self.T_old = {}

    def set_thermal_boundary_conditions(self, V_t_sub, V_t_map=None):
        """
        Apply thermal boundary conditions for both staggered and mixed cases.

        If V_t_map is None → staggered (non-mixed spaces).
        If V_t_map is provided → mixed (collapsed subspace).

        Parameters:
            V_t_sub: FunctionSpace (collapsed temperature space)
            V_t_map: DoF map from mixed to collapsed (only in mixed case)
        """
        print(f"Setting thermal boundary conditions...")
        thermal_bcs_defs = self.boundary_conditions.get("thermal_bcs", {})

        dirichlet_bcs_list = []

        self.bcs_info.setdefault('neumann_thermal', [])
        self.bcs_info['neumann_thermal'].clear()

        for mat_type, bc_list in thermal_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                if region_name is None or bc_type is None:
                    print(f"  [WARNING] Incomplete thermal BC definition for '{mat_type}'. Skipping.")
                    continue
                
                region_id = self.label_map.get(region_name)
                if region_id is None:
                    print(f"  [WARNING] Region '{region_name}' not found in label_map for thermal BC. Skipping.")
                    continue

                if bc_type == "Dirichlet":
                    temperature = bc_info.get("temperature")
                    if temperature is None:
                        print(f"  [WARNING] Dirichlet BC on '{mat_type}' for region '{region_name}' has no temperature. Skipping.")
                        continue

                    T_d = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(temperature))

                    # locate dofs differently depending on mixed or staggered
                    if V_t_map is None:
                        # staggered mode (standard)
                        dofs = dolfinx.fem.locate_dofs_topological(V_t_sub, self.fdim, self.ft.find(region_id))
                    else:
                        # mixed mode
                        dofs_mixed = dolfinx.fem.locate_dofs_topological(
                            (self.W.sub(1), V_t_sub), self.fdim, self.ft.find(region_id)
                        )
                        dofs_mixed = np.array(dofs_mixed, dtype=np.int32)
                        dofs = V_t_map[dofs_mixed[:, 1]]

                    bc = dolfinx.fem.dirichletbc(T_d, dofs, V_t_sub)
                    dirichlet_bcs_list.append(bc)

                    print(f"  [INFO] Dirichlet thermal BC on '{mat_type}' → {temperature} K at region '{region_name}'")

                elif bc_type == "Neumann":
                    flux = bc_info.get("flux")
                    if flux is None:
                        print(f"  [WARNING] Neumann BC on '{mat_type}' for region '{region_name}' has no flux. Skipping.")
                        continue

                    self.bcs_info['neumann_thermal'].append({'id': region_id, 'value': float(flux)})
                    print(f"  [INFO] Neumann thermal BC on '{mat_type}' → {flux} W/m² at region '{region_name}'")

                else:
                    print(f"  [WARNING] Unknown thermal BC type '{bc_type}' for '{mat_type}'. Skipping.")

        return dirichlet_bcs_list

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