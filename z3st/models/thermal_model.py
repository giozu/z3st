# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import sys

import dolfinx
import ufl
from petsc4py import PETSc


class ThermalModel:
    def __init__(self):
        print("[ThermalModel] initializer")
        self.a_thermal = {}
        self.L_thermal = {}
        self.thermal_residuals = {}
        self.neumann_thermal = {}
        self.dirichlet_thermal = {}
        self.robin_thermal = {}

        # --. Thermal model options --..
        self.th_cfg = self.input_file.get("thermal", {})

        if not self.th_cfg:
            raise ValueError("[ThermalModel] 'thermal' missing in input.yaml.")

        print("[ThermalModel] options loaded from input.yaml:")
        for key, value in self.th_cfg.items():
            print(f"  {key:<20}: {value}")

    def set_thermal_boundary_conditions(self, V_t, V_t_map=None):
        """
        Apply thermal boundary conditions for both staggered and mixed cases.

        If V_t_map is None → staggered (non-mixed spaces).
        If V_t_map is provided → mixed (collapsed subspace).

        Parameters:
            V_t_sub: FunctionSpace (collapsed temperature space)
            V_t_map: DoF map from mixed to collapsed (only in mixed case)
        """
        print(f"\nSetting thermal boundary conditions...")
        thermal_bcs_defs = self.boundary_conditions.get("thermal", {})

        seen_regions = {}

        for label, bc_list in thermal_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                key = (region_name, bc_type)
                if key in seen_regions:
                    print(
                        f"[WARNING] Duplicate thermal BC of type '{bc_type}' defined for region '{region_name}' (previously in '{seen_regions[key]}', now in '{label}')."
                    )
                else:
                    seen_regions[key] = label

        for label in self.materials:
            self.neumann_thermal[label] = []
            self.dirichlet_thermal[label] = []
            self.robin_thermal[label] = []

        for label, bc_list in thermal_bcs_defs.items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                bc_type = bc_info.get("type")

                if region_name is None or bc_type is None:
                    print(f"  [ERROR] Incomplete thermal BC definition for '{label}'.")
                    sys.exit(1)

                region_id = self.label_map.get(region_name)
                if region_id is None:
                    print(
                        f"  [ERROR] Region '{region_name}' not found in label_map for thermal BC."
                    )
                    sys.exit(1)

                if bc_type == "Dirichlet":
                    temperature = bc_info.get("temperature")
                    if temperature is None:
                        print(
                            f"  [ERROR] Dirichlet BC on '{label}' for region '{region_name}' has no temperature."
                        )
                        sys.exit(1)

                    T_d = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(temperature))

                    # locate dofs
                    dofs = dolfinx.fem.locate_dofs_topological(
                        V_t, self.fdim, self.facet_tags.find(region_id)
                    )
                    bc = dolfinx.fem.dirichletbc(T_d, dofs, V_t)
                    self.dirichlet_thermal[label].append(bc)

                    print(
                        f"  [INFO] Dirichlet thermal BC on '{label}' → {temperature} K at region '{region_name}'"
                    )
                    # print(f"  {type(bc)}")
                    # print("  Applied dofs:", bc.dof_indices)

                elif bc_type == "Neumann":
                    flux = bc_info.get("flux")
                    if flux is None:
                        print(
                            f"  [ERROR] Neumann BC on '{label}' for region '{region_name}' has no flux."
                        )
                        sys.exit(1)

                    self.neumann_thermal[label].append({"id": region_id, "value": float(flux)})
                    print(
                        f"  [INFO] Neumann thermal BC on '{label}' → {flux} W/m² at region '{region_name}'"
                    )

                elif bc_type == "Robin":

                    pair = bc_info.get("pair")
                    if pair is None:
                        print(f"  [ERROR] Robin BC on '{label}' is missing 'pair'")
                        sys.exit(1)
                    self.robin_thermal[label].append({"id": region_id, "pair": pair})

                    print(
                        f"  [INFO] Robin thermal BC on '{label}' at region '{region_name}' coupled with '{pair}'"
                    )

                else:
                    print(f"  [ERROR] Unknown thermal BC type '{bc_type}' for '{label}'.")
                    print(f"  Available are: Dirichlet, Neumann.")
                    sys.exit(1)

    def heat_flux(self, T):
        """
        Compute average heat flux (magnitude and x-y-z components) per materials

        """

        print("\n--- Average heat flux magnitude per material ---")
        for label, material in self.materials.items():
            tag = self.label_map[label]

            # Measure for integration over the specific material's subdomain
            dx = ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag
            )

            k = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(material["k"]))

            grad_T = ufl.grad(T)

            q_vec = -k * grad_T

            qx = q_vec[0]
            qy = q_vec[1]
            qz = q_vec[2]

            q_mag = ufl.sqrt(ufl.dot(q_vec, q_vec))

            # Directly assemble the integral of the UFL expression over the subdomain
            q_integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_mag * dx))
            qx_integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(qx * dx))
            qy_integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(qy * dx))
            qz_integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(qz * dx))

            # Assemble the volume (or area in 2D) of the subdomain
            volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1.0 * dx))

            # Compute the average
            q_avg = q_integral / volume if volume > 0 else 0.0
            qx_avg = qx_integral / volume if volume > 0 else 0.0
            qy_avg = qy_integral / volume if volume > 0 else 0.0
            qz_avg = qz_integral / volume if volume > 0 else 0.0

            print(f"[INFO] Average heat flux magnitude in {label:<10}: {q_avg:.2f} W/m²")

            print(f"[INFO] Average heat flux -x in {label:<10}: {qx_avg:.2f} W/m²")
            print(f"[INFO] Average heat flux -y in {label:<10}: {qy_avg:.2f} W/m²")
            print(f"[INFO] Average heat flux -z in {label:<10}: {qz_avg:.2f} W/m²")
