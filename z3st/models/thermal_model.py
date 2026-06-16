# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import sys

import dolfinx
import ufl
from mpi4py import MPI
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
        # Default linear-solver backend; user can override in input.yaml.
        self.th_cfg.setdefault("linear_solver", "iterative_hypre")

        print("[ThermalModel] options loaded from input.yaml:")
        for key, value in self.th_cfg.items():
            print(f"  {key:<20}: {value}")

    def set_thermal_boundary_conditions(self, V_t):
        """
        Apply thermal boundary conditions

        Parameters:
            V_t: FunctionSpace (temperature space)
        """
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

                    # Step-dependent temperature: a list of length n_steps gives
                    # a time-varying Dirichlet temperature. The Constant is updated per step by the
                    # solver; a scalar is broadcast to every step.
                    if isinstance(temperature, list):
                        # Step index is capped at the last entry by the solver,
                        # so a length mismatch with n_steps is tolerated (the
                        # final value simply holds). Warn if they differ.
                        if len(temperature) != self.n_steps:
                            print(
                                f"  [WARNING] Thermal Dirichlet list on '{label}' region "
                                f"'{region_name}' has length {len(temperature)} != n_steps "
                                f"{self.n_steps}; the last value will hold for extra steps."
                            )
                        raw_value = [float(t) for t in temperature]
                    else:
                        raw_value = [float(temperature)] * self.n_steps

                    T_d = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(raw_value[0]))

                    # locate dofs
                    dofs = dolfinx.fem.locate_dofs_topological(
                        V_t, self.fdim, self.facet_tags.find(region_id)
                    )
                    bc = dolfinx.fem.dirichletbc(T_d, dofs, V_t)
                    self.dirichlet_thermal[label].append(
                        {"id": region_id, "value": bc, "const": T_d, "raw": raw_value}
                    )

                    print(
                        f"  [INFO] Dirichlet thermal BC on '{label}' → {raw_value[0]} K "
                        f"(first step) at region '{region_name}'"
                    )

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
                    # Two types:
                    #   1. Gap conductance: pair + conductance from model → T_ext from another subdomain
                    #   2. Convective: h_conv + T_ext → fixed external bulk temperature
                    pair = bc_info.get("pair")
                    h_conv = bc_info.get("h_conv")
                    T_ext = bc_info.get("T_ext")

                    if pair is not None:
                        # Gap mode
                        self.robin_thermal[label].append({"id": region_id, "pair": pair})
                        print(f"  [INFO] Robin (gap) thermal BC on '{label}' at region '{region_name}' coupled with '{pair}'")
                    elif h_conv is not None and T_ext is not None:
                        # Convective mode
                        self.robin_thermal[label].append({
                            "id": region_id,
                            "h_conv": float(h_conv),
                            "T_ext": float(T_ext),
                        })
                        print(f"  [INFO] Robin (convective) thermal BC on '{label}' → h={h_conv} W/(m²·K), T_ext={T_ext} K at region '{region_name}'")
                    else:
                        print(f"  [ERROR] Robin BC on '{label}' requires either 'pair' (gap) or 'h_conv'+'T_ext' (convective).")
                        sys.exit(1)

                else:
                    print(f"  [ERROR] Unknown thermal BC type '{bc_type}' for '{label}'.")
                    print(f"  Available are: Dirichlet, Neumann, Robin.")
                    sys.exit(1)

    def heat_flux(self, T):
        """
        Compute the average heat flux (magnitude and per-component) per material.

        Dimension-aware: only the mesh's geometric components are assembled
        (a 2D mesh has no z-flux). Supports both scalar conductivity and
        symbolic k(T) material cards (where ``material["k"]`` is already a
        UFL expression resolved at field initialisation).
        """

        print("\n--- Average heat flux magnitude per material ---")
        gdim = self.mesh.geometry.dim
        comp_labels = ("r", "z") if self.regime == "axisymmetric" else ("x", "y", "z")[:gdim]

        for label, material in self.materials.items():
            tag = self.label_map[label]

            # Measure for integration over the specific material's subdomain
            dx = ufl.Measure(
                "dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag
            )

            k_val = material.get("k")
            if isinstance(k_val, (int, float)):
                k = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(k_val))
            else:
                # symbolic card: spine resolved k(T) to a UFL expression
                k = k_val

            q_vec = -k * ufl.grad(T)
            q_mag = ufl.sqrt(ufl.dot(q_vec, q_vec))

            # Directly assemble the integral of the UFL expression over the subdomain
            q_integral = dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_mag * dx))
            comp_integrals = [
                dolfinx.fem.assemble_scalar(dolfinx.fem.form(q_vec[i] * dx))
                for i in range(gdim)
            ]

            # Assemble the volume (or area in 2D) of the subdomain
            volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1.0 * dx))

            # assemble_scalar is per-rank partial; reduce to global sums
            # so the per-subdomain averages below are physical, not the
            # rank-0 share of the subdomain.
            comm = self.mesh.comm
            q_integral = comm.allreduce(q_integral, op=MPI.SUM)
            comp_integrals = [comm.allreduce(c, op=MPI.SUM) for c in comp_integrals]
            volume = comm.allreduce(volume, op=MPI.SUM)

            # Compute the average
            q_avg = q_integral / volume if volume > 0 else 0.0
            print(f"[INFO] Average heat flux magnitude in {label:<10}: {q_avg:.2f} W/m²")
            for name, integral in zip(comp_labels, comp_integrals):
                c_avg = integral / volume if volume > 0 else 0.0
                print(f"[INFO] Average heat flux -{name} in {label:<10}: {c_avg:.2f} W/m²")
