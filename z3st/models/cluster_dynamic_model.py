# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc


class ClusterDynamicsModel:
    """
    Cluster dynamics model for defect evolution.

    This module implements a 1D cluster dynamics model for simulating
    the evolution of defect cluster size distributions.

    Solves the diffusion equation for cluster density c(n,t) where n is cluster size.
    """
    
    def __init__(self):
        """
        Initialize the cluster dynamics model.
        Uses self.mesh and self.V_c from the parent Spine class.
        
        Physical model:
        ∂c/∂t = -v ∂c/∂n + D ∂²c/∂n²

        where:
        - v is the advection velocity (cluster.advection_velocity, default 1.0)
        - D is the diffusion coefficient (cluster.diffusion_coefficient, default 0.5)
        - c(n,t) is cluster density at size n and time t
        """
        self.cluster_cfg = self.input_file.get("cluster", {})

        # Physical parameters for advection-diffusion equation
        self.v_cluster = self.cluster_cfg.get("advection_velocity", 1.0)  # Advection velocity (cluster growth rate)
        self.D_cluster = self.cluster_cfg.get("diffusion_coefficient", 0.5)  # Diffusion coefficient (cluster size fluctuations)

        print("[ClusterDynamicsModel] Initialized")
        print(f"  [ClusterDynamicsModel] Advection velocity: {self.v_cluster}")
        print(f"  [ClusterDynamicsModel] Diffusion coefficient: {self.D_cluster}")

    def set_cluster_initial_conditions(self):
        """
        Sets the initial condition for the cluster density field.
        Supported types: 'constant', 'gaussian'.
        """
        print("\n[ClusterDynamics] Setting initial distribution...")
        
        # Load configuration from input.yaml
        ic_config = self.cluster_cfg.get("initial_condition", {})
        ic_type = ic_config.get("type", "constant")

        # Constant initial condition
        if ic_type == "constant":
            val = ic_config.get("value")
            region_name = ic_config.get("region")
            
            if val is not None:
                if region_name is not None:
                    # Apply to a specific named region (defined in geometry.yaml)
                    region_id = self.label_map.get(region_name)
                    if region_id is not None:
                        facets = self.facet_tags.find(region_id)
                        dofs = dolfinx.fem.locate_dofs_topological(self.V_c, self.fdim, facets)
                        
                        if len(dofs) > 0:
                            self.c.x.array[dofs] = float(val)
                            self.c_n.x.array[dofs] = float(val)
                            print(f"  [INFO] Cluster IC applied: c = {val} on region '{region_name}' ({len(dofs)} DOFs)")
                        else:
                            print(f"  [WARNING] No DOFs found for region '{region_name}'")
                    else:
                        print(f"  [ERROR] Region '{region_name}' not found in label_map.")
                else:
                    print(f"  [ERROR] No region specified for constant IC.")

                # Mass conservation target initialization
                # Calculate C_tot = ∫ c(n,0) * n dn to be conserved during the simulation
                try:
                    x_coord = ufl.SpatialCoordinate(self.mesh)
                    n_coord = x_coord[0]
                    flux_form = dolfinx.fem.form(self.c * n_coord * ufl.dx)
                    current_mass = self.mesh.comm.allreduce(
                        dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                    )

                    # Normalize to target mass
                    target_mass = float(ic_config.get("value", 1000.0))
                    if current_mass > 0.0:
                        scaling_factor = target_mass / current_mass
                        self.c.x.array[:] *= scaling_factor
                        self.c_n.x.array[:] *= scaling_factor

                    self.C_tot_target = self.mesh.comm.allreduce(
                        dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                    )
                    print(f"  [INFO] Initial mass target (C_tot) normalized to: {self.C_tot_target:.4e}")

                except Exception as e:
                    print(f"  [WARNING] Could not calculate initial C_tot: {e}")

        # Gaussian initial condition
        elif ic_type == "gaussian":
            try:
                # Distribution parameters
                mean = ic_config.get("mean", 5.0)       # Peak position
                std_dev = ic_config.get("std_dev", 1.0) # Width
                target_mass = float(ic_config.get("amplitude", 1000.0)) # Amplitude = total mass

                # Definition and interpolation of the analytical profile (base form)
                def gaussian_expression(x):
                    return np.exp( - (x[0] - mean)**2 / (2 * std_dev**2) )
                
                self.c.interpolate(gaussian_expression)
                self.c_n.interpolate(gaussian_expression)

                # Calculation of the current mass of the interpolated profile (∫ c * n dn)
                x_coord = ufl.SpatialCoordinate(self.mesh)
                n_coord = x_coord[0]
                flux_form = dolfinx.fem.form(self.c * n_coord * ufl.dx)
                current_mass = self.mesh.comm.allreduce(
                    dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                )

                # Normalization
                if current_mass > 0.0:
                    scaling_factor = target_mass / current_mass
                    self.c.x.array[:] *= scaling_factor
                    self.c_n.x.array[:] *= scaling_factor

                    self.C_tot_target = self.mesh.comm.allreduce(
                        dolfinx.fem.assemble_scalar(flux_form), op=MPI.SUM
                    )
                    
                    # Actual peak after normalization (for info)
                    actual_peak = scaling_factor # since we had interpolated with amplitude 1.0
                    print(f"  [INFO] Applied Gaussian IC: mean={mean}, std={std_dev}")
                    print(f"  [INFO] Peak amplitude adjusted to {actual_peak:.2f} to fix total mass at {self.C_tot_target:.2e}")
                else:
                    print("  [WARNING] Initial Gaussian mass is near zero. Check mean and mesh range.")

            except Exception as e:
                print(f"  [ERROR] Error applying Gaussian IC: {e}")