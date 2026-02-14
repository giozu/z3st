# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---


import dolfinx
import ufl
import numpy as np
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
        - v = 1.0 (advection velocity)
        - D = 0.5 (diffusion coefficient)
        - c(n,t) is cluster density at size n and time t
        """
        # Physical parameters for advection-diffusion equation
        self.v_cluster = 1.0  # Advection velocity (cluster growth rate)
        self.D_cluster = 0.5  # Diffusion coefficient (cluster size fluctuations)
        self.C_tot_target = None

        print("[ClusterDynamicsModel] Initialized")

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
            else:
                print(f"  [ERROR] No value specified for constant IC.")

        # Gaussian initial condition
        elif ic_type == "gaussian":
            try:
                # Distribution parameters: mean size, standard deviation, peak amplitude
                mean = ic_config.get("mean", 50.0)
                std_dev = ic_config.get("std_dev", 5.0)
                amplitude = ic_config.get("amplitude", 1000.0)

                # Analytical profile: c(n) = A * exp( - (n - μ)² / (2σ²) )
                def gaussian_expression(x):
                    return amplitude * np.exp( - (x[0] - mean)**2 / (2 * std_dev**2) )
                
                # Interpolate the expression into the finite element function
                self.c.interpolate(gaussian_expression)
                self.c_n.interpolate(gaussian_expression)
                print(f"  [INFO] Applied Gaussian IC: mean={mean}, std={std_dev}, amp={amplitude}")
            except Exception as e:
                print(f"  [ERROR] Error applying Gaussian IC: {e}")

        # Mass conservation target initialization
        # Calculate C_tot = ∫ c(n,0) * n dn to be conserved during the simulation
        try:
            x_coord = ufl.SpatialCoordinate(self.mesh)
            flux_form = dolfinx.fem.form(self.c * x_coord[0] * ufl.dx)
            self.C_tot_target = dolfinx.fem.assemble_scalar(flux_form)
            print(f"  [INFO] Initial mass target (C_tot): {self.C_tot_target:.4e}")
        except Exception as e:
            print(f"  [WARNING] Could not calculate initial C_tot: {e}")