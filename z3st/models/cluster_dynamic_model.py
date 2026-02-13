"""
Cluster Dynamics Model for Z3ST

This module implements a 1D cluster dynamics model for simulating
the evolution of defect cluster size distributions.
"""

import dolfinx
import ufl
import numpy as np
from petsc4py import PETSc


class ClusterDynamicsModel:
    """
    Cluster dynamics model for defect evolution.
    
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
        
        # Legacy parameter (kept for compatibility)
        self.D = self.D_cluster
        
        print("[ClusterDynamicsModel] Initialized")
        print(f"  → Advection velocity v = {self.v_cluster}")
        print(f"  → Diffusion coefficient D = {self.D_cluster}")
