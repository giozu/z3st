# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import basix
import dolfinx



class FiniteElementSetup:
    def __init__(self):
        print("__FiniteElementSetup initializer__")
        self._setup_function_space()

    def _setup_function_space(self):
        """
        Set up the function space for the problem.
        """

        # --. Separate function space --..
        mech_config = self.input_file.get("mechanical", {})
        mech_degree = mech_config.get("order", 1)
        print(f"Mechanical element order: {mech_degree}")

        # --. Temperature --..
        self.V_t = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))
        print("Thermal function space (V_t):", self.V_t)
        
        # --. Displacement --..
        self.V_m = dolfinx.fem.functionspace(self.mesh, ("Lagrange", mech_degree, (self.mesh.topology.dim,)))
        print("Mechanical function space (V_m):", self.V_m)
        
        # --. Damage --..
        if self.on.get("damage", False):
            self.V_d = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))
            print("Scalar function space (V_d):", self.V_d)
        
        # --. Scalar field --..
        self.Q = dolfinx.fem.functionspace(self.mesh, ("DG", 0))
        print("Scalar function space (Q):", self.Q)

        # --. Cluster dynamics --..
        if self.on.get("cluster", False):
            self.V_c = dolfinx.fem.functionspace(self.mesh, ("DG", 1))
            print("Cluster function space (V_c):", self.V_c)
    
        # --. Plasticity --..
        if self.on.get("plasticity", False):
            
            self.q_degree = 2 * mech_degree + 1
            print(f"Plasticity function spaces (V_pl_tensor, Q_pl) initializing with Quadrature degree {self.q_degree}")
            
            # Tensor
            el_tensor = basix.ufl.quadrature_element(self.mesh.topology.cell_name(), value_shape=(3, 3), degree=self.q_degree)
            self.V_pl = dolfinx.fem.functionspace(self.mesh, el_tensor)
            
            # Scalar
            el_scalar = basix.ufl.quadrature_element(self.mesh.topology.cell_name(), value_shape=(), degree=self.q_degree)
            self.Q_pl = dolfinx.fem.functionspace(self.mesh, el_scalar)
            print("Plasticity function space (V_pl_tensor):", self.V_pl)
            print("Plasticity function space (Q_pl):", self.Q_pl)
        else:
            self.q_degree = None

