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
        self.V_t = dolfinx.fem.functionspace(
            self.mesh, ("Lagrange", 1)
        )  # scalar function space (temperature)
        self.V_m = dolfinx.fem.functionspace(
            self.mesh, ("Lagrange", 1, (self.mesh.topology.dim,))
        )  # vector function space (displacement)
        self.Q = dolfinx.fem.functionspace(self.mesh, ("DG", 0))

        # print(f"DoFs number for V_m: {self.V_m.dofmap.index_map.size_global}")
        # print(f"DoFs number for V_t: {self.V_t.dofmap.index_map.size_global}")

        print("Mechanical function space (V_m):", self.V_m)
        print("Thermal function space (V_t):", self.V_t)
        print("Scalar function space (Q):", self.Q)

        # --. Mixed function spaces --..
        # P1_3 = basix.ufl.element("Lagrange", self.mesh.basix_cell(), 1, shape=(self.tdim,), dtype=dolfinx.default_real_type)
        # P1_1 = basix.ufl.element("Lagrange", self.mesh.basix_cell(), 1, dtype=dolfinx.default_real_type)
        # Same P1_3 and P1_1 of above, but for coherence I prefere to define the following from V_m and V_t
        P1_3 = self.V_m.ufl_element()
        P1_1 = self.V_t.ufl_element()

        self.W = dolfinx.fem.functionspace(self.mesh, basix.ufl.mixed_element([P1_3, P1_1]))

        # The following line is currently not used, self.W is used instead
        # self.V_m_mixed, self.V_t_mixed = dolfinx.fem.functionspace(self.mesh, P1_3), dolfinx.fem.functionspace(self.mesh, P1_1)
        # print("Mechanical function space (V_m_mixed):", self.V_m_mixed)
        # print("Thermal function space (V_t_mixed):", self.V_t_mixed)
