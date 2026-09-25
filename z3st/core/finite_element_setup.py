# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
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
        # The cohesive model carries its phase field in the same space and the
        # same Function (self.D), so V_d is needed for either route.
        if self.on.get("damage", False) or self.on.get("cohesive", False):
            self.V_d = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))
            print("Scalar function space (V_d):", self.V_d)

        # --. Cohesive fracture: mixed (u, eigenstrain) space --..
        # The eigenstrain is piecewise constant, matching the constant strain
        # within a linear element; in the multi-axial case it needs only its
        # trace and its deviatoric norm (Vicentini et al. 2026, Sec. 5.1).
        if self.on.get("cohesive", False):
            cell = self.mesh.basix_cell()
            u_el = basix.ufl.element(
                "Lagrange", cell, mech_degree, shape=(self.mesh.topology.dim,)
            )
            eta_el = basix.ufl.element("DG", cell, 0)
            n_eta = 1 if self.mesh.topology.dim == 1 else 2
            self.W = dolfinx.fem.functionspace(
                self.mesh, basix.ufl.mixed_element([u_el] + [eta_el] * n_eta)
            )
            print(f"Cohesive mixed function space (W): {self.W} [1 + {n_eta} blocks]")


        # --. Scalar field --..
        self.Q = dolfinx.fem.functionspace(self.mesh, ("DG", 0))
        print("Scalar function space (Q):", self.Q)

        # --. Cluster dynamics --..
        if self.on.get("cluster", False):
            self.V_c = dolfinx.fem.functionspace(self.mesh, ("DG", 1))
            print("Cluster function space (V_c):", self.V_c)
            
        # --. Porosity migration --..
        # Two discretisations are available, selected by porosity.discretisation:
        #  - "cg"  (default): Lagrange-1, used by the SU/SUPG-stabilised solve
        #    (Barani et al. 2022). Benchmark-validated path.
        #  - "dg": discontinuous Lagrange-1, used by the upwind(+SIPG) solve —
        #    the same operator family as cluster dynamics (V_c). The pore
        #    velocity v = mobility(T) grad(T) is advection-dominated, so the
        #    upwind facet flux supplies the stabilisation directly.
        if self.on.get("porosity", False):
            self.porosity_discretisation = str(
                self.input_file.get("porosity", {}).get("discretisation", "cg")
            ).lower()
            family = "DG" if self.porosity_discretisation == "dg" else "Lagrange"
            self.V_p = dolfinx.fem.functionspace(self.mesh, (family, 1))
            print(f"Porosity function space (V_p): {self.V_p} [{self.porosity_discretisation}]")
    
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

