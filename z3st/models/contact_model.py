# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import ufl
from mpi4py import MPI
from petsc4py import PETSc


class ContactModel:
    """
    Explicit (fixed-point) penalty contact between two facing surfaces
    separated by an initial radial clearance (a pellet-clad gap).

    Mechanism
    --.--.--.
    Each staggered mechanical iteration the current gap is measured from the
    displacement iterate, and a contact pressure is set:

        g(u)  = g0 + <u_r>_clad - <u_r>_fuel        (current mean radial gap)
        p     = k_pen * max(0, -g)                  (penalty, compressive)

    The pressure is applied as a normal traction t = -p * n on BOTH facing
    surfaces (each with its own outward facet normal), so the bodies are
    pushed apart on penetration. Because p is evaluated from the previous
    iterate it enters the linear mechanical solve on the RHS (an updated
    Neumann load); the staggered relaxation loop drives the fixed point.

    This is the explicit counterpart of constraint-based (Lagrange) contact:
    cheap, robust, and idiomatic to the existing staggered solver — intended
    as the conference-demo PCMI engine, not a production contact algorithm.

    Mean radial displacement on a surface is computed as a boundary-integral
    average  <u_r>_Γ = (∫_Γ u_r ds) / (∫_Γ ds)  via assemble_scalar, which is
    unambiguous under blocked vector spaces and MPI-parallel.
    """

    def __init__(self):
        print("[ContactModel] initializer")

        cfg = self.input_file.get("models", {}).get("contact", {})

        # Facing surfaces: 'a' is the inner body's outer face (fuel/pellet),
        # 'b' is the outer body's inner face (cladding).
        self.contact_surface_a = cfg.get("surface_a", "lateral_1")
        self.contact_surface_b = cfg.get("surface_b", "inner_2")
        self.k_pen = float(cfg.get("penalty_stiffness", 5.0e13))  # (Pa/m)

        # Initial radial clearance g0: explicit override, else from geometry.
        g0_cfg = cfg.get("initial_gap", None)
        if g0_cfg is not None:
            self.g0 = float(g0_cfg)
        else:
            r_a = float(self.geometry.get("outer_radius_1"))
            r_b = float(self.geometry.get("inner_radius_2"))
            self.g0 = r_b - r_a

        self.id_a = self.label_map[self.contact_surface_a]
        self.id_b = self.label_map[self.contact_surface_b]

        # Contact pressure as a runtime-updated Constant (enters the weak form).
        self.contact_pressure = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(0.0))

        # Boundary measures and (constant) surface areas for the mean.
        ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.facet_tags)
        self._ds_a = ds(self.id_a)
        self._ds_b = ds(self.id_b)
        self._area_a = self._assemble(dolfinx.fem.form(dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(1.0)) * self._ds_a))
        self._area_b = self._assemble(dolfinx.fem.form(dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(1.0)) * self._ds_b))

        self._last_gap = self.g0
        self._last_pressure = 0.0

        print(f"  surfaces      : '{self.contact_surface_a}' (id {self.id_a}) <-> '{self.contact_surface_b}' (id {self.id_b})")
        print(f"  initial gap g0: {self.g0 * 1e6:.2f} um")
        print(f"  penalty k_pen : {self.k_pen:.2e} Pa/m")

    def _assemble(self, form):
        """MPI-summed scalar assembly."""
        local = dolfinx.fem.assemble_scalar(form)
        return self.mesh.comm.allreduce(local, op=MPI.SUM)

    def update_contact_pressure(self, u):
        """
        Measure the current gap from displacement iterate ``u`` and set the
        penalty contact pressure. Returns (current_gap, pressure).
        """
        # Mean radial displacement (u_r = u[0]) on each facing surface.
        ur_a = self._assemble(dolfinx.fem.form(u[0] * self._ds_a)) / self._area_a
        ur_b = self._assemble(dolfinx.fem.form(u[0] * self._ds_b)) / self._area_b

        current_gap = self.g0 + ur_b - ur_a
        penetration = max(0.0, -current_gap)
        p = self.k_pen * penetration

        self.contact_pressure.value = PETSc.ScalarType(p)
        self._last_gap = current_gap
        self._last_pressure = p

        print(
            f"  [contact] u_r(fuel)={ur_a*1e6:+.2f} um, u_r(clad)={ur_b*1e6:+.2f} um, "
            f"gap={current_gap*1e6:+.2f} um, p={p/1e6:.3f} MPa "
            f"({'CLOSED' if penetration > 0 else 'open'})"
        )
        return current_gap, p

    def contact_traction(self, v):
        """
        Contact contribution to the mechanical weak form, summed over both
        facing surfaces. Traction is t = -p * n with n the OUTWARD facet
        normal of each surface, so penetration pushes the bodies apart.

        Returns the UFL form  sum_Γ  w * dot(-p*n, v) ds  (the external-load
        term to be ADDED to L_m in the linear solve).
        """
        w = self.weight
        p = self.contact_pressure
        n = self.normal  # ufl.FacetNormal, outward on exterior facets
        n_vec = ufl.as_vector([n[0], n[1]])
        term = 0
        for ds_c in (self._ds_a, self._ds_b):
            term += w * ufl.dot(-p * n_vec, v) * ds_c
        return term
