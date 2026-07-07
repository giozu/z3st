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

    Because the mean gap is affine in the applied pressure, the pure
    fixed-point update has loop gain k_pen * C (C = radial interface
    compliance), typically 20-50, and converges only under heavy
    under-relaxation. Once two (pressure, gap) samples are available the
    update therefore switches to a SECANT step: the compliance is estimated
    from the two samples and the pressure jumps directly to the value
    consistent with the penalty equilibrium  g(p) = -p/k_pen. Guards fall
    back to the explicit update whenever the secant slope is unusable
    (no history, near-identical samples, non-physical slope, opening
    contact).

    The pressure is applied as a normal traction t = -p * n on BOTH facing
    surfaces (each with its own outward facet normal), so the bodies are
    pushed apart on penetration. Because p is evaluated from the previous
    iterate it enters the linear mechanical solve on the RHS (an updated
    Neumann load); the staggered relaxation loop drives the fixed point.

    This is the explicit counterpart of constraint-based (Lagrange) contact and
    reuses the existing staggered solver; it is a demo PCMI engine, not a
    production contact algorithm.

    Mean radial displacement on a surface is computed as a boundary-integral
    average  <u_r>_Γ = (∫_Γ u·n ds) / (∫_Γ ds)  via assemble_scalar, which is
    unambiguous under blocked vector spaces and MPI-parallel. Projecting on
    the outward facet normal (with the sign flipped on surface b, whose
    outward normal points radially inward) makes the measure valid both for
    axisymmetric r-z meshes (straight facets, n = ±e_x) and for in-plane
    r-theta disk meshes (curved arcs, n = ±e_r).
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
        # (pressure, gap) sample from the previous iteration for the secant
        # update; reset at every step boundary via reset_contact_history().
        self._sample_prev = None

        print(f"  surfaces      : '{self.contact_surface_a}' (id {self.id_a}) <-> '{self.contact_surface_b}' (id {self.id_b})")
        print(f"  initial gap g0: {self.g0 * 1e6:.2f} um")
        print(f"  penalty k_pen : {self.k_pen:.2e} Pa/m")

    def _assemble(self, form):
        """MPI-summed scalar assembly."""
        local = dolfinx.fem.assemble_scalar(form)
        return self.mesh.comm.allreduce(local, op=MPI.SUM)

    def reset_contact_history(self):
        """Drop the secant history. Call at a step boundary (or rollback):
        the free-expansion part of the gap changes with the step's thermal
        state, so a stale sample pair would produce a bogus slope."""
        self._sample_prev = None

    def update_contact_pressure(self, u):
        """
        Measure the current gap from displacement iterate ``u`` and set the
        penalty contact pressure. Returns (current_gap, pressure).
        """
        # Mean radial displacement <u·n> on each facing surface. On surface a
        # (inner body, outer face) the outward normal is +e_r; on surface b
        # (outer body, inner face) it is -e_r, hence the sign flip.
        n_vec = self._normal_as(u)
        ur_a = self._assemble(dolfinx.fem.form(ufl.dot(u, n_vec) * self._ds_a)) / self._area_a
        ur_b = -self._assemble(dolfinx.fem.form(ufl.dot(u, n_vec) * self._ds_b)) / self._area_b

        current_gap = self.g0 + ur_b - ur_a

        # ``current_gap`` is the response to the pressure applied at the
        # PREVIOUS update (still held by the Constant) — that pair is the
        # newest secant sample.
        p_applied = float(self.contact_pressure.value)
        sample_now = (p_applied, current_gap)

        p, method = self._next_pressure(sample_now)

        self._sample_prev = sample_now
        self.contact_pressure.value = PETSc.ScalarType(p)
        self._last_gap = current_gap
        self._last_pressure = p

        print(
            f"  [contact] u_r(fuel)={ur_a*1e6:+.2f} um, u_r(clad)={ur_b*1e6:+.2f} um, "
            f"gap={current_gap*1e6:+.2f} um, p={p/1e6:.3f} MPa "
            f"({'CLOSED' if p > 0 else 'open'}, {method})"
        )
        return current_gap, p

    def _next_pressure(self, sample_now):
        """Next contact pressure from the newest (pressure, gap) sample.

        Secant: with two samples of the affine map g(p) = g_free + C*p, the
        estimated compliance C gives the pressure consistent with the penalty
        equilibrium g(p*) = -p*/k_pen in one jump:

            p* = (C*p_now - g_now) / (C + 1/k_pen)

        Falls back to the explicit update p = k_pen*max(0, -g) when the slope
        is not usable. Returns (pressure, method_tag).
        """
        p_now, g_now = sample_now
        p_explicit = self.k_pen * max(0.0, -g_now)

        if self._sample_prev is None:
            return p_explicit, "explicit"

        p_old, g_old = self._sample_prev
        dp = p_now - p_old
        # Near-identical pressures give no slope information; a non-positive
        # slope is non-physical (more pressure must open the gap) and means
        # the samples straddle a thermal update or the open/closed switch.
        if abs(dp) < max(1.0, 1.0e-6 * abs(p_now)):
            return p_explicit, "explicit"
        C_est = (g_now - g_old) / dp
        if C_est <= 0.0:
            return p_explicit, "explicit"

        p_star = (C_est * p_now - g_now) / (C_est + 1.0 / self.k_pen)
        if p_star <= 0.0:
            return 0.0, "secant-open"
        return p_star, "secant"

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
        n_vec = self._normal_as(v)
        term = 0
        for ds_c in (self._ds_a, self._ds_b):
            term += w * ufl.dot(-p * n_vec, v) * ds_c
        return term

    def _normal_as(self, u):
        """Outward facet normal truncated to the vector size of ``u`` — the
        displacement space has tdim components while the gmsh-read mesh keeps
        gdim = 3, so the normal may carry a spurious trailing component."""
        n = self.normal
        return ufl.as_vector([n[i] for i in range(u.ufl_shape[0])])
