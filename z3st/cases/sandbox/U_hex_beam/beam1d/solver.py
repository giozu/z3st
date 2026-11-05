# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Mesh, assembly, boundary conditions, solution and evaluation (Oechsner
Sects. 5.2.3 and 5.3.3).

One bending plane a = x or y: unknown u_a(z), the centroid displacement,
global dofs (u_i, phi_i) of node i at 2i, 2i + 1, phi = u_a'. The fibre strain
is eps = u_z' - x u_x'' - y u_y'', so a free bow has u_a'' = -kappa*_a: the
plane is solved with the free curvature k_a = -kappa*_a (beam1d.element) and
its bending moment is M_a = EI (u_a'' - k_a). Axial part: unknown u_z, dof i
at node i (linear rod).

Pads are unilateral contacts on one dof: the restraint sits at side * gap and
forbids side * u > gap. They are solved by an active set (Ma et al. 2019,
Fig. 13, for any number of pads): an active pad is a prescribed displacement
u = side * gap, its reaction r = (K U - f)[dof] the force the restraint puts
on the beam; a pad that penetrates is activated, one that would pull
(side * r > 0) is released, until nothing changes.
"""

import numpy as np

from beam1d import element


def mesh(Lz, nz, z_nodes=()):
    """nz uniform elements on [0, Lz], plus the forced nodes (the stations)."""
    z = np.concatenate([np.linspace(0.0, Lz, nz + 1), np.asarray(z_nodes, dtype=float)])
    z = np.unique(np.round(z, 12))
    return z[(z >= 0.0) & (z <= Lz)]


def node_at(z, zq, tol=1e-9):
    """Index of the node at height zq (a station is always a node)."""
    i = int(np.argmin(np.abs(z - zq)))
    if abs(z[i] - zq) > tol:
        raise ValueError(f"no node at z = {zq} m (nearest {z[i]} m)")
    return i


def gauss_points(z, n_gauss):
    """Heights (n_el, n_gauss) of the element Gauss points, with xi and weights on [0, 1]."""
    xi, w = element.gauss(n_gauss)
    L = np.diff(z)
    return z[:-1, None] + xi[None, :] * L[:, None], xi, w


def assemble(z, EI, k_g, xi, w, q_g=0.0):
    """Global K and f of a bending plane for the free curvature k_g (n_el, n_gauss) and the load q_g (N/m)."""
    n = 2 * len(z)
    K, f = np.zeros((n, n)), np.zeros(n)
    q_g = np.broadcast_to(q_g, np.shape(k_g))
    for e, L in enumerate(np.diff(z)):
        d = slice(2 * e, 2 * e + 4)
        K[d, d] += element.stiffness(EI, L)
        f[d] += element.load(EI, L, k_g[e], xi, w) + element.load_q(L, q_g[e], xi, w)
    return K, f


def assemble_rod(z, EA, e0_g, xi, w, q_g=0.0):
    """Global K and f of the axial rod for the free strain e0_g (n_el, n_gauss) and the load q_g (N/m)."""
    n = len(z)
    K, f = np.zeros((n, n)), np.zeros(n)
    q_g = np.broadcast_to(q_g, np.shape(e0_g))
    for e, L in enumerate(np.diff(z)):
        d = slice(e, e + 2)
        K[d, d] += element.rod_stiffness(EA, L)
        f[d] += element.rod_load(EA, L, e0_g[e], xi, w, q_g[e])
    return K, f


def solve(K, f, fixed, values=None):
    """U with U[fixed] = values (default 0), and the reactions K U - f at the fixed dofs."""
    fixed = np.asarray(fixed, dtype=int)
    free = np.setdiff1d(np.arange(len(f)), fixed)
    U = np.zeros(len(f))
    if values is not None:
        U[fixed] = values
    U[free] = np.linalg.solve(K[np.ix_(free, free)], f[free] - K[np.ix_(free, fixed)] @ U[fixed])
    return U, (K @ U - f)[fixed]


def solve_contact(K, f, fixed, values, contacts, max_iters=20):
    """Active-set solution with unilateral pads.

    contacts: list of dicts with ``dof``, ``gap`` (m) and ``side`` (+1 or -1).
    Returns U, the reactions at ``fixed``, the pad forces (0 when open), the
    active flags and the number of iterations.
    """
    fixed, values = list(fixed), list(values)
    active = np.zeros(len(contacts), dtype=bool)
    for it in range(1, max_iters + 1):
        on = [c for c, a in zip(contacts, active) if a]
        U, r = solve(K, f, fixed + [c["dof"] for c in on], values + [c["side"] * c["gap"] for c in on])
        R, r_pad = r[:len(fixed)], iter(r[len(fixed):])
        force = np.array([next(r_pad) if a else 0.0 for a in active])
        tol_u = 1e-12 * max([1.0] + [abs(c["gap"]) for c in contacts])
        tol_f = 1e-12 * max(1.0, float(np.max(np.abs(r))) if len(r) else 1.0)
        changed = False
        for j, c in enumerate(contacts):
            if active[j] and c["side"] * force[j] > tol_f:
                active[j], changed = False, True        # the restraint would pull: release
            elif not active[j] and c["side"] * U[c["dof"]] - c["gap"] > tol_u:
                active[j], changed = True, True         # the pad penetrates: close it
        if not changed:
            return U, R, force, active, it
    raise RuntimeError(f"contact active set not converged in {max_iters} iterations (active {active.tolist()})")


def evaluate(z, U, zq):
    """u, u', u'' at the heights zq from the nodal (u, phi) by Hermite interpolation."""
    zq = np.atleast_1d(np.asarray(zq, dtype=float))
    e = np.clip(np.searchsorted(z, zq, side="right") - 1, 0, len(z) - 2)
    L = z[e + 1] - z[e]
    N, dN, d2N = element.shape((zq - z[e]) / L, L)
    Ue = np.stack([U[2 * e], U[2 * e + 1], U[2 * e + 2], U[2 * e + 3]])
    return (N * Ue).sum(0), (dN * Ue).sum(0), (d2N * Ue).sum(0)


def evaluate_rod(z, U, zq):
    """u_z, u_z' at the heights zq from the nodal u_z (linear rod; u_z' is constant per element)."""
    zq = np.atleast_1d(np.asarray(zq, dtype=float))
    e = np.clip(np.searchsorted(z, zq, side="right") - 1, 0, len(z) - 2)
    L = z[e + 1] - z[e]
    N, dN = element.rod_shape((zq - z[e]) / L, L)
    Ue = np.stack([U[e], U[e + 1]])
    return (N * Ue).sum(0), (dN * Ue).sum(0)
