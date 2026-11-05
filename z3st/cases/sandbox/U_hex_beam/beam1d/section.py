# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Thin-walled regular hexagon: section properties and free strains.

The section is the undeformed wall Omega_0, the same at every z: the region
between the outer hexagon (circumradius a_o = D/2, D the outer
corner-to-corner diameter of geometry.yaml) and the inner one,
a_i = a_o - 2 t / sqrt(3) (inner apothem = outer apothem - t). These are the
numbers z3st meshes the 3D wall from, corners sharp. Corner k sits at
phi_0 + k pi/3, phi_0 = pi/2 for 'vertical' (corners on +-y), 0 for 'flat';
flat k runs from corner k to corner k + 1.

Each flat is split at the fractions SPLITS (mid-flat and quarter points) into
trapezoids, integrated by Gauss-Legendre along the flat (xi) and through the
wall (eta). The lines of constant xi are radial, so every quadrature point
carries its flat index and its fraction f along the flat, the coordinates a
temperature table on the contour needs (beam1d.thermal). Integrals of
polynomials of low degree are exact, and so are those of a T piecewise linear
along the contour with kinks at the split points:

    A_w = 3 sqrt(3) / 2 (a_o^2 - a_i^2) = 6 a_m t
    I_w = 5 sqrt(3) / 16 (a_o^4 - a_i^4)     (any centroidal axis, I_xy = 0)
    S_x = S_y = 0                            (the axis is the centroid)

An eigenstrain eps*(x, y, z) gives the free axial strain and the free
curvature vector (log, Eq. eigencurv; Ma et al. 2019, Eqs. 30-31):

    eps0*(z)  = 1 / A_w  int_{A_w} eps*(x, y, z) dA
    kappa*(z) = 1 / I_w  int_{A_w} eps*(x, y, z) (x, y) dA

For the hexagon the two bending planes are uncoupled; the centroid moves with
u_a'' = -kappa*_a in a free bow (beam1d.solver).
"""

import numbers

import numpy as np


class HexSection:
    """Hexagonal wall from geometry.yaml (D, t, orientation)."""

    # fractions along each flat where the wall is split: the quarter points are
    # the sector edges of the 'ma12' table model, the mid-flat its key point S
    SPLITS = (0.25, 0.5, 0.75)

    def __init__(self, D, t, orientation="flat", n_gauss=4, splits=SPLITS):
        self.t = float(t)
        self.a_o = 0.5 * float(D)
        self.a_i = self.a_o - 2.0 * self.t / np.sqrt(3.0)
        self.a_m = 0.5 * (self.a_o + self.a_i)
        self.orientation = orientation
        # first corner on +y for 'vertical' (corners on +-y), on +x for 'flat'
        self.phi0 = np.pi / 2 if orientation == "vertical" else 0.0
        ang = self.phi0 + np.arange(7) * np.pi / 3
        self.ring = np.column_stack([np.cos(ang), np.sin(ang)])   # unit corners, closed
        P_o, P_i = self.a_o * self.ring, self.a_i * self.ring

        g, wg = np.polynomial.legendre.leggauss(n_gauss)
        xi, eta = (a.ravel() for a in np.meshgrid(g, g, indexing="ij"))
        w2 = np.outer(wg, wg).ravel()
        # bilinear map of each trapezoid: inner face eta = -1, outer face eta = +1
        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
        dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])

        f_edges = np.concatenate([[0.0], np.asarray(splits, dtype=float), [1.0]])
        xy, dA, flat, frac = [], [], [], []
        for k in range(6):
            for f0, f1 in zip(f_edges[:-1], f_edges[1:]):
                pi = [P_i[k] + f * (P_i[k + 1] - P_i[k]) for f in (f0, f1)]
                po = [P_o[k] + f * (P_o[k + 1] - P_o[k]) for f in (f0, f1)]
                corners = np.array([pi[0], pi[1], po[1], po[0]])
                J_xi, J_eta = dN_dxi.T @ corners, dN_deta.T @ corners
                xy.append(N.T @ corners)
                dA.append(w2 * np.abs(J_xi[:, 0] * J_eta[:, 1] - J_xi[:, 1] * J_eta[:, 0]))
                flat.append(np.full(len(w2), k))
                frac.append(f0 + (f1 - f0) * 0.5 * (xi + 1.0))
        self.xy = np.vstack(xy)          # (n_pts, 2) quadrature points (m)
        self.dA = np.concatenate(dA)     # (n_pts,) weights (m^2)
        self.flat = np.concatenate(flat) # (n_pts,) flat index k
        self.frac = np.concatenate(frac) # (n_pts,) fraction along the flat, 0 at corner k

        self.A = float(self.dA.sum())
        self.S_x = float(np.sum(self.dA * self.xy[:, 0]))
        self.S_y = float(np.sum(self.dA * self.xy[:, 1]))
        self.I_xx = float(np.sum(self.dA * self.xy[:, 1] ** 2))
        self.I_yy = float(np.sum(self.dA * self.xy[:, 0] ** 2))
        self.I_xy = float(np.sum(self.dA * self.xy[:, 0] * self.xy[:, 1]))
        self.I = self.I_yy

    def closed_form(self):
        """A_w and I_w of the hexagonal annulus."""
        return {
            "A": 1.5 * np.sqrt(3.0) * (self.a_o**2 - self.a_i**2),
            "I": 5.0 * np.sqrt(3.0) / 16.0 * (self.a_o**4 - self.a_i**4),
        }

    def radius_of_gyration(self):
        """r = sqrt(I_w / A_w) (m)."""
        return float(np.sqrt(self.I / self.A))

    def modulus_weighted(self, E):
        """B0 = int E dA, B1 = int E (x, y) dA, B2 = int E (x, y)(x, y)^T dA (Ma et al. 2019, Eq. 15).

        Only a uniform E is accepted: then B0 = E A_w, B1 = E (S_x, S_y) = 0 and
        B2 = E I_w 1, so the axial and the two bending problems decouple
        (Ma Eq. 17). An E(T) card would shift the modulus-weighted centroid and
        couple them: not implemented.
        """
        if not isinstance(E, numbers.Real):
            raise NotImplementedError(f"E = {E!r}: an E varying over the section couples axial and bending "
                                      "(B1 != 0), not implemented in beam1d")
        E = float(E)
        return {"B0": E * self.A, "B1": E * np.array([self.S_x, self.S_y]),
                "B2": E * np.array([[self.I_yy, self.I_xy], [self.I_xy, self.I_xx]])}

    def contour_coords(self, x, y):
        """Flat index k and fraction f along the flat of the points (x, y), along the radial line.

        Exact for the quadrature points (their xi lines are radial); any point of
        the wall maps to the contour point at the same polar angle.
        """
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        rel = np.mod(np.arctan2(y, x) - self.phi0, 2.0 * np.pi)
        k = np.minimum((rel // (np.pi / 3)).astype(int), 5)
        C0, C1 = self.ring[k], self.ring[k + 1]
        e = C1 - C0
        d0, d1 = np.cos(rel + self.phi0), np.sin(rel + self.phi0)
        # C0 + f e on the ray of direction d: cross(C0 + f e, d) = 0
        f = -(C0[..., 0] * d1 - C0[..., 1] * d0) / (e[..., 0] * d1 - e[..., 1] * d0)
        return k, np.clip(f, 0.0, 1.0)

    def field(self, fun, z):
        """fun(x) at the quadrature points and the heights z, shape (n_pts, len(z)).

        ``fun`` takes x = (x, y, z) as numpy arrays, like vp1_temperature.T_vp1.
        """
        z = np.atleast_1d(np.asarray(z, dtype=float))
        x = (self.xy[:, 0][:, None], self.xy[:, 1][:, None], z[None, :])
        return np.broadcast_to(fun(x), (len(self.dA), len(z)))

    def mean(self, eps, z):
        """eps0* (-): (1/A) int eps* dA at the heights z, shape (len(z),)."""
        return self.dA @ self.field(eps, z) / self.A

    def curvature(self, eps, z):
        """kappa* (1/m), shape (len(z), 2): (1/I) int eps* (x, y) dA at the heights z."""
        return ((self.dA[:, None] * self.xy).T @ self.field(eps, z)).T / self.I

    def fibre_stress(self, E, eps_pts, duz, d2ux, d2uy):
        """sigma_zz = E (u_z' - x u_x'' - y u_y'' - eps*) at the quadrature points, shape (n_pts, n_z).

        eps_pts is eps* at the points (self.field), the other arguments are the
        beam derivatives at the same heights.
        """
        x, y = self.xy[:, 0][:, None], self.xy[:, 1][:, None]
        return E * (np.asarray(duz)[None, :] - x * np.asarray(d2ux)[None, :]
                    - y * np.asarray(d2uy)[None, :] - eps_pts)

    def resultants(self, sigma):
        """N = int sigma dA and (int sigma x dA, int sigma y dA) per height.

        int sigma x dA = -M_x with M_x = EI (u_x'' + kappa*_x) the bending
        moment of the x plane (beam1d.solver), likewise for y.
        """
        return self.dA @ sigma, (self.dA[:, None] * self.xy).T @ sigma

    def thermal_shear(self, n=2000):
        """c (m^2) of the thermal shear angle gamma = w' - psi = -c dkappa*/dz (thin wall, mid-wall contour).

        Not used by the Euler-Bernoulli beam (psi = w'): it is the 3D effect of a
        transverse T difference that varies along z, theta = g(z) h. The in-plane
        free thermal displacement per unit alpha g, u_th = ((h^2 - p^2)/2, h p),
        changes with z; with no shear flow in the wall (V = 0) the warping
        Phi' = -u_th . t along the contour closes, and its plane fit on h shifts
        the section rotation psi from the slope w' of the section mean:
            c = oint Phi h ds / oint h^2 ds .
        Circle of radius R: c = R^2 / 2. Regular hexagon, R the mid-wall
        circumradius: c = 0.410 R^2, for any orientation.
        """
        R = self.a_m
        ang = np.arange(7) * np.pi / 3
        P = R * np.column_stack([np.cos(ang), np.sin(ang)])
        s = (np.arange(n) + 0.5) / n
        X = np.vstack([a + np.outer(s, b - a) for a, b in zip(P[:-1], P[1:])])
        t = np.repeat((P[1:] - P[:-1]) / R, n, axis=0)   # side length = R
        ds = R / n
        h, p = X[:, 0], X[:, 1]
        Phi = -np.cumsum((0.5 * (h**2 - p**2) * t[:, 0] + h * p * t[:, 1]) * ds)
        Phi -= Phi.mean()
        return float(np.sum(Phi * h) / np.sum(h**2))
