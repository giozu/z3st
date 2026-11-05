# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Two-node elements: Euler-Bernoulli beam (cubic Hermite) for bending, linear
rod for the axial part.

Reference: Oechsner, "Elasto-Plasticity of Frame Structure Elements"
(Springer, 2014). Rod, Sect. 5.2.1.1: weak form Eq. (5.62), stiffness
Eqs. (5.73) and (5.79). Beam, Sect. 5.3: weak form Eq. (5.213), shape
functions Eqs. (5.222)-(5.225), element stiffness Eq. (5.236), equivalent
nodal loads Eq. (5.255).

Beam. Nodal unknowns (u_1, phi_1, u_2, phi_2), phi = du/dz. With x = xi L,
xi in [0, 1]:

    N_1u   = 1 - 3 xi^2 + 2 xi^3        N_2u   = 3 xi^2 - 2 xi^3
    N_1phi = L (xi - 2 xi^2 + xi^3)     N_2phi = L (-xi^2 + xi^3)

Free curvature. With the fibre strain eps = -eta u'' (Oechsner Eq. 4.15) and
sigma = E (eps - eps*), the beam bends freely with u'' = k, k the free
curvature of the plane (beam1d.section). The energy
1/2 int EI (u'' - k)^2 dz - int q u dz gives K_e u_e = f_e with

    K_e = int EI N'' N''^T dz,     f_e = int EI k N'' dz + int q N dz .

For k constant, f_e = EI k (0, -1, 0, +1): two opposite nodal moments; for q
constant, int q N dz = q L (1/2, L/12, 1/2, -L/12), Eq. (5.255).

Rod. Nodal unknowns (u_1, u_2), N_1 = 1 - xi, N_2 = xi. The energy
1/2 int EA (u' - eps0*)^2 dz - int q_z u dz gives

    K_e = EA / L [[1, -1], [-1, 1]],     f_e = EA int eps0* N' dz + int q_z N dz .
"""

import numpy as np


def gauss(n):
    """Gauss-Legendre points and weights on [0, 1]."""
    g, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (g + 1.0), 0.5 * w


def shape(xi, L):
    """N, dN/dz, d2N/dz2 at xi in [0, 1], each of shape (4, len(xi))."""
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    N = np.array([1 - 3 * xi**2 + 2 * xi**3,
                  L * (xi - 2 * xi**2 + xi**3),
                  3 * xi**2 - 2 * xi**3,
                  L * (-xi**2 + xi**3)])
    dN = np.array([(-6 * xi + 6 * xi**2) / L,
                   1 - 4 * xi + 3 * xi**2,
                   (6 * xi - 6 * xi**2) / L,
                   -2 * xi + 3 * xi**2])
    d2N = np.array([(-6 + 12 * xi) / L**2,
                    (-4 + 6 * xi) / L,
                    (6 - 12 * xi) / L**2,
                    (-2 + 6 * xi) / L])
    return N, dN, d2N


def stiffness(EI, L, n_gauss=2):
    """K_e = int_0^L EI N'' N''^T dz; N'' is linear, so 2 points are exact."""
    xi, w = gauss(n_gauss)
    B = shape(xi, L)[2]
    return EI * L * (B * w) @ B.T


def stiffness_closed(EI, L):
    """K_e of Oechsner Eq. (5.236)."""
    return EI / L**3 * np.array([[12.0, 6 * L, -12.0, 6 * L],
                                 [6 * L, 4 * L**2, -6 * L, 2 * L**2],
                                 [-12.0, -6 * L, 12.0, -6 * L],
                                 [6 * L, 2 * L**2, -6 * L, 4 * L**2]])


def load(EI, L, k, xi, w):
    """f_e = int_0^L EI k N'' dz, with k the free curvature at the points xi (weights w)."""
    B = shape(xi, L)[2]
    return EI * L * B @ (w * k)


def load_q(L, q, xi, w):
    """f_e = int_0^L q N dz, q the transverse load per unit length (N/m) at the points xi."""
    N = shape(xi, L)[0]
    return L * N @ (w * np.broadcast_to(q, np.shape(xi)))


def rod_shape(xi, L):
    """N, dN/dz of the linear rod at xi in [0, 1], each of shape (2, len(xi))."""
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    return np.array([1 - xi, xi]), np.array([-np.ones_like(xi), np.ones_like(xi)]) / L


def rod_stiffness(EA, L):
    """K_e of Oechsner Eq. (5.79)."""
    return EA / L * np.array([[1.0, -1.0], [-1.0, 1.0]])


def rod_load(EA, L, e0, xi, w, q=0.0):
    """f_e = int_0^L (EA eps0* N' + q_z N) dz, eps0* and q_z at the points xi."""
    N, dN = rod_shape(xi, L)
    return L * (EA * dN @ (w * e0) + N @ (w * np.broadcast_to(q, np.shape(xi))))
