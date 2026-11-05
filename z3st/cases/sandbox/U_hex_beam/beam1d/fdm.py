# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
The finite-difference scheme of Ma et al. (2019), Sect. 3.2, for a cantilever
in free bowing, as a reference for the Hermite elements (beam1d.solver).
Uniform grid x_0 = 0, ..., x_N = L, h = L / N.

Bending, Eq. (24) w'' = -M_T / EI = k, central differences, Eq. (27):

    w_{i-1} - 2 w_i + w_{i+1} = h^2 k_i ,     i = 1, ..., N - 1

The clamp needs w_0 = 0 and w'(0) = 0, but the scheme has no rotation
unknown. Eq. (28) writes w'(x_0) with the one-sided second-order formula, the
slope at x_0 of the parabola through the first three nodes,

    w'(x_0) = (-w_2 + 4 w_1 - 3 w_0) / (2 h) = 0   ->   4 w_1 = w_2  (Eq. 29),

and Eq. (27) at i = 1 then gives w_1 = h^2 k_1 / 2. From there the scheme
marches up: w_{i+1} = 2 w_i - w_{i-1} + h^2 k_i. Its error is O(h^2); it is
exact for a cubic w (k linear) away from kinks of k.

Axial, Eq. (25), first-order backward differences:

    u_i = u_{i-1} + h (eps0*(x_i) + N(x_i) / EA),     u_0 = 0,

with N the axial force (0 in free bowing without gravity). Error O(h).
"""

import numpy as np


def _grid(x):
    x = np.asarray(x, dtype=float)
    h = x[1] - x[0]
    if not np.allclose(np.diff(x), h, rtol=1e-10, atol=0.0):
        raise ValueError("the scheme of Ma Sect. 3.2 needs a uniform grid")
    return x, h


def bending(x, k):
    """w at the nodes x for w'' = k (k at the nodes), clamped at x_0 by Eqs. (28)-(29)."""
    x, h = _grid(x)
    k = np.broadcast_to(np.asarray(k, dtype=float), x.shape)
    w = np.zeros_like(x)
    w[1] = 0.5 * h**2 * k[1]
    for i in range(1, len(x) - 1):
        w[i + 1] = 2.0 * w[i] - w[i - 1] + h**2 * k[i]
    return w


def axial(x, e0, N_over_EA=0.0):
    """u at the nodes x from eps0* (and N / EA) at the nodes, u_0 = 0, Eq. (25)."""
    x, h = _grid(x)
    rhs = np.broadcast_to(np.asarray(e0, dtype=float) + N_over_EA, x.shape)
    u = np.zeros_like(x)
    u[1:] = np.cumsum(h * rhs[1:])
    return u
