# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Section properties of a z3st 3D mesh, to check that the beam section is the
undeformed cross-section of the 3D model: A, S_x, S_y, I_xx, I_yy, I_xy
integrated over the facets of one physical surface (e.g. 'bottom'), read
with the gmsh Python API. Linear triangles and quadrilaterals, 3 x 3 Gauss
points per facet (exact for the polynomials of degree 2 on straight-sided
facets).
"""

import numpy as np


def section_from_msh(path, surface="bottom"):
    """dict A, S_x, S_y, I_xx = int y^2, I_yy = int x^2, I_xy of the physical surface of a .msh file."""
    import gmsh

    started = not gmsh.isInitialized()
    if started:
        gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    try:
        gmsh.open(path)
        tags = [t for d, t in gmsh.model.getPhysicalGroups(2) if gmsh.model.getPhysicalName(2, t) == surface]
        if not tags:
            raise KeyError(f"no physical surface '{surface}' in {path}")
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        X = dict(zip(node_tags, coords.reshape(-1, 3)))
        g, wg = np.polynomial.legendre.leggauss(3)
        acc = np.zeros(6)
        for ent in gmsh.model.getEntitiesForPhysicalGroup(2, tags[0]):
            for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2, ent)):
                name, _, _, n_nodes, _, _ = gmsh.model.mesh.getElementProperties(etype)
                for nodes in enodes.reshape(-1, n_nodes):
                    P = np.array([X[n][:2] for n in nodes])
                    if n_nodes == 4:      # bilinear quadrilateral on [-1, 1]^2
                        xi, eta = (a.ravel() for a in np.meshgrid(g, g, indexing="ij"))
                        w = np.outer(wg, wg).ravel()
                        N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                                             (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
                        dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
                        deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
                        Jx, Je = dxi.T @ P, deta.T @ P
                        dA = w * np.abs(Jx[:, 0] * Je[:, 1] - Jx[:, 1] * Je[:, 0])
                    elif n_nodes == 3:    # linear triangle, Gauss on the collapsed square
                        u, v = (0.5 * (a.ravel() + 1) for a in np.meshgrid(g, g, indexing="ij"))
                        w = 0.25 * np.outer(wg, wg).ravel() * (1 - u)
                        r, s = u, v * (1 - u)
                        N = np.array([1 - r - s, r, s])
                        dA = w * abs(np.cross(P[1] - P[0], P[2] - P[0]))
                    else:
                        raise NotImplementedError(f"facet type {name} on '{surface}'")
                    x, y = (N.T @ P).T
                    acc += [dA.sum(), dA @ x, dA @ y, dA @ y**2, dA @ x**2, dA @ (x * y)]
        return dict(zip(("A", "S_x", "S_y", "I_xx", "I_yy", "I_xy"), map(float, acc)))
    finally:
        if started:
            gmsh.finalize()
