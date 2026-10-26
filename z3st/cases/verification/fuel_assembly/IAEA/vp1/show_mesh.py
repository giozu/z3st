#!/usr/bin/env python3
"""
Display mesh.msh with pyvista, read through dolfinx as z3st reads it.

    python show_mesh.py              # boundary facets coloured by tag
    python show_mesh.py flats/mesh.msh   # another mesh (a variant's)
    python show_mesh.py --zmax 0.3   # only the bottom 0.3 m, to see the elements
    python show_mesh.py --volume     # hexahedra instead of tagged facets
    python show_mesh.py --save mesh.png

Tags: 10 bottom, 11 top, 12 inner, 13 outer (facets); 1 steel (cells).
"""

import argparse
import os

import numpy as np
import pyvista as pv
from dolfinx import plot
from mpi4py import MPI

try:
    from dolfinx.io import gmsh as gmshio
except ImportError:
    from dolfinx.io import gmshio

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = {10: "bottom", 11: "top", 12: "inner", 13: "outer"}


def linear(topology, celltypes, x):
    """dolfinx exports VTK Lagrange cells; at degree 1 the node order equals the
    linear VTK cell, which pyvista draws without triangulating the faces."""
    lagrange_to_linear = {pv.CellType.LAGRANGE_HEXAHEDRON: pv.CellType.HEXAHEDRON,
                          pv.CellType.LAGRANGE_QUADRILATERAL: pv.CellType.QUAD}
    types = np.array([lagrange_to_linear.get(c, c) for c in celltypes], dtype=np.uint8)
    return pv.UnstructuredGrid(topology, types, x)


def load(path=os.path.join(HERE, "mesh.msh")):
    r = gmshio.read_from_msh(path, MPI.COMM_WORLD, 0, gdim=3)
    mesh, ft = (r[0], r[2]) if isinstance(r, tuple) else (r.mesh, r.facet_tags)
    cells = linear(*plot.vtk_mesh(mesh))
    mesh.topology.create_connectivity(2, 0)
    facets = linear(*plot.vtk_mesh(mesh, 2, ft.indices))
    facets.cell_data["tag"] = ft.values
    return cells, facets


def below(grid, zmax):
    return grid.extract_cells(np.flatnonzero(grid.cell_centers().points[:, 2] <= zmax))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("msh", nargs="?", default=os.path.join(HERE, "mesh.msh"))
    ap.add_argument("--zmax", type=float, help="show only z <= zmax (m)")
    ap.add_argument("--volume", action="store_true", help="show hexahedra, not facets")
    ap.add_argument("--save", help="write a screenshot instead of opening a window")
    a = ap.parse_args()

    cells, facets = load(a.msh)
    print(f"{cells.n_cells} hexahedra, {cells.n_points} nodes")
    for tag, n in zip(*np.unique(facets["tag"], return_counts=True)):
        print(f"  facet tag {tag} {NAMES.get(tag, '?'):6s} {n}")

    grid = cells if a.volume else facets
    if a.zmax is not None:
        grid = below(grid, a.zmax)

    p = pv.Plotter(off_screen=bool(a.save))
    if a.volume:
        p.add_mesh(grid, color="lightsteelblue", show_edges=True)
    else:
        grid.cell_data["region"] = np.array([NAMES[t] for t in grid["tag"]])
        p.add_mesh(grid, scalars="region", show_edges=True, cmap="Set2")
    p.add_axes()
    p.view_isometric()
    if a.save:
        p.screenshot(a.save)
        print(f"[INFO] {a.save}")
    else:
        p.show()
