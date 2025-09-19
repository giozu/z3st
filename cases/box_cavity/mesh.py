import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile
import pyvista
import dolfinx

# Facet markers
INNER_SPHERE = 1
XMIN, XMAX = 11, 12
YMIN, YMAX = 13, 14
ZMIN, ZMAX = 15, 16

# ID Maps
FACET_LABELS = {
    INNER_SPHERE: "inner_sphere",
    XMIN: "xmin",
    XMAX: "xmax",
    YMIN: "ymin",
    YMAX: "ymax",
    ZMIN: "zmin",
    ZMAX: "zmax",
}

def build_mesh(L=1.0, a=0.1, h_box=0.05, h_sphere=0.01,
                                    order=1, comm=MPI.COMM_WORLD,
                                    save=False, basename="box_bubble"):
    """
    Box with edge (length L) and spherical cavity (radius a) at the center.
    Mesh sizes: h_box (box), h_sphere (sphere).
    Order: polynomial order of the mesh elements.
    save=True saves mesh and markers in XDMF (for Paraview)
    """
    gmsh.initialize()
    gmsh.model.add("box_bubble")

    half = L / 2

    # Box + Sphere
    box = gmsh.model.occ.addBox(-half, -half, -half, L, L, L)
    sph = gmsh.model.occ.addSphere(0, 0, 0, a)

    # subtract sphere from box
    vol, _ = gmsh.model.occ.cut([(3, box)], [(3, sph)], removeObject=True, removeTool=False)
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(3, [vol[0][1]], tag=100, name="solid")

    # surfaces
    surfs = gmsh.model.getEntities(dim=2)

    # Mesh size
    gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), h_box)
    gmsh.model.mesh.setSize(
        gmsh.model.getEntitiesInBoundingBox(-a, -a, -a, a, a, a, dim=0), h_sphere
    )

    # Finding the inner spherical surface
    inner_surf = None
    for (dim, s) in surfs:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, s)
        if (abs((xmax - xmin) - 2 * a) < 1e-3 and
            abs((ymax - ymin) - 2 * a) < 1e-3 and
            abs((zmax - zmin) - 2 * a) < 1e-3):
            inner_surf = s
            break
    assert inner_surf is not None, "Inner spherical surface not found."

    gmsh.model.addPhysicalGroup(2, [inner_surf], tag=INNER_SPHERE, name="inner_sphere")

    # Bounding boxes for other facets
    eps = 1e-6 # tolerance (decrease if labels are not found)
    tags_map = {
        XMIN: (-half-eps, -half-eps, -half-eps, -half+eps,  half+eps,  half+eps),
        XMAX: ( half-eps, -half-eps, -half-eps,  half+eps,  half+eps,  half+eps),
        YMIN: (-half-eps, -half-eps, -half-eps,  half+eps, -half+eps,  half+eps),
        YMAX: (-half-eps,  half-eps, -half-eps,  half+eps,  half+eps,  half+eps),
        ZMIN: (-half-eps, -half-eps, -half-eps,  half+eps,  half+eps, -half+eps),
        ZMAX: (-half-eps, -half-eps,  half-eps,  half+eps,  half+eps,  half+eps),
    }

    for tag, bbox in tags_map.items():
        candidates = gmsh.model.getEntitiesInBoundingBox(*bbox, dim=2)
        ids = [s for (_, s) in candidates]
        if ids:
            gmsh.model.addPhysicalGroup(2, ids, tag=tag, name=FACET_LABELS[tag])

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    msh, cell_tags, facet_tags = gmshio.model_to_mesh(gmsh.model, comm, 0, gdim=3)
    gmsh.finalize()

    # Save
    if save:
        if comm.rank == 0:
            print(f"[INFO] Saving mesh to {basename}.xdmf and {basename}.vtu")

        # XDMF (mesh + markers)
        with XDMFFile(comm, f"{basename}.xdmf", "w") as xdmf:
            xdmf.write_mesh(msh)
            xdmf.write_meshtags(cell_tags, msh.geometry)
            xdmf.write_meshtags(facet_tags, msh.geometry)

        # VTU (per ParaView)
        # with VTXWriter(comm, f"{basename}.vtu", [msh.geometry]) as vtu:
        #     pass

    return msh, cell_tags, facet_tags

def plot_mesh(mesh, cell_tags=None, facet_tags=None, show_facets=True):
    """
    Visualizza una mesh dolfinx con pyvista.
    """
    topology, cell_types, geometry = dolfinx.plot.vtk_mesh(mesh, mesh.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=True, color="lightgray", opacity=0.5)

    if cell_tags is not None:
        grid.cell_data["cell_tags"] = cell_tags.values
        plotter.add_mesh(grid, scalars="cell_tags", show_edges=True, opacity=0.6)

    if facet_tags is not None and show_facets:
        topology_f, cell_types_f, geometry_f = dolfinx.plot.vtk_mesh(
            mesh, mesh.topology.dim - 1
        )
        grid_f = pyvista.UnstructuredGrid(topology_f, cell_types_f, geometry_f)

        # Costruisci array lungo quanto il numero di celle superficiali
        ncells = grid_f.n_cells
        values = np.zeros(ncells, dtype=np.int32)
        for idx, val in zip(facet_tags.indices, facet_tags.values):
            if idx < ncells:
                values[idx] = val

        grid_f.cell_data["facet_tags"] = values
        plotter.add_mesh(grid_f, scalars="facet_tags", show_edges=True, line_width=2)

        # Etichette leggibili
        unique = np.unique(values)
        print("\nFacet tags presenti:")
        for tag in unique:
            name = FACET_LABELS.get(int(tag), f"tag {tag}")
            print(f"  {tag} → {name}")

    plotter.add_axes()
    plotter.show()

if __name__ == "__main__":
    msh, cell_tags, facet_tags = build_mesh(L=1.0, a=0.2, save=True)
    plot_mesh(msh, cell_tags=cell_tags, facet_tags=facet_tags)
