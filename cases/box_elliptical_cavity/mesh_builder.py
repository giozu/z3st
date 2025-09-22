import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile
import dolfinx
import pyvista
import yaml

# Facet markers (now defined in geometry.yaml)
INNER_SPHERE = 1
XMIN, XMAX = 11, 12
YMIN, YMAX = 13, 14
ZMIN, ZMAX = 15, 16

FACET_LABELS = {
    INNER_SPHERE: "inner_sphere",
    XMIN: "xmin",
    XMAX: "xmax",
    YMIN: "ymin",
    YMAX: "ymax",
    ZMIN: "zmin",
    ZMAX: "zmax",
}

class MeshBuilder:
    def __init__(self, geometry_file: str, comm=MPI.COMM_WORLD):
        """
        Build a box with spherical cavity from YAML description.
        """
        with open(geometry_file, "r") as f:
            self.geometry = yaml.safe_load(f)

        self.comm = comm

        # Dimensions
        self.Lx = float(self.geometry.get("Lx", 1.0))
        self.Ly = float(self.geometry.get("Ly", 1.0))
        self.Lz = float(self.geometry.get("Lz", 1.0))
        self.a = float(self.geometry.get("a", 0.1))

        # Mesh params
        self.h_box = float(self.geometry.get("h_box", 0.05))
        self.h_sphere = float(self.geometry.get("h_sphere", 0.01))
        self.order = int(self.geometry.get("order", 1))

        # Labels
        self.label_map = self.geometry.get("labels", {})

        self.mesh = None
        self.cell_tags = None
        self.facet_tags = None

    def build(self):
        gmsh.initialize()
        gmsh.model.add("box_bubble")

        # Box centered at origin
        box = gmsh.model.occ.addBox(-self.Lx/2, -self.Ly/2, -self.Lz/2,
                                    self.Lx, self.Ly, self.Lz)

        # Sezione: cavità (sfera o ellissoide)
        self.ax = float(self.geometry.get("ax", self.a))
        self.ay = float(self.geometry.get("ay", self.a))
        self.az = float(self.geometry.get("az", self.a))

        # Creo sfera unitaria e la dilato
        sph = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
        gmsh.model.occ.dilate([(3, sph)], 0, 0, 0, self.ax, self.ay, self.az)

        # Cut: box - ellissoide
        vol, _ = gmsh.model.occ.cut([(3, box)], [(3, sph)],
                                    removeObject=True, removeTool=False)
        gmsh.model.occ.synchronize()

        gmsh.model.addPhysicalGroup(3, [vol[0][1]],
                                    tag=self.label_map["solid"], name="solid")

        # Mesh size
        gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), self.h_box)
        gmsh.model.mesh.setSize(
            gmsh.model.getEntitiesInBoundingBox(-self.ax, -self.ay, -self.az,
                                                self.ax, self.ay, self.az, dim=0),
            self.h_sphere
        )

        inner_surf = None
        for (dim, s) in gmsh.model.getEntities(dim=2):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, s)
            if xmin < 0 and xmax > 0 and ymin < 0 and ymax > 0 and zmin < 0 and zmax > 0:
                inner_surf = s
                break

        assert inner_surf is not None, "Inner ellipsoidal surface not found."

        gmsh.model.addPhysicalGroup(2, [inner_surf],
                                    tag=self.label_map["inner_sphere"],
                                    name="inner_sphere")

        # Add box face groups from label_map
        eps = 1e-6
        tags_map = {
            self.label_map["xmin"]: (-self.Lx/2-eps, -self.Ly/2-eps, -self.Lz/2-eps,
                                    -self.Lx/2+eps,  self.Ly/2+eps,  self.Lz/2+eps),
            self.label_map["xmax"]: ( self.Lx/2-eps, -self.Ly/2-eps, -self.Lz/2-eps,
                                    self.Lx/2+eps,  self.Ly/2+eps,  self.Lz/2+eps),
            self.label_map["ymin"]: (-self.Lx/2-eps, -self.Ly/2-eps, -self.Lz/2-eps,
                                    self.Lx/2+eps, -self.Ly/2+eps,  self.Lz/2+eps),
            self.label_map["ymax"]: (-self.Lx/2-eps,  self.Ly/2-eps, -self.Lz/2-eps,
                                    self.Lx/2+eps,  self.Ly/2+eps,  self.Lz/2+eps),
            self.label_map["zmin"]: (-self.Lx/2-eps, -self.Ly/2-eps, -self.Lz/2-eps,
                                    self.Lx/2+eps,  self.Ly/2+eps, -self.Lz/2+eps),
            self.label_map["zmax"]: (-self.Lx/2-eps, -self.Ly/2-eps,  self.Lz/2-eps,
                                    self.Lx/2+eps,  self.Ly/2+eps,  self.Lz/2+eps),
        }

        for tag, bbox in tags_map.items():
            candidates = gmsh.model.getEntitiesInBoundingBox(*bbox, dim=2)
            ids = [s for (_, s) in candidates]
            if ids:
                gmsh.model.addPhysicalGroup(2, ids, tag=tag,
                                            name=[k for k, v in self.label_map.items() if v == tag][0])

        gmsh.model.mesh.setOrder(self.order)
        gmsh.model.mesh.generate(3)

        self.mesh, self.cell_tags, self.facet_tags = gmshio.model_to_mesh(
            gmsh.model, self.comm, 0, gdim=3)

        return self.mesh, self.cell_tags, self.facet_tags

    def save(self, basename="mesh", export_msh=True, export_xdmf=True):
        """Save mesh in .msh (Gmsh) and/or .xdmf (FEniCSx) format."""
        if export_msh:
            if self.comm.rank == 0:
                gmsh.write(f"{basename}.msh")
                print(f"[INFO] Saved native Gmsh mesh → {basename}.msh")

        if export_xdmf:
            if self.comm.rank == 0:
                print(f"[INFO] Saving mesh to {basename}.xdmf")

            with XDMFFile(self.comm, f"{basename}.xdmf", "w") as xdmf:
                xdmf.write_mesh(self.mesh)
                xdmf.write_meshtags(self.cell_tags, self.mesh.geometry)
                xdmf.write_meshtags(self.facet_tags, self.mesh.geometry)

        gmsh.finalize()

    def plot(self, show_facets=True):
        """Quick PyVista visualization."""
        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(self.mesh, self.mesh.topology.dim)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        plotter = pyvista.Plotter()
        plotter.add_mesh(grid, show_edges=True, color="lightgray", opacity=0.5)

        if self.facet_tags is not None and show_facets:
            topology_f, cell_types_f, geometry_f = dolfinx.plot.vtk_mesh(
                self.mesh, self.mesh.topology.dim - 1
            )
            grid_f = pyvista.UnstructuredGrid(topology_f, cell_types_f, geometry_f)

            ncells = grid_f.n_cells
            values = np.zeros(ncells, dtype=np.int32)
            for idx, val in zip(self.facet_tags.indices, self.facet_tags.values):
                if idx < ncells:
                    values[idx] = val

            grid_f.cell_data["facet_tags"] = values
            plotter.add_mesh(grid_f, scalars="facet_tags", show_edges=True, line_width=2)

            unique = np.unique(values)
            print("\nFacet tags present:")
            for tag in unique:
                name = FACET_LABELS.get(int(tag), f"tag {tag}")
                print(f"  {tag} → {name}")

        plotter.add_axes()
        plotter.show()

if __name__ == "__main__":
    builder = MeshBuilder("geometry.yaml")
    msh, cell_tags, facet_tags = builder.build()
    builder.save("mesh", export_msh=True, export_xdmf=True)
    builder.plot()
