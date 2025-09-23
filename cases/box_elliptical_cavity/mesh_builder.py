import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile
import dolfinx
import pyvista
import yaml

class MeshBuilder:
    def __init__(self, geometry_file: str, comm=MPI.COMM_WORLD):
        """
        Build a box with spherical or ellipsoidal cavity from YAML description.
        """
        with open(geometry_file, "r") as f:
            self.geometry = yaml.safe_load(f)

        self.comm = comm

        # Dimensions
        self.Lx = float(self.geometry.get("Lx", 1.0))
        self.Ly = float(self.geometry.get("Ly", 1.0))
        self.Lz = float(self.geometry.get("Lz", 1.0))

        # Cavità (sfera di default, ellissoide se ax, ay, az specificati)
        self.a = float(self.geometry.get("a", 0.1))
        self.ax = float(self.geometry.get("ax", self.a))
        self.ay = float(self.geometry.get("ay", self.a))
        self.az = float(self.geometry.get("az", self.a))

        # Mesh params
        self.h_box = float(self.geometry.get("h_box", 0.05))
        self.h_sphere = float(self.geometry.get("h_sphere", 0.01))
        self.order = int(self.geometry.get("order", 1))

        # Labels
        self.label_map = self.geometry.get("labels", {})
        self.inv_label_map = {v: k for k, v in self.label_map.items()}

        # Containers
        self.mesh = None
        self.cell_tags = None
        self.facet_tags = None

    def build(self):
        gmsh.initialize()
        gmsh.model.add("box_cavity")

        # Box centrato all'origine
        box = gmsh.model.occ.addBox(
            -self.Lx/2, -self.Ly/2, -self.Lz/2,
             self.Lx,   self.Ly,   self.Lz
        )

        # Creo sfera unitaria e la dilato in ellissoide
        sph = gmsh.model.occ.addSphere(0, 0, 0, 1.0)
        gmsh.model.occ.dilate([(3, sph)], 0, 0, 0, self.ax, self.ay, self.az)

        # Cut: box - ellissoide
        vol, _ = gmsh.model.occ.cut([(3, box)], [(3, sph)],
                                    removeObject=True, removeTool=False)
        gmsh.model.occ.synchronize()

        # Volume solido
        gmsh.model.addPhysicalGroup(3, [vol[0][1]],
                                    tag=self.label_map["solid"],
                                    name="solid")

        # Mesh size
        gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), self.h_box)
        gmsh.model.mesh.setSize(
            gmsh.model.getEntitiesInBoundingBox(
                -self.ax, -self.ay, -self.az,
                 self.ax,  self.ay,  self.az, dim=0),
            self.h_sphere
        )

        # Cavità: prendo le superfici dal tag dell’ellissoide
        cavity_surfaces = gmsh.model.getBoundary([(3, sph)], oriented=False, recursive=False)
        cavity_surface_tags = [s for (dim, s) in cavity_surfaces if dim == 2]
        assert cavity_surface_tags, "Ellipsoidal surface not found."

        gmsh.model.addPhysicalGroup(2, cavity_surface_tags,
                                    tag=self.label_map["inner_sphere"],
                                    name="inner_sphere")

        # Facce del box
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
                                            name=self.inv_label_map[tag])

        # Mesh
        gmsh.model.mesh.setOrder(self.order)
        gmsh.model.mesh.generate(3)

        self.mesh, self.cell_tags, self.facet_tags = gmshio.model_to_mesh(
            gmsh.model, self.comm, 0, gdim=3)

        return self.mesh, self.cell_tags, self.facet_tags

    def save(self, basename="mesh", export_msh=True, export_xdmf=True):
        """Save mesh in .msh (Gmsh) and/or .xdmf (FEniCSx) format."""
        if export_msh and self.comm.rank == 0:
            gmsh.write(f"{basename}.msh")
            print(f"[INFO] Saved native Gmsh mesh → {basename}.msh")

        if export_xdmf:
            if self.comm.rank == 0:
                print(f"[INFO] Saving mesh to {basename}.xdmf")

            with XDMFFile(self.comm, f"{basename}.xdmf", "w") as xdmf:
                xdmf.write_mesh(self.mesh)
                xdmf.write_meshtags(self.cell_tags, self.mesh.geometry)
                xdmf.write_meshtags(self.facet_tags, self.mesh.geometry)

        gmsh.finalize()   # chiudo gmsh qui

    def plot(self, show_facets=True):
        """Quick PyVista visualization (rank 0 only)."""
        if self.comm.rank != 0:
            return

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
                name = self.inv_label_map.get(int(tag), f"tag {tag}")
                print(f"  {tag} → {name}")

        plotter.add_axes()
        plotter.show()


if __name__ == "__main__":
    builder = MeshBuilder("geometry.yaml")
    msh, cell_tags, facet_tags = builder.build()
    builder.save("mesh", export_msh=True, export_xdmf=True)
    builder.plot()
