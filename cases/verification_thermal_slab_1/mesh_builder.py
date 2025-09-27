import gmsh
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile
import dolfinx
import yaml
import pyvista
import numpy as np


class MeshBuilder:
    def __init__(self, geometry_file: str, comm=MPI.COMM_WORLD):
        with open(geometry_file, "r") as f:
            self.geometry = yaml.safe_load(f)

        self.comm = comm
        self.Lx = float(self.geometry.get("Lx", 1.0))
        self.Ly = float(self.geometry.get("Ly", 0.1))
        self.Lz = float(self.geometry.get("Lz", 0.1))
        self.h_box = float(self.geometry.get("h_box", 0.05))
        self.order = int(self.geometry.get("order", 1))
        self.label_map = self.geometry.get("labels", {})
        self.inv_label_map = {v: k for k, v in self.label_map.items()}
        self.plot_mesh = self.geometry.get("plot", False)

        self.mesh = None
        self.cell_tags = None
        self.facet_tags = None

    def build(self):
        gmsh.model.add("box")

        # Box centrato all'origine
        box = gmsh.model.occ.addBox(
            -self.Lx / 2, -self.Ly / 2, -self.Lz / 2,
             self.Lx, self.Ly, self.Lz
        )
        gmsh.model.occ.synchronize()

        # n_cells = 20
        # for e in gmsh.model.getEntities(1):
        #     gmsh.model.mesh.setTransfiniteCurve(e[1], n_cells+1)

        # for s in gmsh.model.getEntities(2):
        #     gmsh.model.mesh.setTransfiniteSurface(s[1])
        #     gmsh.model.mesh.setRecombine(2, s[1])

        # for v in gmsh.model.getEntities(3):
        #     gmsh.model.mesh.setTransfiniteVolume(v[1])
        #     gmsh.model.mesh.setRecombine(3, v[1])

        # gmsh.model.mesh.setTransfiniteVolume(box)   # usa direttamente il tag del box
        # gmsh.model.mesh.setRecombine(3, box)

        # Volume solido
        gmsh.model.addPhysicalGroup(3, [box], tag=self.label_map["solid"], name="solid")

        # Bounding boxes per tutte le facce
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
                gmsh.model.addPhysicalGroup(2, ids, tag=tag, name=self.inv_label_map[tag])

        # Mesh
        gmsh.model.mesh.setSize(gmsh.model.getEntities(dim=0), self.h_box)
        gmsh.model.mesh.setOrder(self.order)

        # gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        # gmsh.option.setNumber("Mesh.RecombineAll", 1)

        gmsh.model.mesh.generate(3)

        types = gmsh.model.mesh.getElementTypes(3)
        print("Element types in 3D volume:", types)

        # Converto in mesh Dolfinx
        self.mesh, self.cell_tags, self.facet_tags = gmshio.model_to_mesh(
            gmsh.model, self.comm, 0, gdim=3
        )

        return self.mesh, self.cell_tags, self.facet_tags

    def save(self, basename="mesh"):
        if self.comm.rank == 0:
            gmsh.write(f"{basename}.msh")
            print(f"[INFO] Saved mesh → {basename}.msh")

        with XDMFFile(self.comm, f"{basename}.xdmf", "w") as xdmf:
            xdmf.write_mesh(self.mesh)
            xdmf.write_meshtags(self.cell_tags, self.mesh.geometry)
            xdmf.write_meshtags(self.facet_tags, self.mesh.geometry)

    def plot(self):

        if not self.plot_mesh:
            return

        if self.comm.rank != 0:
            return

        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(self.mesh, self.mesh.topology.dim)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        plotter = pyvista.Plotter()
        plotter.add_mesh(grid, show_edges=True, color="lightgray", opacity=0.3)

        if self.facet_tags is not None:
            topology_f, cell_types_f, geometry_f = dolfinx.plot.vtk_mesh(
                self.mesh, self.mesh.topology.dim - 1
            )
            grid_f = pyvista.UnstructuredGrid(topology_f, cell_types_f, geometry_f)

            ncells = grid_f.n_cells
            values = np.zeros(ncells, dtype=np.int32)
            for idx, val in zip(self.facet_tags.indices, self.facet_tags.values):
                if idx < ncells:
                    values[idx] = val

            color_map = {
                self.label_map["xmin"]: "red",
                self.label_map["xmax"]: "blue",
                self.label_map["ymin"]: "green",
                self.label_map["ymax"]: "yellow",
                self.label_map["zmin"]: "magenta",
                self.label_map["zmax"]: "cyan",
            }

            for tag, color in color_map.items():
                mask = values == tag
                if np.any(mask):
                    cells = grid_f.extract_cells(np.where(mask)[0])
                    plotter.add_mesh(
                        cells, color=color, show_edges=True,
                        label=self.inv_label_map.get(tag, f"tag {tag}")
                    )

            plotter.add_legend()

            unique = np.unique(values)
            print("\nFacet tags present:")
            for tag in unique:
                name = self.inv_label_map.get(int(tag), f"tag {tag}")
                print(f"  {tag} → {name}")

        plotter.add_axes()
        plotter.show()


if __name__ == "__main__":
    gmsh.initialize()
    builder = MeshBuilder("geometry.yaml")
    msh, cell_tags, facet_tags = builder.build()
    builder.save("mesh")
    builder.plot()
    gmsh.finalize()
