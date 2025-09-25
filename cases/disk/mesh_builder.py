import gmsh
import numpy as np
import yaml
from mpi4py import MPI
from dolfinx.io import gmshio, XDMFFile
import dolfinx
import pyvista


class CylinderMeshBuilder:
    def __init__(self, geometry_file: str, comm=MPI.COMM_WORLD):
        """
        Structured cylinder mesh builder.
        Reads parameters from a YAML geometry file.
        """
        print("\n[CylinderMeshBuilder] Initializing from", geometry_file)

        with open(geometry_file, "r") as f:
            self.geometry = yaml.safe_load(f)

        self.comm = comm

        # Dimensions
        self.r = float(self.geometry.get("r", 0.004))
        self.Lz = float(self.geometry.get("Lz", 0.01))

        # Mesh params
        self.n_theta = int(self.geometry.get("n_theta", 64))
        self.n_z = int(self.geometry.get("n_z", 20))
        self.order = int(self.geometry.get("order", 1))

        # Labels
        self.label_map = self.geometry.get("labels", {})
        self.inv_label_map = {v: k for k, v in self.label_map.items()}

        # Containers
        self.mesh = None
        self.cell_tags = None
        self.facet_tags = None

        print(f"  → radius       : {self.r:.4e} m")
        print(f"  → height       : {self.Lz:.4e} m")
        print(f"  → n_theta      : {self.n_theta}")
        print(f"  → n_z          : {self.n_z}")
        print(f"  → order        : {self.order}")
        print(f"  → labels       : {self.label_map}")

    def build(self):
        print("\n[CylinderMeshBuilder] Building geometry...")

        gmsh.initialize()
        gmsh.model.add("cylinder")

        # Base points
        p0 = gmsh.model.occ.addPoint(0, 0, 0)
        p1 = gmsh.model.occ.addPoint(self.r, 0, 0)
        p2 = gmsh.model.occ.addPoint(0, self.r, 0)
        p3 = gmsh.model.occ.addPoint(-self.r, 0, 0)
        p4 = gmsh.model.occ.addPoint(0, -self.r, 0)

        # Circle arcs
        c1 = gmsh.model.occ.addCircleArc(p1, p0, p2)
        c2 = gmsh.model.occ.addCircleArc(p2, p0, p3)
        c3 = gmsh.model.occ.addCircleArc(p3, p0, p4)
        c4 = gmsh.model.occ.addCircleArc(p4, p0, p1)

        loop = gmsh.model.occ.addCurveLoop([c1, c2, c3, c4])
        surf = gmsh.model.occ.addPlaneSurface([loop])

        # Extrusion along z
        print(f"  Extruding surface along z = {self.Lz:.4e} m with {self.n_z} layers...")
        out = gmsh.model.occ.extrude([(2, surf)], 0, 0, self.Lz,
                                     numElements=[self.n_z], recombine=True)
        gmsh.model.occ.synchronize()

        # --- Physical groups ---
        print("\n  Defining physical groups...")
        if "bottom" in self.label_map:
            gmsh.model.addPhysicalGroup(2, [surf],
                                        tag=self.label_map["bottom"],
                                        name="bottom")
            print(f"    bottom   → tag {self.label_map['bottom']}")

        if "top" in self.label_map:
            gmsh.model.addPhysicalGroup(2, [out[0][1]],
                                        tag=self.label_map["top"],
                                        name="top")
            print(f"    top      → tag {self.label_map['top']}")

        lateral_surfaces = [ent[1] for ent in out[2:] if ent[0] == 2]
        if "lateral" in self.label_map:
            gmsh.model.addPhysicalGroup(2, lateral_surfaces,
                                        tag=self.label_map["lateral"],
                                        name="lateral")
            print(f"    lateral  → tag {self.label_map['lateral']} (n={len(lateral_surfaces)} lateral faces)")

        if "fuel" in self.label_map:
            gmsh.model.addPhysicalGroup(3, [out[1][1]],
                                        tag=self.label_map["fuel"],
                                        name="fuel")
            print(f"    fuel     → tag {self.label_map['fuel']}")

        # --- Mesh controls ---
        print("\n  Applying transfinite/recombine settings...")
        n_quarter = self.n_theta // 4
        # for c in [c1, c2, c3, c4]:
        #     gmsh.model.mesh.setTransfiniteCurve(c, n_quarter + 1)
        # gmsh.model.mesh.setTransfiniteSurface(surf, cornerTags=[p1, p2, p3, p4])
        gmsh.model.mesh.setRecombine(2, surf)

        print(f"  Generating mesh (order {self.order})...")
        gmsh.model.mesh.setOrder(self.order)
        gmsh.model.mesh.generate(3)

        # gmsh.model.mesh.refine()
        # gmsh.model.mesh.refine()

        # Convert to dolfinx
        self.mesh, self.cell_tags, self.facet_tags = gmshio.model_to_mesh(
            gmsh.model, self.comm, 0, gdim=3
        )

        print("[CylinderMeshBuilder] Mesh built successfully.")
        return self.mesh, self.cell_tags, self.facet_tags

    def save(self, basename="mesh", export_msh=True, export_xdmf=True):
        print(f"\n[CylinderMeshBuilder] Saving mesh with basename '{basename}'")

        if export_msh and self.comm.rank == 0:
            gmsh.write(f"{basename}.msh")
            print(f"  → Saved gmsh mesh → {basename}.msh")

        if export_xdmf:
            with XDMFFile(self.comm, f"{basename}.xdmf", "w") as xdmf:
                xdmf.write_mesh(self.mesh)
                xdmf.write_meshtags(self.cell_tags, self.mesh.geometry)
                xdmf.write_meshtags(self.facet_tags, self.mesh.geometry)
            print(f"  → Saved XDMF mesh → {basename}.xdmf")

        gmsh.finalize()

    def plot(self):
        if self.comm.rank != 0:
            return
        print("\n[CylinderMeshBuilder] Plotting mesh (pyvista)...")
        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(self.mesh, self.mesh.topology.dim)
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        plotter = pyvista.Plotter()
        plotter.add_mesh(grid, show_edges=True, color="lightgray", opacity=0.5)
        plotter.show()


if __name__ == "__main__":
    builder = CylinderMeshBuilder("geometry.yaml")
    msh, cell_tags, facet_tags = builder.build()
    builder.save("mesh", export_msh=True, export_xdmf=True)
    builder.plot()
