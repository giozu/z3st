from mpi4py import MPI
import dolfinx
import ufl
import pyvista
import matplotlib.cm as cm
import numpy as np

class Mesh:
    def __init__(self, mesh_file, geometry, rank=0, gdim=3):

        """
        Initialize the class with the mesh file.

        Parameters:
        mesh_file (str): Path to the mesh file.
        rank (int): Rank of the process (default is 0).
        gdim (int): Geometric dimension of the mesh (default is 3).
        """
        print("__Mesh initializer__")

        self.mesh_file = mesh_file
        self.rank = rank
        self.gdim = gdim
        self.geometry = geometry
        
        self._read_mesh()

    def _read_mesh(self):
        """
        Read the mesh and facet tags from the specified mesh file.
        """
        print(f"Reading mesh file {self.mesh_file}")

        try:
            self.mesh, self.tags, self.ft = dolfinx.io.gmshio.read_from_msh(self.mesh_file, MPI.COMM_WORLD, rank=self.rank, gdim=self.gdim)
            self.tdim = self.mesh.topology.dim
            print(f"Info    : Mesh topology dimension d={self.tdim}")
            self.fdim = self.tdim - 1
            self.mesh.topology.create_connectivity(self.fdim, self.tdim)
            self.mesh.topology.create_connectivity(self.tdim, self.tdim) 
            self.boundary_facets = dolfinx.mesh.exterior_facet_indices(self.mesh.topology)

            print("Available volume tags (dx):")
            tag_values = self.tags.values
            unique_tags = set(tag_values)
            for tag in sorted(unique_tags):
                print(f"  Tag ID: {tag}")

            # Labels
            self.label_map = self.geometry["labels"]

            # print("\n--- Verifying facet tags ---")
            # unique_tags = np.unique(self.ft.values)
            # print(f"Unique tags found in mesh facet data (self.ft): {unique_tags}")
            # print(f"Interface tags from label_map: fuel_side={self.label_map['face_xx_1_1']}, clad_side={self.label_map['face_xx_0_2']}")

            print("Geometry entries:")
            for k in sorted(self.geometry.keys()):
                if k == "labels":
                    continue
                v = self.geometry[k]
                print(f"  {k:<15} → {v} ({type(v).__name__})")

            self.geometry_type = self.geometry.get("geometry_type", "").lower()

            self.height = self.geometry.get("height", None)
            print(f"  height = {self.height:.2f} m")

            if self.geometry_type == "rect":
                self.edge = self.geometry.get("edge", None)
                self.width = self.geometry.get("width", None)
                self.perimeter = self.edge * self.width
                self.area = self.edge * self.width
                print(f"  edge = {self.edge:.3f} m")
            elif self.geometry_type == "cyl":
                self.inner_radius = self.geometry.get("inner_radius", None)
                self.outer_radius = self.geometry.get("outer_radius", None)
                self.perimeter = 2. * np.pi * self.outer_radius
                self.area = np.pi * (self.outer_radius**2 - self.inner_radius**2)
                print(f"  inner_radius = {self.inner_radius:.2e} m")
                print(f"  outer_radius = {self.outer_radius:.2e} m")
            elif self.geometry_type == "cyl-cyl":
                self.fuel_inner_radius = self.geometry.get("fuel_inner_radius", None)
                self.fuel_outer_radius = self.geometry.get("fuel_outer_radius", None)
                self.cladding_inner_radius = self.geometry.get("cladding_inner_radius", None)
                self.cladding_outer_radius = self.geometry.get("cladding_outer_radius", None)
                self.perimeter = 2. * np.pi * self.fuel_outer_radius
                self.area = np.pi * (self.fuel_outer_radius**2 - self.fuel_inner_radius**2)
                print(f"  fuel_inner_radius = {self.fuel_inner_radius:.2e} m")
                print(f"  fuel_outer_radius = {self.fuel_outer_radius:.2e} m")
                print(f"  cladding_inner_radius = {self.cladding_inner_radius:.2e} m")
                print(f"  cladding_outer_radius = {self.cladding_outer_radius:.2e} m")

            print(f"  area = {self.area:.3e} m2")
            print(f"Label map from '{self.config["geometry_path"]}'")
            for label, tag in sorted(self.label_map.items(), key=lambda item: item[1]):
                print(f"  {label:<12} → {tag}")
            
            self.normal = ufl.FacetNormal(self.mesh)
            self.dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
            self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft)

            # fuel_face_tag = self.label_map["face2_1"]      # fuel side of the interface
            # clad_face_tag = self.label_map["face4_2"]      # cladding side

            # Define measures for each surface
            # self.ds_fuel_face = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=fuel_face_tag)
            # self.ds_clad_face = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft, subdomain_id=clad_face_tag)

            # self.ds_ci = self.ds(self.cladding_inner_label)
            # self.ds_co = self.ds(self.cladding_outer_label)
            # self.ds_fo = self.ds(self.fuel_outer_label)

            # msh = meshio.read(self.mesh_file)

            # facet_tags = self.tags
            # boundary_facets = np.where(facet_tags.values > 0)[0]
            # boundary_values = facet_tags.values[boundary_facets]

            # print("Boundary facets:", boundary_facets)
            # print("Physical tags:", boundary_values)

            # print(msh.cells_dict)
            # print(msh.cell_data_dict)

            # for cell_block in msh.cells:
            #     if cell_block.type == "triangle":
            #         print("Triangular faces:", cell_block.data)
            #         print("Face tags:", msh.cell_data_dict["gmsh:physical"]["triangle"])

        except Exception as e:
            raise ValueError(f"Error reading mesh file {self.mesh_file}: {e}")
        
    def locateFacetDofs(self, label: int, V: dolfinx.fem.functionspace) -> np.ndarray:
        """Locate degrees of freedom (DoFs) on a specific facet boundary.

        This is a convenience function for finding all DoFs associated with
        a boundary marked by a specific integer label. It's commonly used
        to apply boundary conditions (e.g., DirichletBC).

        Args:
            label: The integer marker for the desired facet set (e.g., a
                   specific boundary like an "inlet" or "wall"). This label must
                   exist in the mesh's facet tags.
            V: The function space from which to extract the DoFs.

        Returns:
            A 1D numpy array of DoF indices that are located on the
            specified facets.
        """
        # Find the indices of all facets (e.g., edges in 2D, faces in 3D)
        # that are marked with the given integer 'label'.
        # self.ft is the MeshTags object for facets.
        facets = self.ft.find(label)

        # Locate the DoFs that are topologically associated with these facets.
        # - V: The function space defines the DoFs.
        # - self.fdim: The topological dimension of the entities to search on (dim-1).
        # - facets: The list of entity indices to search.
        return dolfinx.fem.locate_dofs_topological(V, self.fdim, facets)

    def locateDomainDofs(self, label_map: dict, material_label: str, V: dolfinx.fem.functionspace) -> np.ndarray:
        """Locate degrees of freedom (DoFs) within a specific subdomain.

        This function identifies all DoFs that belong to cells of a certain
        type (e.g., a specific material). It is useful for defining
        material-dependent properties, source terms, or for subdomain-specific
        post-processing.

        Args:
            label_map: A dictionary mapping human-readable material labels
                       (e.g., "water", "steel") to their corresponding integer
                       tags in the mesh's cell data.
            material_label: The human-readable key (e.g., "water") to look
                            up the integer tag in `label_map`.
            V: The function space from which to extract the DoFs.

        Returns:
            A 1D numpy array of DoF indices located within the cells of the
            specified subdomain.
        """
        # Use the map to get the integer tag corresponding to the material label.
        # This decouples the physics definition from the mesh generation details.
        domain_tag = label_map[material_label]

        # Find the indices of all cells (e.g., triangles in 2D, tetrahedra in 3D)
        # that are marked with that integer domain tag.
        # self.tags is the MeshTags object for cells.
        cells = self.tags.find(domain_tag)

        # Locate the DoFs that are topologically associated with these cells.
        # - V: The function space defines the DoFs.
        # - self.tdim: The topological dimension of the entities (the cells themselves).
        # - cells: The list of cell indices to search.
        dofs = dolfinx.fem.locate_dofs_topological(V, self.tdim, cells)
        return dofs

    def volume(self, log_path="output/simulation_log.txt"):

        # Integration measures the domain
        for label in ["fuel", "cladding"]:
            if label in self.label_map:
                tag = self.label_map[label]
                setattr(self, f"dx_{label}", self.dx(tag))
                print(f"[INFO] Defined dx_{label} for tag ID {tag}")

        # Compute mesh volumes
        volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx))
        volume_1 = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_fuel))
        volume_2 = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_cladding))
        
        with open(log_path, "a") as f:
            f.write("Mesh Volumes\n")
            f.write("──────────────────────────────────────────────\n")
            f.write(f"{'Volume  (from mesh)':<32}: {volume  *1e9:.3e} mm³\n")
            f.write(f"{'Volume1 (from mesh)':<32}: {volume_1*1e9:.3e} mm³\n")
            f.write(f"{'Volume2 (from mesh)':<32}: {volume_2*1e9:.3e} mm³\n")

    # def volume(self, edge, log_path="output/simulation_log.txt"):
    #     self.edge = edge
    #     self.area = edge**2
    #     self.height = edge

    #     mesh_volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx))
    #     analytical_volume = self.area * self.height
    #     rel_error = abs(mesh_volume - analytical_volume) / analytical_volume * 100

    #     with open(log_path, "a") as f:
    #         f.write("\nVolume Comparison\n")
    #         f.write(f"Mesh volume: {mesh_volume:.6e} m3\n")
    #         f.write(f"Analytical volume: {analytical_volume:.6e} m3\n")
    #         f.write(f"Relative error: {rel_error:.3e} %\n")

    # def volume(self, fuel_outer_radius, cladding_inner_radius, cladding_outer_radius, height, log_path="output/simulation_log.txt"):
    #     # Integration measures the domain
    #     dx_subdomain = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
    #     self.dx_fuel = dx_subdomain(1)      # Marker ID for fuel region
    #     self.dx_cladding = dx_subdomain(2)  # Marker ID for cladding region

    #     # Compute mesh volumes
    #     volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx_subdomain))
    #     volume_fuel = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_fuel))
    #     volume_cladding = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * self.dx_cladding))
        
    #     analytical_volume = (self.area_fuel + self.area_cladding)*self.height
    #     rel_error = abs(volume - analytical_volume) / analytical_volume * 100

    #     with open(log_path, "a") as f:
    #         f.write("Mesh Volumes\n")
    #         f.write("──────────────────────────────────────────────\n")
    #         f.write(f"{'Fuel volume (from mesh)':<32}: {volume_fuel*1e9:.3e} mm³\n")
    #         f.write(f"{'Cladding volume (from mesh)':<32}: {volume_cladding*1e9:.3e} mm³\n")
    #         f.write(f"{'Fuel volume (analytical)':<32}: {self.area_fuel * self.height * 1e9:.3e} mm³\n")
    #         f.write(f"{'Cladding volume (analytical)':<32}: {self.area_cladding * self.height * 1e9:.3e} mm³\n")
    #         f.write(f"{'Relative error':<32}: {rel_error:.3e} %\n\n")

    def _mesh_plot(self, mesh_plot_mode):
        # from dolfinx.io import XDMFFile
        
        # self.mesh.topology.create_connectivity(self.fdim , self.tdim)
        # with XDMFFile(self.mesh.comm, "facet_tags.xdmf", "w") as xdmf:
        #     xdmf.write_mesh(self.mesh)
        #     xdmf.write_meshtags(self.tags, self.geometry)

        if mesh_plot_mode is False:
            return
        else:
            topology, cell_types, _ = dolfinx.plot.vtk_mesh(self.mesh, self.mesh.topology.dim - 1)
            surface = pyvista.UnstructuredGrid(topology, cell_types, self.mesh.geometry.x)

            unique_tags = np.unique(self.ft.values)
            print("[INFO] Face labels present:", unique_tags)

            colors = cm.get_cmap("tab10", len(unique_tags))
            label_map = {}
            for i, tag in enumerate(unique_tags):
                label_name = list(self.label_map.keys())[list(self.label_map.values()).index(tag)]
                label_map[tag] = (label_name, colors(i)[:3])

            plotter = pyvista.Plotter()
            for tag, (name, color) in label_map.items():
                facets = self.ft.find(tag)
                cells = surface.extract_cells(facets)
                plotter.add_mesh(cells, color=color, show_edges=True, label=name)

            plotter.add_legend()
            plotter.show()
