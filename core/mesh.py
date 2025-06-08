from mpi4py import MPI
import dolfinx
import ufl

class Mesh:
    def __init__(self, mesh_file, geometry, rank=0, gdim=3):
        """
        Initialize the class with the mesh file.

        Parameters:
        mesh_file (str): Path to the mesh file.
        rank (int): Rank of the process (default is 0).
        gdim (int): Geometric dimension of the mesh (default is 3).
        """
        self.mesh_file = mesh_file
        self.rank = rank
        self.gdim = gdim
        self.geometry = geometry
        
        self._read_mesh()
        self._setup_function_space()
        self._initialize_functions()

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
            self.boundary_facets = dolfinx.mesh.exterior_facet_indices(self.mesh.topology)

            print("[INFO] Available volume tags (dx):")
            tag_values = self.tags.values
            unique_tags = set(tag_values)
            for tag in sorted(unique_tags):
                print(f"  Tag ID: {tag}")


            # Labels
            self.label_map = self.geometry["labels"]

            print("[INFO] Geometry entries:")
            for k in sorted(self.geometry.keys()):
                if k == "labels":
                    continue
                v = self.geometry[k]
                print(f"  {k:<15} → {v} ({type(v).__name__})")

            self.edge = self.geometry.get("edge", None)
            self.height = self.geometry.get("height", None)
            self.geometry_type = self.geometry.get("geometry_type", "").lower()

            if self.geometry_type == "rect":
                self.area = self.edge * 2
            else:
                self.area = self.edge ** 2
            print(f"[INFO] edge = {self.edge} m, height = {self.height} m")
            print(f"[INFO] area = {self.area} m2")

            print("[INFO] Label map:")
            for label, tag in sorted(self.label_map.items(), key=lambda item: item[1]):
                print(f"  {label:<12} → {tag}")
            
            self.normal = ufl.FacetNormal(self.mesh)
            self.dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
            self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft)

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
        
    def locateDofs(self, label):
        return dolfinx.fem.locate_dofs_topological(self.V_m, self.fdim, self.ft.find(label))
    
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