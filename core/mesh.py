from mpi4py import MPI
import dolfinx
import ufl

class Mesh:
    def __init__(self, mesh_file, rank=0, gdim=3):
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

            # Labels
            # --- Cube ---
            self.label_map = {
                    "Bottom": 1,
                    "Face 1": 2,
                    "Face 2": 3,
                    "Face 3": 4,
                    "Face 4": 5,
                    "Top": 6,
                }

            # --- Cylindrical inclusion ---
            label_map = {
                    "Domain": 1,
                    "Inclusion": 2,
                    "Boundary": 3,
                }
            
            self.normal = ufl.FacetNormal(self.mesh)

            self.dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
            self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft)

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