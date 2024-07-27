from __future__ import print_function

from mpi4py import MPI
from dolfinx.io import XDMFFile, gmshio
from dolfinx.mesh import meshtags_from_entities
from dolfinx.fem.petsc import LinearProblem
from dolfinx import fem, plot
from dolfinx import default_scalar_type

import ufl
import matplotlib.pyplot as plt
import numpy as np
import ufl.measure

try:
    import gmsh
except ImportError:
    print("This demo requires gmsh to be installed")
    exit(0)

Re = 11.
Ri = 9.

def gmsh_hollow_circle(model: gmsh.model, name: str, option: gmsh.option) -> gmsh.model:
    model.add(name)
    model.setCurrent(name)

    # Create outer and inner circles
    outer_circle = model.occ.addDisk(0, 0, 0, Re, Re)
    inner_circle = model.occ.addDisk(0, 0, 0, Ri, Ri)
    
    # Cut the inner circle from the outer circle to create a hollow circle
    hollow_circle = model.occ.cut([(2, outer_circle)], [(2, inner_circle)])[0]
    model.occ.synchronize()

    # Add physical groups
    gdim = 2
    model.addPhysicalGroup(gdim, [hollow_circle[0][1]], tag=1)
    model.setPhysicalName(2, 1, "Hollow circle")
    
    boundary_entities = model.getEntities(2)
    boundary_ids = [entity[1] for entity in boundary_entities if entity[1] != hollow_circle[0][1]]
    model.addPhysicalGroup(2, boundary_ids, tag=2)
    model.setPhysicalName(2, 2, "Boundary")

    option.setNumber("Mesh.CharacteristicLengthMin", 0.5)
    option.setNumber("Mesh.CharacteristicLengthMax", 0.5)
    model.mesh.generate(gdim)
    return model

def create_mesh(comm: MPI.Comm, model: gmsh.model, name: str, filename: str, mode: str):
    msh, ct, ft = gmshio.model_to_mesh(model, comm, rank=0)
    msh.name = name
    ct.name = f"{msh.name}_cells"
    ft.name = f"{msh.name}_facets"

    tdim = msh.topology.dim
    fdim = tdim - 1

    with XDMFFile(msh.comm, filename, mode) as file:
        msh.topology.create_connectivity(fdim, tdim)
        file.write_mesh(msh)
        file.write_meshtags(
            ct, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )
        file.write_meshtags(
            ft, msh.geometry, geometry_xpath=f"/Xdmf/Domain/Grid[@Name='{msh.name}']/Geometry"
        )

# Initialize Gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)

# Create model
model = gmsh.model()
option = gmsh.option()

# Create a Gmsh model of a 2D hollow circle
model = gmsh_hollow_circle(model, "HollowCircle", option)
model.setCurrent("HollowCircle")

# Create a DOLFINx mesh from the Gmsh model and save to an XDMF file
create_mesh(MPI.COMM_SELF, model, "hollow_circle", "out_gmsh/hollow_circle.xdmf", "w")

# We can use this to clear all the model data:
# gmsh.clear()
# gmsh.finalize()

# The generated mesh can be visualized using ParaView.
gdim = 2
gmsh_model_rank = 0
mesh_comm = MPI.COMM_WORLD

domain, cell_markers, facet_markers = gmshio.model_to_mesh(gmsh.model, mesh_comm, gmsh_model_rank, gdim=gdim)

# Material properties
E = 1e6 # Pa
nu = 0.33
G = E / (2*(1+nu))
lambda_ = nu * E /((1+nu)*(1-2*nu))
mu = G

x = ufl.SpatialCoordinate(domain)

# Functional spaces
V = fem.functionspace(domain, ("CG", 2, (domain.geometry.dim,)))

# Test and trial functions
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_markers)
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

p = fem.Constant(domain, default_scalar_type(100.0))

def on_boundary(x):
    return np.isclose(np.sqrt(x[0]**2 + x[1]**2), Re)

boundary_dofs = fem.locate_dofs_geometrical(V, on_boundary)

n = ufl.FacetNormal(domain)

def sigma(u, lambda_, mu):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*mu*epsilon(u)

def epsilon(u):
    return ufl.sym(ufl.grad(u))

a = ufl.inner(sigma(u, lambda_, mu), epsilon(v)) * x[0] * ufl.dx
L = ufl.inner(-p*n, v) * x[0] * ds

# Define the boundary condition
u_bc = np.array((0,) * domain.geometry.dim, dtype=default_scalar_type)
boundary_dofs = fem.locate_dofs_geometrical(V, on_boundary)
bc = fem.dirichletbc(u_bc, boundary_dofs, V)

problem = LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

uh = problem.solve()

import pyvista as pv
from dolfinx.plot import vtk_mesh
import numpy as np

# Extract mesh data from DOLFINx
topology, cell_types, geometry = vtk_mesh(V)
grid = pv.UnstructuredGrid(topology, cell_types, geometry)

# Check the length and dimensions of uh.x.array
num_points = geometry.shape[0]
if len(uh.x.array) != num_points * 2:
    raise ValueError(f"Length of uh.x.array ({len(uh.x.array)}) does not match expected number of components ({num_points * 2})")

# Reshape the array for 2D vectors
vector_data = uh.x.array.reshape((num_points, 2))

# Compute magnitudes of the vectors for warping
magnitude = np.sqrt(np.sum(vector_data**2, axis=1))

# Add vector magnitudes to the grid as scalar data
grid.point_data["magnitude"] = magnitude

# Warp the grid by scalar magnitude
warped = grid.warp_by_scalar("magnitude", factor=1.5)

# Create a PyVista Plotter object for standalone GUI window
plotter = pv.Plotter()
plotter.add_mesh(grid, style="wireframe", color="k")
plotter.add_mesh(warped, show_edges=True, scalars="magnitude", cmap="viridis")

# Show the plot in a GUI window
plotter.show()
