import numpy as np

import ufl
from ufl import ds, dx, grad, inner

from dolfinx import fem, io, mesh, plot
from dolfinx.io import gmshio # to import gmsh functions
from dolfinx import geometry # to define the line plot
from dolfinx import default_scalar_type
import dolfinx.fem.petsc

from mpi4py import MPI

# from petsc4py import PETSc
from petsc4py.PETSc import ScalarType

import matplotlib.pyplot as plt

class linear_elasticty():
    def __init__(self):
        # Domain
        # self.domain = mesh.create_box(MPI.COMM_WORLD, [np.array([0,0,0]), np.array([1, 0.2, 0.2])], [20, 6, 6], cell_type=mesh.CellType.hexahedron)
        self.domain, self.tags, self.ft = gmshio.read_from_msh("pipe.msh", MPI.COMM_WORLD, gdim = 3)
        self.tdim = self.domain.topology.dim
        self.fdim = self.tdim - 1
        self.domain.topology.create_connectivity(self.fdim, self.tdim)

        # Functional spaces
        self.V = fem.functionspace(self.domain, ("Lagrange", 1, (self.domain.geometry.dim,)))

        # Test and trial functions
        self.ds = ufl.Measure("ds", domain=self.domain)
        self.u = ufl.TrialFunction(self.V)
        self.v = ufl.TestFunction(self.V)

    def parameters(self, rho, g, lambda_, mu):
        self.rho = rho
        self.g = g
        self.lambda_ = lambda_
        self.mu = mu

    def boundaryConditions(self, bc):
        # Dirichlet boundary conditions
        # self.facets = mesh.locate_entities_boundary(self.domain, self.fdim, boundary)
        self.dofs_D = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(11))
        self.bc = fem.dirichletbc(value=bc, dofs=self.dofs_D, V=self.V)
        
        # Neumann boundary conditions
        self.traction = fem.Constant(self.domain, ScalarType((0, 0, 0)))

    def assemble(self):
        def sigma(u, lambda_, mu):
            return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*mu*epsilon(u)

        def epsilon(u):
            return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

        self.f = fem.Constant(self.domain, ScalarType((0.0, 0.0, - self.rho * self.g)))
        self.a = ufl.inner(sigma(self.u, self.lambda_, self.mu), epsilon(self.v)) * ufl.dx
        self.L = ufl.dot(self.f, self.v) * ufl.dx + ufl.dot(self.traction, self.v) * self.ds

    def solve(self):
        self.problem = fem.petsc.LinearProblem(self.a, self.L, bcs=[self.bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        self.solution = self.problem.solve()

        def sigma(u, lambda_, mu):
            return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*mu*epsilon(u)

        def epsilon(u):
            return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
        
        self.stressField = sigma(self.solution, self.lambda_, self.mu)
        self.strainField = epsilon(self.solution)                

        return self.solution
    

# Scaled variable
E = 100e9 # Pa
nu = 0.33

G = E / (2*(1+nu))

lame1 = nu * E /((1+nu)*(1-2*nu))
lame2 = G

print(f"Young's modulus = {E*1e-9:.2f} GPa")
print(f"Poisson's ratio = ", nu)
print(f"Shear modulus = {G*1e-9:.2f} GPa")
print(f"Lamé first parameter = {lame1*1e-9:.2f} GPa")
print(f"Lamé second parameter = {lame2*1e-9:.2f} GPa")

u_D = np.array([0, 0, 0], dtype=default_scalar_type)

rho = 6560 # kg/m3
g = 9.81

problem = linear_elasticty()

# Volume of the mesh = volume integral
dx_subdomain = ufl.Measure("dx", domain=problem.domain)
volume = fem.assemble_scalar(fem.form(1 * dx_subdomain))
print("---------")
print(f"Volume = ", volume)

problem.parameters(rho, g / volume, lame1, lame2)
problem.boundaryConditions(u_D)
problem.assemble()

uh = problem.solve()

import pyvista
pyvista.start_xvfb()

# Create plotter and pyvista grid
p = pyvista.Plotter()
topology, cell_types, geometry = plot.vtk_mesh(problem.V)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

# Attach vector values to grid and warp grid by vector
grid["u"] = uh.x.array.reshape((geometry.shape[0], 3))
actor_0 = p.add_mesh(grid, style="wireframe", color="k")
warped = grid.warp_by_vector("u", factor=1.5)
actor_1 = p.add_mesh(warped, show_edges=True)
p.show_axes()
if not pyvista.OFF_SCREEN:
   p.show()
else:
   figure_as_array = p.screenshot("deflection.png")

# We could also use Paraview for visualizing this. We save the solution with XDMFFile.
# After opening the file deformation.xdmf in Paraview and pressing Apply, one can press the Warp by vector button Warp
# by vector or go through the top menu (Filters->Alphabetical->Warp by Vector) and press Apply.
# We can also change the color of the deformed beam by changing the value in the color menu color from Solid Color to Deformation.
with io.XDMFFile(problem.domain.comm, "deformation.xdmf", "w") as xdmf:
    xdmf.write_mesh(problem.domain)
    uh.name = "Deformation"
    xdmf.write_function(uh)

deviatoric_stress = problem.stressField -1./3*ufl.tr(problem.stressField)*ufl.Identity(len(problem.solution))
von_Mises = ufl.sqrt(3./2*ufl.inner(deviatoric_stress, deviatoric_stress))

V_von_mises = fem.functionspace(problem.domain, ("DG", 0))
stress_expr = fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
stresses = fem.Function(V_von_mises)
stresses.interpolate(stress_expr)

warped.cell_data["VonMises"] = stresses.vector.array
warped.set_active_scalars("VonMises")
p = pyvista.Plotter()
p.add_mesh(warped)
p.show_axes()
if not pyvista.OFF_SCREEN:
   p.show()
else:
   stress_figure = p.screenshot(f"stresses.png")

