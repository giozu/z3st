import numpy as np

import ufl
from ufl import ds, dx, grad, inner

from dolfinx import fem, io, mesh, plot
from dolfinx.io import gmshio # to import gmsh functions
from dolfinx import geometry # to define the line plot
from dolfinx import default_scalar_type
import dolfinx.fem.petsc

from mpi4py import MPI

from petsc4py.PETSc import ScalarType

import matplotlib.pyplot as plt

class linear_elasticty():
    def __init__(self):
        # Domain
        # self.domain = mesh.create_box(MPI.COMM_WORLD, [np.array([0,0,0]), np.array([1, 0.2, 0.2])], [20, 6, 6], cell_type=mesh.CellType.hexahedron)
        self.domain, self.tags, self.ft = gmshio.read_from_msh("cylindrical_shell.msh", MPI.COMM_WORLD, gdim = 3)

        self.tdim = self.domain.topology.dim
        self.fdim = self.tdim - 1
        self.domain.topology.create_connectivity(self.fdim, self.tdim)

        # Functional spaces
        self.V = fem.functionspace(self.domain, ("Lagrange", 1, (self.domain.geometry.dim,)))

        # Test, trial functions and measure
        self.ds = ufl.Measure("ds", domain=self.domain)
        self.u = ufl.TrialFunction(self.V)
        self.v = ufl.TestFunction(self.V)

    def parameters(self, rho, g, lmbda, mu):
        self.rho = rho
        self.g = g
        self.lmbda = lmbda
        self.mu = mu

    def boundaryConditions(self):
        # Dirichlet boundary conditions
        # self.facets = mesh.locate_entities_boundary(self.domain, self.fdim, boundary)
        # self.dofs_D_top = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(13))
        # self.dofs_D_bottom = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(14))
        # self.dofs_D_outer = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(15))
        # self.dofs_D_inner = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(16))

        # Plain strain on top and bottom
        # self.bc_top = fem.dirichletbc(np.array([0, 0, 0], dtype=default_scalar_type), dofs=self.dofs_D_top, V=self.V)
        # self.bc_bottom = fem.dirichletbc(np.array([0, 0, 0], dtype=default_scalar_type), dofs=self.dofs_D_bottom, V=self.V)
        # self.bc = [self.bc_top, self.bc_bottom]

        # self.dofs_D = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(11)) # bottom
        self.dofs_D = fem.locate_dofs_topological(V=self.V, entity_dim=self.fdim, entities=self.ft.find(12)) # bottom

        self.bc = [fem.dirichletbc(np.array([0, 0, 0], dtype=default_scalar_type), dofs=self.dofs_D, V=self.V)]

        # Neumann boundary conditions
        self.traction = fem.Constant(self.domain, ScalarType((0, 0, 0)))

    def assemble(self):
        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*mu*epsilon(u)

        def epsilon(u):
            return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

        self.f = fem.Constant(self.domain, ScalarType((0.0, 0.0, - self.rho * self.g)))
        self.a = ufl.inner(sigma(self.u, self.lmbda, self.mu), epsilon(self.v)) * ufl.dx
        self.L = ufl.dot(self.f, self.v) * ufl.dx + ufl.dot(self.traction, self.v) * self.ds

    def solve(self):
        self.problem = fem.petsc.LinearProblem(self.a, self.L, bcs=self.bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        self.displacement = self.problem.solve()

        def sigma(u, lmbda, mu):
            return lmbda * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*mu*epsilon(u)

        def epsilon(u):
            return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
        
        self.stressField = sigma(self.displacement, self.lmbda, self.mu)
        self.strainField = epsilon(self.displacement)                

# Scaled variable
E = 1e6 # Pa
nu = 0.33

G = E / (2*(1+nu))

lame1 = nu * E /((1+nu)*(1-2*nu))
lame2 = G

print(f"Young's modulus = {E*1e-9:.2f} GPa")
print(f"Poisson's ratio = ", nu)
print(f"Shear modulus = {G*1e-9:.2f} GPa")
print(f"Lamé first parameter = {lame1*1e-9:.2f} GPa")
print(f"Lamé second parameter = {lame2*1e-9:.2f} GPa")

rho = 6560 # kg/m3
g = 9.81

problem = linear_elasticty()

# Volume of the mesh = volume integral
dx_subdomain = ufl.Measure("dx", domain=problem.domain)
volume = fem.assemble_scalar(fem.form(1 * dx_subdomain))
print("---------")
print(f"Volume = ", volume)

problem.parameters(rho, g / volume, lame1, lame2)
problem.boundaryConditions()
problem.assemble()
problem.solve()

uh = problem.displacement
