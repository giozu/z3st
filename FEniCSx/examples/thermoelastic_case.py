import gmsh
from mpi4py import MPI
import sys
from dolfinx.io.gmshio import model_to_mesh
import ufl
from dolfinx.fem.petsc import LinearProblem
import pyvista as pv

# Geometry and Mesh
gmsh.initialize(sys.argv)
gmsh.model.add("boolean")
model_rank = 0
mesh_comm = MPI.COMM_WORLD
gdim = 2

L = 100.
l = 25.
R = l/3./2.
plate = gmsh.model.occ.addRectangle(0, 0, 0, L, l)
hole = gmsh.model.occ.addDisk(L/2, l/2, 0, R, R)
holed_plate = gmsh.model.occ.cut([(2, plate)],[(2, hole)])
gmsh.model.occ.synchronize()

surf = gmsh.model.getEntities(2)
gmsh.model.addPhysicalGroup(surf[0][0], [surf[0][1]], 3)

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 3)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3)
gmsh.model.mesh.generate(gdim)

if(False):
    gmsh.fltk.run()
domain, cell_markers, facet_markers = model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=gdim)
gmsh.finalize()

from dolfinx import fem, plot, default_scalar_type, mesh
import numpy as np

V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))
#V = (domain, ("CG", 1, (domain.geometry.dim, )))

def clamped_boundary(x):
    return np.isclose(x[0], 0)
def imposed_disp_boundary(x):
    return np.isclose(x[0], L)
fdim = domain.topology.dim - 1
boundary_facets_I = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
boundary_facets_C = mesh.locate_entities_boundary(domain, fdim, imposed_disp_boundary)
u_C = np.array([0, 0], dtype=default_scalar_type)
u_I = np.array([0, 0], dtype=default_scalar_type)
bc_I = fem.dirichletbc(u_C, fem.locate_dofs_topological(V, fdim, boundary_facets_I), V)
bc_C = fem.dirichletbc(u_I, fem.locate_dofs_topological(V, fdim, boundary_facets_C), V)

E_young = 1000
nu = 0.3

lambda_lame = E_young * nu / ((1 + nu) * (1 - 2 *nu))
mu_lame = E_young /(2 * (1 + nu))

alpha = 1e-5
Temp = 100.

T = fem.Constant(domain, default_scalar_type((0, 0)))
ds = ufl.Measure("ds", domain=domain)

def epsilon(u):
    return ufl.sym(ufl.grad(u))  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

def sigma(v, dT):
    eps_th = alpha * dT
    eps_el = epsilon(v) - eps_th * ufl.Identity(2)
    return lambda_lame * ufl.tr(eps_el) * ufl.Identity(2) + 2.0 * mu_lame * eps_el

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, default_scalar_type((0, 0)))
a = ufl.inner(sigma(u, Temp), epsilon(v)) * ufl.dx
Lin = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

F= ufl.inner(sigma(u, Temp), epsilon(v)) * ufl.dx- ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds
a, Lin  = ufl.system(F)

problem = LinearProblem(a, Lin, bcs=[bc_C, bc_I], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
Uf = problem.solve()