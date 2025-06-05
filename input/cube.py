import dolfinx
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from finiteElementSettings import FiniteElementSetting
from models.thermal_model_cube import ThermalSolver
from models.mechanical_model_cube import MechanicalSolver
import pyvista

class Cube(FiniteElementSetting):
    def volume(self, edge, log_path="output/simulation_log.txt"):
        self.edge = edge
        self.area = edge**2
        self.height = edge

        dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.tags)
        mesh_volume = dolfinx.fem.assemble_scalar(dolfinx.fem.form(1 * dx))
        analytical_volume = self.area * self.height
        rel_error = abs(mesh_volume - analytical_volume) / analytical_volume * 100

        with open(log_path, "a") as f:
            f.write("\nVolume Comparison\n")
            f.write(f"Mesh volume: {mesh_volume:.6e} m3\n")
            f.write(f"Analytical volume: {analytical_volume:.6e} m3\n")
            f.write(f"Relative error: {rel_error:.3e} %\n")

    def parameters(self, material, g, lhr, T_ref):
        self.g = g
        self.lhr = lhr
        self.T_ref = T_ref
        self.material = material
        self.material["lmbda"] = material["E"] * material["nu"] / ((1 + material["nu"]) * (1 - 2 * material["nu"]))
        self.material["mu"] = material["E"] / (2 * (1 + material["nu"]))
        self.material["bulk_modulus"] = material["E"] / (3 * (1 - 2 * material["nu"]))
        self.q_third = lhr / self.area

    def boundary_conditions(self, pressure, temperature):
        self.bottom_label = 1
        self.top_label = 6

        self.ds = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.ft)
        self.dofs_bottom = dolfinx.fem.locate_dofs_topological(self.V_m, self.fdim, self.ft.find(self.bottom_label))

        self.bc_M = [
            dolfinx.fem.dirichletbc(dolfinx.fem.Constant(self.mesh, (0.0, 0.0, 0.0)), self.dofs_bottom, self.V_m)
        ]

        self.T_bottom = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(temperature))
        self.bc_T = [
            dolfinx.fem.dirichletbc(self.T_bottom, self.dofs_bottom, self.V_t)
        ]

        self.normal = ufl.FacetNormal(self.mesh)
        self.traction = dolfinx.fem.Constant(self.mesh, pressure) * self.normal

        topology, cell_types, _ = dolfinx.plot.vtk_mesh(self.mesh, self.mesh.topology.dim - 1)
        surface = pyvista.UnstructuredGrid(topology, cell_types, self.mesh.geometry.x)

        plotter = pyvista.Plotter()
        unique_labels = np.unique(self.ft.values)

        label_map = {
                1: ("Bottom", "red"),
                2: ("Face 1", "blue"),
                3: ("Face 2", "green"),
                4: ("Face 3", "orange"),
                5: ("Face 4", "purple"),
                6: ("Top", "black")
            }

        print("Face labels present:", unique_labels)
        for label in unique_labels:
            facets = self.ft.find(label)
            print(f"Label {label}: {len(facets)} facets")

        for label in np.unique(self.ft.values):
            facets = self.ft.find(label)
            surface_cells = surface.extract_cells(facets)
            name, color = label_map.get(label, (f"Face {label}", "white"))
            plotter.add_mesh(surface_cells, show_edges=True, opacity=1.0, color=color, label=name)

        plotter.add_legend()
        plotter.show()

    def update_state(self):
        self.q_third = self.lhr / self.area
        self.thermal_solver.q_third = self.q_third

    def initialize_fields(self):
        self.T = dolfinx.fem.Function(self.V_t, name="Temperature")
        self.T.x.array[:] = float(self.T_bottom.value)

        self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
        self.u.x.array[:] = 0.0

        self.thermal_solver = ThermalSolver(
            mesh=self.mesh,
            V_t=self.V_t,
            q_third=self.q_third,
            k=self.material["k"],
            bc=self.bc_T,
            source_domain=ufl.dx
        )

        self.mechanical_solver = MechanicalSolver(
            mesh=self.mesh,
            V_m=self.V_m,
            bc=self.bc_M,
            lmbda=self.material["lmbda"],
            mu=self.material["mu"],
            alpha=self.material["alpha"],
            T=self.T,
            T_ref=self.T_ref,
            traction=self.traction,
            rho=self.material["rho"],
            g=self.g,
            ds=self.ds
        )

    def initialize_geometry(self, atol=1e-6):
        """
        Store original node indices for selected surfaces of the cube.
        """
        coords = self.mesh.geometry.x
        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
        edge = self.edge

        self.idx_bottom = np.where(np.isclose(z, 0.0, atol=atol))[0]
        self.idx_top = np.where(np.isclose(z, edge, atol=atol))[0]
        self.idx_left = np.where(np.isclose(x, 0.0, atol=atol))[0]
        self.idx_right = np.where(np.isclose(x, edge, atol=atol))[0]
        self.idx_front = np.where(np.isclose(y, 0.0, atol=atol))[0]
        self.idx_back = np.where(np.isclose(y, edge, atol=atol))[0]

    def solve(self):
        self.thermal_solver.solve(self.T)
        self.mechanical_solver.solve(self.u)

    def get_results(self):

        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u, lmbda, mu):
            return lmbda * ufl.tr(epsilon(u)) * ufl.Identity(3) + 2 * mu * epsilon(u)
        
        self.strain = epsilon(self.u)

        self.stress = (
            sigma(self.u, self.material["lmbda"], self.material["mu"]) -
            (3 * self.material["lmbda"] + 2 * self.material["mu"]) * self.material["alpha"] *
            (self.T - self.T_ref) * ufl.Identity(3)
        )
        return self.u, self.strain, self.stress, self.T
