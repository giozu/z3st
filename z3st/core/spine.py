# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import importlib

import dolfinx
import numpy as np
import yaml
from mpi4py import MPI

from z3st.core.config import Config
from z3st.core.finite_element_setup import FiniteElementSetup
from z3st.core.mesh import load_mesh
from z3st.core.mesh.manager import MeshManager
from z3st.core.solver import Solver
from z3st.models.damage_model import DamageModel
from z3st.models.gap_model import GapModel
from z3st.models.mechanical_model import MechanicalModel
from z3st.models.thermal_model import ThermalModel


class Spine(
    Config, FiniteElementSetup, Solver, ThermalModel, MechanicalModel, GapModel, DamageModel
):
    """Main Z3ST simulation driver."""

    def __init__(self, input_file, mesh_file, geometry, gdim=3):

        mesh, cell_tags, facet_tags = load_mesh(mesh_file, comm=MPI.COMM_WORLD, gdim=gdim)

        self.mgr = MeshManager(mesh, cell_tags, facet_tags, geometry=geometry)
        self.mgr.summary()

        # Store references
        self.mesh = self.mgr.mesh
        self.cell_tags = self.mgr.cell_tags
        self.facet_tags = self.mgr.facet_tags
        self.label_map = self.mgr.label_map
        self.geometry = self.mgr.geometry
        self.geometry_type = self.mgr.geometry_type
        self.area = getattr(self.mgr, "area", 0.0)
        self.perimeter = getattr(self.mgr, "perimeter", 0.0)
        self.inner_radius = getattr(self.mgr, "inner_radius", 0.0)
        self.normal = self.mgr.normal
        self.tdim = self.mgr.tdim
        self.fdim = self.mgr.fdim

        # Initialize modules
        Config.__init__(self, input_file)
        FiniteElementSetup.__init__(self)
        Solver.__init__(self)
        ThermalModel.__init__(self)
        MechanicalModel.__init__(self)
        GapModel.__init__(self)
        if self.on.get("damage", False):
            DamageModel.__init__(self)

    def parameters(self, lhr):
        self.g = 0.0  # m/s2
        self.lhr = lhr

    def resolve_function(self, path: str):
        module_path, func_name = path.rsplit(".", 1)
        mod = importlib.import_module(module_path)
        return getattr(mod, func_name)

    def load_materials(self, **materials):
        print("[LOADING MATERIALS]")
        self.materials = {}

        for name, mat in materials.items():
            print(f"{name.capitalize()} material loaded:")

            if "E" in mat and "nu" in mat:
                mat["E"] = float(mat["E"])
                mat["nu"] = float(mat["nu"])
                mat["lmbda"] = mat["E"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
                mat["G"] = mat["E"] / (2 * (1 + mat["nu"]))
                mat["bulk_modulus"] = mat["E"] / (3 * (1 - 2 * mat["nu"]))
            else:
                print(
                    f"  [INFO] '{name}' has no elasticity parameters — skipping mechanical properties."
                )

            if "k" in mat:
                if isinstance(mat["k"], str):
                    print(f"  → k defined as symbolic function: {mat['k']}")
                    k_func = self.resolve_function(mat["k"])
                    mat["_k_func"] = k_func
                else:
                    print(f"  → k defined as constant: {mat['k']}")
            else:
                print(f"  → k not defined for {name}")

            constitutive_mode = mat.get("constitutive", "lame").lower()
            mat["constitutive_mode"] = constitutive_mode
            print(f"  → constitutive model: {constitutive_mode}")
            if constitutive_mode == "voigt" and mat.get("C_matrix") is not None:
                print("    using user-provided C_matrix (6x6)")

            self.materials[name] = mat
            for k in sorted(mat.keys()):
                v = mat[k]
                print(f"  {k:<15} → {v} ({type(v).__name__})")

    def set_boundary_conditions(self):
        print(f"[SETTING BOUNDARY CONDITIONS]")

        with open(self.input_file["boundary_conditions_path"], "r") as f:
            print(
                f"Loading boundary conditions from '{self.input_file['boundary_conditions_path']}'"
            )
            self.boundary_conditions = yaml.safe_load(f)

        if self.coupling == "staggered":
            self.set_thermal_boundary_conditions(self.V_t)
            self.set_mechanical_boundary_conditions(self.V_m)
        else:
            V_t_sub, V_t_map = self.W.sub(1).collapse()
            V_t_map = np.array(V_t_map, dtype=np.int32)
            V_u_sub, V_u_map = self.W.sub(0).collapse()
            V_u_map = np.array(V_u_map, dtype=np.int32)
            self.set_thermal_boundary_conditions(V_t_sub, V_t_map)
            self.set_mechanical_boundary_conditions(V_u_sub, V_u_map)

    def initialize_fields(self):
        print(f"[INITIALIZING FIELDS]")

        self.q_third = dolfinx.fem.Function(self.V_t, name="q_third")
        self.q_third.x.array[:] = 0.0
        self.set_power()

        # Temperature
        print("\nInitializing the temperature field...")
        self.T = dolfinx.fem.Function(self.V_t, name="Temperature")
        for name, mat in self.materials.items():
            print(f"  → Setting initial temperature for material: '{name}'")
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_t)
            self.T.x.array[dofs] = mat["T_ref"]
            print(f"    Set {len(dofs)} DOFs to {mat['T_ref']:.2f} K")

        self.T.x.scatter_forward()
        T_vals = self.T.x.array
        print(
            f"  Initial T: min={T_vals.min():.2f} K, max={T_vals.max():.2f} K, mean={T_vals.mean():.2f} K"
        )

        # Displacement
        print("\nInitializing the displacement field...")
        self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
        self.u.x.array[:] = 0.0
        self.u.x.scatter_forward()
        u_vals = self.u.x.array
        print(
            f"  Initial u: min={u_vals.min():.2e} m, max={u_vals.max():.2e} m, mean={u_vals.mean():.2e} m"
        )

        # Temperature/displacement (mixed)
        if not hasattr(self, "sol_mixed"):
            print("Initializing self.sol_mixed")
            self.sol_mixed = dolfinx.fem.Function(self.W, name="MixedSolution")
            try:
                print("Initialize with the current state (displacement and temperature fields)")
                self.sol_mixed.sub(0).interpolate(self.u)
                self.sol_mixed.sub(1).interpolate(self.T)
            except Exception:
                self.sol_mixed.x.array[:] = 0.0

        # Damage variables:
        if self.on.get("damage", False):
            print("\nInitializing the damage field...")
            self.D = dolfinx.fem.Function(self.V_d, name="Damage")
            self.D.x.array[:] = 0.0  # undamaged initial state
            self.H = dolfinx.fem.Function(self.Q, name="CrackDrivingForce")
            self.H.x.array[:] = 0.0

        print("\nComparing solution function spaces:")
        print(f"  Displacement space (self.u):          {self.u.function_space.ufl_element()}")
        print(f"  Mixed-space displacement (W.sub(0)):  {self.W.sub(0).ufl_element()}")
        print(f"  Temperature space (self.T):           {self.T.function_space.ufl_element()}")
        print(f"  Mixed-space temperature (W.sub(1)):   {self.W.sub(1).ufl_element()}")

        for name, mat in self.materials.items():
            if "_k_func" in mat:
                k_func = mat["_k_func"]
                mat["k"] = k_func(self.T)
                print("\nk expression for", name, "→", mat["k"])

    def set_power(self):
        print(f"[UPDATING q_third]")
        self.q_third.x.array[:] = 0.0

        for name, mat in self.materials.items():
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_t)

            if mat.get("fissile", False):
                print("Fissile material")
                q_val = self.lhr / self.area
                self.q_third.x.array[dofs] = q_val
                print(f"  q_third = {q_val:.3e} W/m³ (fissile: {mat.get('fissile', False)})")
                print(f"  Heat flux = {self.lhr / self.perimeter:.3e} W/m2")

            if float(mat["gamma_heating"]) > 0.0:
                q_third_0 = float(mat["gamma_heating"])
                mu = float(mat["mu_gamma"])

                def f(x, q_third_0=q_third_0, mu=mu):
                    if self.geometry_type == "rect":
                        return q_third_0 * np.exp(-x[0] * mu)
                    elif self.geometry_type in ["cyl", "cylinder"]:
                        import scipy.special as sp

                        return (
                            q_third_0
                            * sp.k0(mu * np.sqrt(x[0] ** 2 + x[1] ** 2))
                            / sp.k0(mu * self.inner_radius)
                        )
                    elif self.geometry_type == "sphere":
                        r = np.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)
                        return (
                            q_third_0
                            * (self.inner_radius / r)
                            * np.exp(-mu * (r - self.inner_radius))
                        )

                f_func = dolfinx.fem.Function(self.V_t)
                f_func.interpolate(f)
                self.q_third.x.array[dofs] = f_func.x.array[dofs]

        self.q_third.x.scatter_forward()

    def solve(self, max_iters=100):
        print("\nSolver settings:")
        print(f"  → Coupling : {self.coupling}")

        if self.coupling == "monolithic":
            self.solve_monolithic(tol=1e-6, max_iter=max_iters)
        elif self.coupling == "staggered":
            self.solve_staggered(
                max_iter=max_iters,
                rtol_th=self.th_cfg["rtol"],
                rtol_mech=self.rtol_mech,
                stag_tol_th=self.th_cfg["stag_tol"],
                stag_tol_mech=self.stag_tol_mech,
            )
        else:
            raise ValueError(f"Unknown coupling strategy: {self.coupling}")

    def get_results(self):
        print("Computing symbolic result fields (strain, stress, ...)")
        self.stress = {}
        self.stress_mech = {}
        self.stress_th = {}
        self.energy_density = {}
        self.strain = self.epsilon(self.u)

        # damage fields
        if self.on.get("damage", False):
            self.damage_field = {}

        for name, mat in self.materials.items():
            self.energy_density[name] = self.elastic_energy_density(self.u, mat)
            self.stress_mech[name] = self.sigma_mech(self.u, mat)
            self.stress_th[name] = self.sigma_th(self.T, mat)
            self.stress[name] = self.stress_mech[name] + self.stress_th[name]

            if self.on.get("damage", False):
                self.damage_field[name] = self.D

        return True
