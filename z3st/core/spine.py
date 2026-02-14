# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import importlib

import dolfinx
import ufl
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
from z3st.models.cluster_dynamic_model import ClusterDynamicsModel


class Spine(
    Config, FiniteElementSetup, Solver, ThermalModel, MechanicalModel, GapModel, DamageModel, ClusterDynamicsModel
):
    """Main Z3ST simulation driver."""

    def __init__(self, input_file, mesh_file, geometry):

        self.current_step = 0

        mesh, cell_tags, facet_tags = load_mesh(mesh_file, comm=MPI.COMM_WORLD)

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

        # Initialize models
        if self.on.get("thermal", False):
            ThermalModel.__init__(self)
        if self.on.get("mechanical", False):
            MechanicalModel.__init__(self)
        if self.on.get("gap", False):
            GapModel.__init__(self)
        if self.on.get("damage", False):
            DamageModel.__init__(self)
        if self.on.get("cluster", False):
            ClusterDynamicsModel.__init__(self)

    def parameters(self, lhr):
        self.g = 0.0  # m/s2
        self.lhr = lhr

    def resolve_function(self, path: str):
        module_path, func_name = path.rsplit(".", 1)
        mod = importlib.import_module(module_path)
        return getattr(mod, func_name)

    def load_materials(self, **materials):

        print("\n")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print(f"--. spine - load_materials --..")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print("\n")

        self.materials = {}
        
        lc = getattr(self, "dmg_cfg", {}).get("lc")

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

            sigma_c = mat.get("sigma_c")
            Gc = mat.get("Gc")

            if "Gc" in mat:
                if isinstance(mat["Gc"], str):
                    print(f"  → Gc defined as symbolic function: {mat['Gc']}")
                    Gc_func = self.resolve_function(mat["Gc"])
                    mat["_Gc_func"] = Gc_func
                else:
                    print(f"  → Gc defined as constant: {mat['Gc']}")
            else:
                print(f"  → Gc not defined for {name}")

            dmg_type = getattr(self, "dmg_cfg", {}).get("type")

            if lc:
                if sigma_c is not None:
                    if dmg_type == "AT2":
                        Gc = (256/27) * lc * (sigma_c**2) / mat["E"]
                        print(f"  - Material '{name}': Gc (AT2) from sigma_c = {sigma_c:.2e} Pa")

                    elif dmg_type == "AT1":
                        Gc = (8/3) * lc * (sigma_c**2) / mat["E"]
                        print(f"  - Material '{name}': Gc (AT1) from sigma_c = {sigma_c:.2e} Pa")

                    mat["Gc"] = float(Gc)

                elif Gc is not None and type(Gc) == float:
                    if dmg_type == "AT2":
                        sigma_c = ((27 * mat["E"] * Gc) / (256 * lc))**0.5
                        print(f"  - Material '{name}': sigma_c (AT2) from Gc = {Gc:.2f} J/m2")

                    elif dmg_type == "AT1":
                        sigma_c = ((3 * mat["E"] * Gc) / (8 * lc))**0.5
                        print(f"  - Material '{name}': sigma_c (AT1) from Gc = {Gc:.2f} J/m2")

                    mat["sigma_c"] = float(sigma_c)

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
        print("\n")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print(f"--. spine - set_boundary_conditions --..")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print("\n")

        with open(self.input_file["boundary_conditions_path"], "r") as f:
            print(
                f"Loading boundary conditions from '{self.input_file['boundary_conditions_path']}'"
            )
            self.boundary_conditions = yaml.safe_load(f)

        if self.coupling == "staggered":
            if self.on.get("thermal", False): self.set_thermal_boundary_conditions(self.V_t)
            if self.on.get("mechanical", False): self.set_mechanical_boundary_conditions(self.V_m)
            if self.on.get("damage", False): self.set_damage_boundary_conditions(self.V_d)
        else:
            V_t_sub, V_t_map = self.W.sub(1).collapse()
            V_t_map = np.array(V_t_map, dtype=np.int32)
            V_u_sub, V_u_map = self.W.sub(0).collapse()
            V_u_map = np.array(V_u_map, dtype=np.int32)
            self.set_thermal_boundary_conditions(V_t_sub, V_t_map)
            self.set_mechanical_boundary_conditions(V_u_sub, V_u_map)

    def initialize_fields(self):
        print(f"[INITIALIZING FIELDS]")

        self.q_third = None
        self.T = None
        self.u = None
        self.D = None
        self.c = None

        # Temperature
        if self.on.get("thermal", False):
            self.q_third = dolfinx.fem.Function(self.V_t, name="q_third")
            self.q_third.x.array[:] = 0.0
            self.set_power()

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
        if self.on.get("mechanical", False):
            print("\nInitializing the displacement field...")
            self.u = dolfinx.fem.Function(self.V_m, name="Displacement")
            self.u.x.array[:] = 0.0
            self.u.x.scatter_forward()
            u_vals = self.u.x.array
            print(
                f"  Initial u: min={u_vals.min():.2e} m, max={u_vals.max():.2e} m, mean={u_vals.mean():.2e} m"
            )

        # Temperature/displacement (mixed)
        if self.on.get("mechanical", False) and self.on.get("thermal", False):
            if not hasattr(self, "sol_mixed"):
                print("Initializing self.sol_mixed")
                self.sol_mixed = dolfinx.fem.Function(self.W, name="MixedSolution")
                try:
                    print("Initialize with the current state (displacement and temperature fields)")
                    self.sol_mixed.sub(0).interpolate(self.u)
                    self.sol_mixed.sub(1).interpolate(self.T)
                except Exception:
                    self.sol_mixed.x.array[:] = 0.0

        # Damage variables
        if self.on.get("damage", False):
            print("\nInitializing the damage field...")
            self.D = dolfinx.fem.Function(self.V_d, name="Damage")
            self.D.x.array[:] = 0.0  # undamaged initial state
            self.H = dolfinx.fem.Function(self.Q, name="CrackDrivingForce")
            self.H.x.array[:] = 0.0

        # CD variables
        if self.on.get("cluster", False):
            print("\nInitializing the cluster density field...")
            self.c = dolfinx.fem.Function(self.V_c, name="ClusterDensity")
            self.c.x.array[:] = 0.0
            self.c_n = dolfinx.fem.Function(self.V_c, name="ClusterDensity_old")
            self.c_n.x.array[:] = 0.0

            print("\nSetting cluster initial conditions...")
            self.set_cluster_initial_conditions()

        # if self.u:
        #      print(f"  Displacement space (self.u):          {self.u.function_space.ufl_element()}")
        # if self.T:
        #      print(f"  Temperature space (self.T):           {self.T.function_space.ufl_element()}")

        # Material properties
        for name, mat in self.materials.items():
            if "_k_func" in mat and self.T:
                k_func = mat["_k_func"]
                mat["k"] = k_func(self.T)
                print("\nk expression for", name, "→", mat["k"])

            if "_Gc_func" in mat:
                Gc_func = mat["_Gc_func"]
                mat["Gc"] = Gc_func(self.mesh)
                print("\nGc expression for", name, "→", mat["Gc"])

                # Calculate sigma_c from Gc using UFL/dolfinx compatible operations
                lc = getattr(self, "dmg_cfg", {}).get("lc")
                dmg_type = getattr(self, "dmg_cfg", {}).get("type")
                
                if dmg_type == "AT2" and "E" in mat:
                    mat["sigma_c"] = ufl.sqrt((27 * mat["E"] * mat["Gc"]) / (256 * lc))
                    print(f"  - Material '{name}': sigma_c (AT2) evaluated from Gc expression")

                elif dmg_type == "AT1" and "E" in mat:
                    mat["sigma_c"] = ufl.sqrt((3 * mat["E"] * mat["Gc"]) / (8 * lc))
                    print(f"  - Material '{name}': sigma_c (AT1) evaluated from Gc expression")


    def set_power(self):
        if self.q_third is None:
            return

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

                        if (
                            self.mech_cfg["mechanical_regime"].lower() == "axisymmetric"
                            or self.mech_cfg["mechanical_regime"].lower() == "2d"
                        ):
                            # 2D axisymmetric case, x[0] = r, x[1] = z
                            # 2D axisymmetric case, x[0] = x, x[1] = y
                            radius = x[0]
                        elif self.mech_cfg["mechanical_regime"].lower() == "3d":
                            # 3D cartesian case, x[0] = x, x[1] = y
                            radius = np.sqrt(x[0] ** 2 + x[1] ** 2)

                        return q_third_0 * sp.k0(mu * radius) / sp.k0(mu * self.inner_radius)
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

    def solve(self, max_iters=100, dt=0.0):
        print("\n")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print(f"--. spine - solve --..")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print("\n")

        print(f"Current step = {self.current_step} | dt = {dt:.2e} s")
        print(f"Coupling = {self.coupling}")

        if self.coupling == "staggered":
            self.solve_staggered(
                max_iter=max_iters,
                dt=dt,
                rtol_th=self.th_cfg.get("rtol", 1e-6) if self.on.get("thermal") else 1e-6,
                rtol_mech=self.mech_cfg.get("rtol", 1e-6) if self.on.get("mechanical") else 1e-6,
                rtol_dmg=self.dmg_cfg.get("rtol", 1e-6) if hasattr(self, 'dmg_cfg') else 1e-6,
                stag_tol_th=self.th_cfg.get("stag_tol", 1e-4) if self.on.get("thermal") else 1e-4,
                stag_tol_mech=self.mech_cfg.get("stag_tol", 1e-4) if self.on.get("mechanical") else 1e-4,
                stag_tol_dmg=self.dmg_cfg.get("stag_tol", 1e-4) if hasattr(self, 'dmg_cfg') else 1e-4
            )
        else:
            raise ValueError(f"Unknown coupling strategy: {self.coupling}. Only staggered coupling is supported.")

    def get_results(self):
        if not (self.on.get("mechanical", False) or self.on.get("thermal", False)):
            return

        print("Computing symbolic result fields (strain, stress, ...)")
        self.stress = {}
        self.stress_mech = {}
        self.stress_th = {}
        self.energy_density = {}
        
        if self.on.get("mechanical", False):
            self.strain = self.epsilon(self.u)
        else:
            self.strain = None

        for name, mat in self.materials.items():
            if self.on.get("mechanical", False):
                self.energy_density[name] = self.elastic_energy_density(self.u, mat)
                self.stress_mech[name] = self.sigma_mech(self.u, mat)
            
            if self.on.get("thermal", False):
                self.stress_th[name] = self.sigma_th(self.T, mat)      
                      
            if name in self.stress_mech and name in self.stress_th:
                self.stress[name] = self.stress_mech[name] + self.stress_th[name]
            elif name in self.stress_mech:
                self.stress[name] = self.stress_mech[name]
            elif name in self.stress_th:
                self.stress[name] = self.stress_th[name]

