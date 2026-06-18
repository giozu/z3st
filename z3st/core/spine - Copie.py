# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.1.0 (2025)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import importlib
import dolfinx.nls.petsc

from petsc4py import PETSc
import dolfinx
import dolfinx.fem as fem
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
from z3st.models.plasticity_model import PlasticityModel
from z3st.models.contact_model import ContactModel


class Spine(
    Config, FiniteElementSetup, Solver, ThermalModel, MechanicalModel, GapModel, DamageModel, ClusterDynamicsModel, PlasticityModel
):
    """Main Z3ST simulation driver."""

    def __init__(self, input_file, mesh_file, geometry):

        self.current_step = 0
        # Initialize dictionary-based fields for contact/multi-domain
        self.subdomains = {}
        self.u_new = {}
        self.u_old = {}

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
        self.subdomain_areas = getattr(self.mgr, "subdomain_areas", {})
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
        if self.on.get("plasticity", False):
            PlasticityModel.__init__(self)

        # Initialize Contact Model if enabled
        self.contact_model = None
        if self.input_file.get("models", {}).get("contact", False):
            labels = self.geometry["labels"]
            # We expect 'steel_in' and 'steel_o' or similar for contact cases
            for key in ["steel_in", "steel_o"]:
                if key in labels:
                    self.subdomains[key] = labels[key]
            
            self.contact_model = ContactModel(
                mesh=self.mesh,
                subdomains=self.subdomains,
                V_m={name: self.V_m for name in self.subdomains}, # Shared V_m on same mesh
                facet_tags=self.facet_tags,
                geometry_config=self.geometry,
                contact_config=self.input_file.get("contact", {})
            )

    def parameters(self, lhr):
        self.g = 0.0  # m/s2
        self.lhr = lhr

    def resolve_function(self, path: str):
        module_path, func_name = path.rsplit(".", 1)
        mod = importlib.import_module(module_path)
        return getattr(mod, func_name)

    def load_materials(self, **materials):
        print(f"[spine.load_materials]")

        self.materials = {}
        
        lc = getattr(self, "dmg_cfg", {}).get("lc")

        for name, mat in materials.items():
            print(f"Material loaded: {name}")

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

            # YAML quirk: scientific notation without explicit +/- in the exponent
            # (e.g. `1.0e9`) is parsed as a *string*, not a float. Coerce defensively
            # so users don't hit a confusing TypeError downstream. A genuine symbolic
            # Gc string (e.g. "module.func") will fail the float() and keep its string
            # form for the resolver below.
            if isinstance(sigma_c, str):
                try:
                    sigma_c = float(sigma_c)
                    mat["sigma_c"] = sigma_c
                except ValueError:
                    pass  # leave as-is (no symbolic sigma_c path is currently used)
            if isinstance(Gc, str):
                try:
                    Gc = float(Gc)
                    mat["Gc"] = Gc
                except ValueError:
                    pass  # genuine "module.func" symbolic path — handled below

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

                elif Gc is not None and isinstance(Gc, (int, float, np.floating, np.integer)):
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

            if self.on.get("plasticity", False) and constitutive_mode == "lame" \
                    and "yield_strength" in mat:
                mat["constitutive_mode"] = "plasticity"
                constitutive_mode = "plasticity"
                print(f"  → constitutive model promoted to: plasticity (yield_strength present)")

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

        if self.on.get("thermal", False): self.set_thermal_boundary_conditions(self.V_t)
        if self.on.get("mechanical", False): self.set_mechanical_boundary_conditions(self.V_m)
        if self.on.get("damage", False): self.set_damage_boundary_conditions(self.V_d)

    def initialize_fields(self):
        print(f"[spine.initialize_fields]")

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
                T_init = mat.get("T_initial", mat["T_ref"])
                self.T.x.array[dofs] = T_init
                print(f"    Set {len(dofs)} DOFs to {T_init:.2f} K")

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

            # Initialize sub-domain fields for contact logic
            for name in self.subdomains:
                self.u_new[name] = dolfinx.fem.Function(self.V_m, name=f"u_{name}")
                self.u_old[name] = dolfinx.fem.Function(self.V_m, name=f"u_old_{name}")

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
                # Use subdomain specific area if available, fallback to global area
                area_val = self.subdomain_areas.get(name, self.area)
                q_val = self.lhr / area_val if area_val > 0 else 0.0
                # Accumulate: if multiple sources are configured on the same
                # material (e.g. fissile + gamma_heating below), they should
                # add — not overwrite.
                self.q_third.x.array[dofs] += q_val
                print(f"  q_third += {q_val:.3e} W/m³ (fissile: {mat.get('fissile', False)})")
                if self.perimeter > 0:
                    print(f"  Heat flux = {self.lhr / self.perimeter:.3e} W/m2")

            if float(mat.get("gamma_heating", 0.0)) > 0.0:
                # Cylindrical and spherical gamma-decay correlations use
                # `inner_radius` as the reference surface. If it is zero
                # we'd hit K_0(0) = +inf (cyl) or 1/r at r=0 (sphere); the
                # spurious result is silent zero or NaN heating. Surface up.
                if (
                    self.geometry_type in ("cyl", "cylinder", "sphere")
                    and float(getattr(self, "inner_radius", 0.0) or 0.0) == 0.0
                ):
                    raise ValueError(
                        f"Material '{name}' has gamma_heating > 0 with "
                        f"geometry_type='{self.geometry_type}' and inner_radius == 0. "
                        f"The decay correlation requires a non-zero inner radius "
                        f"as the reference surface; set inner_radius > 0 in geometry.yaml "
                        f"or use geometry_type='rect'."
                    )

                q_third_0 = float(mat["gamma_heating"])
                mu = float(mat["mu_gamma"])

                def f(x, q_third_0=q_third_0, mu=mu):
                    if self.geometry_type == "rect":
                        return q_third_0 * np.exp(-x[0] * mu)
                    elif self.geometry_type in ["cyl", "cylinder"]:
                        import scipy.special as sp

                        if (
                            self.regime == "axisymmetric"
                            or self.regime == "2d"
                        ):
                            # 2D axisymmetric case, x[0] = r, x[1] = z
                            # 2D axisymmetric case, x[0] = x, x[1] = y
                            radius = x[0]
                        elif self.regime == "3d":
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
                # Accumulate (see CODE-P1-10): allows fissile + gamma_heating
                # on the same material to combine rather than overwrite.
                self.q_third.x.array[dofs] += f_func.x.array[dofs]

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
                # Elastic energy uses the elastic strain (eps - alpha*(T - T_ref)*I)
                # so uniform thermal expansion does not appear as stored elastic energy.
                T_field = getattr(self, "T", None) if self.on.get("thermal", False) else None
                self.energy_density[name] = self.elastic_energy_density(self.u, mat, T=T_field)
                self.stress_mech[name] = self.sigma_mech(self.u, mat)
            
            if self.on.get("thermal", False):
                self.stress_th[name] = self.sigma_th(self.T, mat)      
                      
            if name in self.stress_mech and name in self.stress_th:
                self.stress[name] = self.stress_mech[name] + self.stress_th[name]
            elif name in self.stress_mech:
                self.stress[name] = self.stress_mech[name]
            elif name in self.stress_th:
                self.stress[name] = self.stress_th[name]

    def _mechanical_step(self, u_new_ref, u_old_ref, bcs_mech_ref, rtol_mech, stag_tol_mech, prev_res_u, T_current=None):
        """
        Etape mécanique adaptée pour le contact multi-domaine.
        Note: On utilise les dictionnaires self.u_new/u_old si contact_model est présent.
        """
        if self.contact_model is None:
            # Logique standard mono-domaine (fallback vers MechanicalModel)
            return super()._mechanical_step(u_new_ref, u_old_ref, bcs_mech_ref, rtol_mech, stag_tol_mech, prev_res_u, T_current) # type: ignore

        # --- Début logique Contact ---
        # u_old_ref holds the global displacement from the previous staggered iteration
        u_old_ref.x.array[:] = u_new_ref.x.array

        # Synchronize subdomain functions with global reference before assembly
        for name in self.u_new:
            self.u_old[name].x.array[:] = u_new_ref.x.array # Store for relaxation
            self.u_new[name].x.array[:] = u_new_ref.x.array # Use current global state for assembly

        print("\n[INFO] Assembling multi-domain mechanical problem (Contact)...")

        # Mise à jour des conditions aux limites dépendantes du temps (Dirichlet)
        for _, bc_list in self.dirichlet_mechanical.items():
            for bc in bc_list:
                if isinstance(bc, dict) and "raw" in bc:
                    idx = min(self.current_step, len(bc["raw"]) - 1)
                    bc["const"].value = np.array(bc["raw"][idx], dtype=dolfinx.default_scalar_type)

        # Mise à jour des tractions (Neumann)
        for domain, bc_list in self.traction.items():
            for bc in bc_list:
                if "raw" in bc:
                    idx = min(self.current_step, len(bc["raw"]) - 1)
                    bc["const"].value = np.array(bc["raw"][idx], dtype=dolfinx.default_scalar_type)

        # Assemblage des formes faibles pour chaque sous-domaine
        F_m = {name: 0 for name in self.subdomains}
        v_m = ufl.TestFunction(self.V_m)

        for name, tag in self.subdomains.items():
            dx = ufl.dx(subdomain_data=self.cell_tags, domain=self.mesh, subdomain_id=tag)
            material = self.materials[name]
            
            # Utiliser u_new_ref pour que la Jacobienne (dérivée de F par rapport à u) soit non nulle
            sigma = self.sigma_mech(u_new_ref, material)
            # Ajout force de volume (gravité) - Note: sigma_mech should use u_new_ref
            rho = dolfinx.default_scalar_type(material["rho"])
            body_force = fem.Constant(self.mesh, (0.0, -rho * 9.81))
            
            F_m[name] += ufl.inner(sigma, self.epsilon(v_m)) * dx - ufl.dot(body_force, v_m) * dx

            # Appliquer les tractions sur ce domaine
            if name in self.traction:
                for bc in self.traction[name]:
                    ds = ufl.ds(subdomain_data=self.facet_tags, domain=self.mesh, subdomain_id=bc["id"])
                    normal = self._get_normal_vector(bc["id"])
                    F_m[name] -= ufl.dot(ufl.as_vector(bc["const"] * normal), v_m) * ds

        # Ajouter le résidu de contact
        F_m = self.contact_model.add_contact_residual(F_m, u_new_ref)

        # Collecter toutes les conditions de Dirichlet
        bcs_mech = []
        for domain, bc_list in self.dirichlet_mechanical.items():
            for bc in bc_list:
                if isinstance(bc, dict):
                    bcs_mech.append(bc["value"]) # Extraire l'objet DirichletBC réel
                else:
                    bcs_mech.append(bc)

        # Résolution non-linéaire couplée (Somme des résidus)
        F_total = sum(F_m.values())

        # Calcul de la Jacobienne (forme tangente)
        # Indispensable pour que NonlinearProblem expose l'attribut '.a' attendu par NewtonSolver
        J_total = ufl.derivative(F_total, u_new_ref)
        
        petsc_opts = {
            "snes_type": "newtonls",
            "snes_atol": rtol_mech, # Use rtol_mech for SNES absolute tolerance
            "snes_rtol": rtol_mech, # Use rtol_mech for SNES relative tolerance
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        }

        # NonlinearProblem ne prend pas 'petsc_options' (le dict) dans son constructeur.
        # On passe tous les arguments en mots-clés pour satisfaire nanobind.
        problem = dolfinx.fem.petsc.NonlinearProblem(
            fem.form(F_total),
            u_new_ref,
            bcs=bcs_mech,
            J=fem.form(J_total),
            petsc_options_prefix="mechanical_"
        )
        
        # Injecter la mise à jour des contacts dans la boucle Newton
        def update_problem():
            self.contact_model.compute_contact_forces({"steel_in": u_new_ref, "steel_o": u_new_ref})
        problem.update = update_problem

        # Utiliser le solveur de Newton de DOLFINx (wrapper PETSc).
        solver = dolfinx.nls.petsc.NewtonSolver(self.mesh.comm, problem)
        solver.atol = rtol_mech * 100.0 # Relaxer la tolérance absolue pour le contact
        solver.rtol = rtol_mech
        solver.max_it = 500
        solver.report = True
        
        # Appliquer les options PETSc dans la base de données globale via le préfixe
        opts = PETSc.Options()
        for k, v in petsc_opts.items():
            opts.setValue(f"mechanical_{k}", str(v))
        
        # Configure the underlying linear solver (KSP) to use MUMPS
        ksp = solver.krylov_solver
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverType("mumps")

        solver.solve(u_new_ref)  # Résoudre pour le champ de déplacement global

        # Relaxation et mise à jour du champ principal pour Solver.solve_staggered
        norm_du = 0.0
        for name in self.u_new:
            self.u_new[name].x.array[:] = (
                self.relax_u * self.u_new[name].x.array +
                (1 - self.relax_u) * self.u_old[name].x.array
            )
            # Calcul de la norme de changement
            diff = self.u_new[name].x.array - self.u_old[name].x.array
            norm_du += np.linalg.norm(diff)**2
        
        # Apply relaxation to the global field u_new_ref
        u_new_ref.x.array[:] = (
            self.relax_u * u_new_ref.x.array +
            (1 - self.relax_u) * u_old_ref.x.array
        )
        dolfinx.fem.set_bc(u_new_ref.x.array, bcs_mech) # Apply BCs to the relaxed global field

        norm_du = np.sqrt(norm_du)
        rel_norm_du = norm_du # Simplifié

        # Relaxation adaptative
        if self.input_file["solver_settings"].get("relax_adaptive", False):
            if prev_res_u is not None:
                if norm_du < prev_res_u:
                    self.relax_u = min(self.relax_u * self.relax_growth, self.relax_max)
                else:
                    self.relax_u = max(self.relax_u * self.relax_shrink, self.relax_min)
            prev_res_u = norm_du

        return converged and (norm_du < stag_tol_mech), norm_du, rel_norm_du, prev_res_u

    def _get_normal_vector(self, region_id):
        """Récupère le vecteur normal basé sur les labels de géométrie."""
        labels = self.geometry["labels"]
        # Support specific contact labels from geometry.yaml
        if region_id in [labels.get("interface"), labels.get("steel_in_right"), labels.get("steel_o_left")]:
            return np.array([1.0, 0.0])
        elif region_id == labels.get("outer_radius"):
            return np.array([1.0, 0.0])
        elif region_id == labels.get("inner_radius"):
            return np.array([-1.0, 0.0])
        return np.array([0.0, 0.0])
