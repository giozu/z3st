# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import importlib
import sys

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
from z3st.models.contact_model import ContactModel
from z3st.models.cracking_model import CrackingModel
from z3st.models.creep_model import CreepModel
from z3st.models.damage_model import DamageModel
from z3st.models.gap_model import GapModel
from z3st.models.mechanical_model import MechanicalModel
from z3st.models.thermal_model import ThermalModel
from z3st.models.cluster_dynamic_model import ClusterDynamicsModel
from z3st.models.plasticity_model import PlasticityModel
from z3st.models.porosity_migration_model import PorosityMigrationModel


class Spine(
    Config, FiniteElementSetup, Solver, ThermalModel, MechanicalModel, GapModel, ContactModel, DamageModel, ClusterDynamicsModel, PlasticityModel, CreepModel, CrackingModel, PorosityMigrationModel
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
        if self.on.get("contact", False):
            ContactModel.__init__(self)
        if self.on.get("damage", False):
            DamageModel.__init__(self)
        if self.on.get("cluster", False):
            ClusterDynamicsModel.__init__(self)
        if self.on.get("plasticity", False):
            PlasticityModel.__init__(self)
        if self.on.get("porosity", False):
            PorosityMigrationModel.__init__(self)

    def parameters(self, lhr):
        self.g = 0.0  # m/s2
        self.lhr = lhr
        # Rescale the elastic constants
        if getattr(self, "materials", None):
            self.update_cracking()

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
                # E and/or nu may be a symbolic "module.func" card; the derived
                # elastic constants are then built as UFL expressions in the T
                # field by initialize_fields (T does not exist yet here).
                # NOTE: symbolic E/nu is not yet compatible with the per-step
                # cracking rescale or the numpy creep predictor.
                def _is_symbolic(v):
                    if not isinstance(v, str):
                        return False
                    try:
                        float(v)
                        return False
                    except ValueError:
                        return True
                E_sym = _is_symbolic(mat["E"])
                nu_sym = _is_symbolic(mat["nu"])
                if E_sym:
                    print(f"  → E defined as symbolic function: {mat['E']}")
                    mat["_E_func"] = self.resolve_function(mat["E"])
                else:
                    mat["E"] = float(mat["E"])
                if nu_sym:
                    print(f"  → nu defined as symbolic function: {mat['nu']}")
                    mat["_nu_func"] = self.resolve_function(mat["nu"])
                else:
                    mat["nu"] = float(mat["nu"])
                if E_sym or nu_sym:
                    print(f"  → '{name}': elastic constants deferred to UFL(T)")
                else:
                    mat["lmbda"] = mat["E"] * mat["nu"] / ((1 + mat["nu"]) * (1 - 2 * mat["nu"]))
                    mat["G"] = mat["E"] / (2 * (1 + mat["nu"]))
                    mat["bulk_modulus"] = mat["E"] / (3 * (1 - 2 * mat["nu"]))
            else:
                print(
                    f"  [INFO] '{name}' has no elasticity parameters — skipping mechanical properties."
                )

            if "k" in mat:
                if isinstance(mat["k"], dict) and \
                        str(mat["k"].get("type", "")).lower() in ("neural_network", "nn"):
                    print(f"  → k defined as neural network: {mat['k'].get('weights')}")
                    from z3st.models.nn_conductivity import load_from_card
                    mat["_k_nn"] = load_from_card(mat["k"])
                elif isinstance(mat["k"], str):
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
            # so users don't hit a confusing TypeError downstream. A symbolic
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
                    pass  # "module.func" symbolic path, handled below

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

            # Material inelastic eigenstrain
            # A material card may expose an ``eigenstrain`` callable "module.func"
            # It is resolved here and consumed by MechanicalModel.eigenstrain
            if isinstance(mat.get("eigenstrain"), str):
                print(f"  → eigenstrain defined as callable: {mat['eigenstrain']}")
                mat["_eigenstrain_func"] = self.resolve_function(mat["eigenstrain"])

            # Radial power form factor f(r, bu), source-bus analogue of the eigenstrain callable
            # Resolved here and consumed by set_power
            if isinstance(mat.get("radial_profile"), str):
                print(f"  → radial_profile defined as callable: {mat['radial_profile']}")
                mat["_radial_profile_func"] = self.resolve_function(mat["radial_profile"])

            # Axial power form factor f(z)
            if isinstance(mat.get("axial_profile"), str):
                print(f"  → axial_profile defined as callable: {mat['axial_profile']}")
                mat["_axial_profile_func"] = self.resolve_function(mat["axial_profile"])

            if mat.get("cracking") is not None:
                if str(mat["cracking"]).lower() != "isotropic":
                    raise ValueError(
                        f"Material '{name}': cracking model '{mat['cracking']}' unknown "
                        f"(only 'isotropic' is implemented)."
                    )
                for key in ("cracking_lhr0", "cracking_n0", "cracking_n_inf", "cracking_tau"):
                    if key in mat:
                        mat[key] = float(mat[key])
                if "lmbda" in mat:
                    mat["lmbda"] = dolfinx.fem.Constant(
                        self.mesh, dolfinx.default_scalar_type(mat["lmbda"]))
                    mat["G"] = dolfinx.fem.Constant(
                        self.mesh, dolfinx.default_scalar_type(mat["G"]))
                print(f"  → cracking: Isotropic softening "
                      f"(LHR0 = {float(mat.get('cracking_lhr0', 5.0e3))/1e3:.1f} kW/m, "
                      f"n_inf = {float(mat.get('cracking_n_inf', 12.0)):.0f})")

            if mat.get("creep") is not None:
                if str(mat["creep"]).lower() != "norton":
                    raise ValueError(
                        f"Material '{name}': creep model '{mat['creep']}' unknown "
                        f"(only 'norton' is implemented)."
                    )
                for key in ("creep_A0", "creep_n", "creep_Q"):
                    if key not in mat:
                        raise ValueError(f"Material '{name}': creep requires '{key}'.")
                    mat[key] = float(mat[key])
                # Optional irradiation creep ε̇_irr = B·φ·σ_eq: both keys or neither.
                has_B, has_phi = "creep_irr_B" in mat, "fast_flux" in mat
                if has_B != has_phi:
                    raise ValueError(
                        f"Material '{name}': irradiation creep requires BOTH "
                        f"'creep_irr_B' and 'fast_flux' (got only one)."
                    )
                if has_B:
                    mat["creep_irr_B"] = float(mat["creep_irr_B"])
                    mat["fast_flux"] = float(mat["fast_flux"])
                    print(f"  → irradiation creep: B = {mat['creep_irr_B']:.3e} Pa^-1/(n/m^2), "
                          f"phi = {mat['fast_flux']:.3e} n/(m^2.s)")
                if constitutive_mode != "lame":
                    raise ValueError(
                        f"Material '{name}': creep is only supported with the "
                        f"'lame' constitutive route (got '{constitutive_mode}')."
                    )
                if self.on.get("damage", False) or self.on.get("plasticity", False):
                    raise ValueError(
                        "Creep cannot yet be combined with damage or plasticity "
                        "in the same run."
                    )
                print(f"  → creep: Norton, A0 = {mat['creep_A0']:.3e} Pa^-n/s, "
                      f"n = {mat['creep_n']:.2f}, Q = {mat['creep_Q']:.3e} J/mol")

            mat["__label__"] = name
            self.materials[name] = mat
            # Full per-key material dump only under --debug (CODE-P2-4)
            if "--debug" in sys.argv:
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
        self.burnup = None
        self.gas_swelling = None      # SCIANTIX total gaseous swelling ΔV/V (eigenstrain bus)
        self.fg_fields = None         # dict of per-dof FG concentration Functions (at/m^3, output)
        self.sciantix_field = None    # the per-dof SciantixField driver
        self._sciantix_dofs = None    # V_t dof indices the field covers
        self.T = None
        self.u = None
        self.D = None
        self.c = None
        self.porosity = None

        # Temperature
        if self.on.get("thermal", False):
            self.q_third = dolfinx.fem.Function(self.V_t, name="q_third")
            self.q_third.x.array[:] = 0.0
            self.set_power()

            # Burnup state field (MWd/kgU). A fissile material accumulates its own
            # local burnup from the deposited fission power (see update_state).
            # Created once, here, and *never* zeroed — it is history, not a per-step
            # source. Absent when no material is fissile, so non-fuel runs pay
            # nothing.
            if any(m.get("fissile", False) for m in self.materials.values()):
                self.burnup = dolfinx.fem.Function(self.V_t, name="Burnup")
                self.burnup.x.array[:] = 0.0
                print("Initialized burnup field (fissile material present).")

                # SCIANTIX gaseous-swelling field (fission-gas coupling, opt-in via
                # models.fission_gas; default off).
                if self.on.get("fission_gas", False):
                    self._init_sciantix_field()

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

        # Porosity variables
        if self.on.get("porosity", False):
            print("\nInitializing the porosity field...")
            self.porosity = dolfinx.fem.Function(self.V_p, name="Porosity")
            self.porosity.x.array[:] = 0.0
            self.porosity_n = dolfinx.fem.Function(self.V_p, name="Porosity_old")
            self.porosity_n.x.array[:] = 0.0

            print("\nSetting porosity initial conditions...")
            self.set_porosity_initial_conditions()

        # Material properties
        for name, mat in self.materials.items():
            if "_k_func" in mat and self.T:
                k_func = mat["_k_func"]
                mat["k"] = k_func(self.T)
                print("\nk expression for", name, "→", mat["k"])

            # Neural-network conductivity
            if "_k_nn" in mat and self.T is not None:
                k_fn = dolfinx.fem.Function(self.V_t, name=f"k_nn_{name}")
                k_fn.x.array[:] = mat["_k_nn"](self.T.x.array)
                k_fn.x.scatter_forward()
                mat["k"] = k_fn
                print(f"\nk neural-network field for {name} → Function on V_t "
                      f"(min={k_fn.x.array.min():.3f}, max={k_fn.x.array.max():.3f} W/m/K)")

            # Porosity-dependent conductivity
            if mat.get("thermal_conductivity_model") == "kato_porosity" and self.T is not None:
                k_fn = dolfinx.fem.Function(self.V_t, name=f"k_porosity_{name}")
                mat["k"] = k_fn
                print(f"\nInitialized porosity-dependent thermal conductivity field for {name}")

            # Temperature-dependent elastic constants: build lmbda/G/bulk_modulus
            # as UFL expressions in the live T field, so the per-iteration T
            # propagates by reference into both the mechanical form and the
            # (pre-compiled) output-writer stress expression.
            if "_E_func" in mat or "_nu_func" in mat:
                if getattr(self, "T", None) is None:
                    raise ValueError(
                        f"Material '{name}': temperature-dependent E/nu requires an "
                        f"active thermal field (set models.thermal: true)."
                    )
                E_T = mat["_E_func"](self.T) if "_E_func" in mat else mat["E"]
                nu_T = mat["_nu_func"](self.T) if "_nu_func" in mat else mat["nu"]
                mat["lmbda"] = E_T * nu_T / ((1 + nu_T) * (1 - 2 * nu_T))
                mat["G"] = E_T / (2 * (1 + nu_T))
                mat["bulk_modulus"] = E_T / (3 * (1 - 2 * nu_T))
                print(f"\nE/nu expression for {name} → lmbda, G built as UFL(T)")

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

                if dmg_type == "AT1" and "E" in mat:
                    mat["sigma_c"] = ufl.sqrt((3 * mat["E"] * mat["Gc"]) / (8 * lc))
                    print(f"  - Material '{name}': sigma_c (AT1) evaluated from Gc expression")

        if self.on.get("porosity", False):
            self.update_porosity_dependent_properties(self.T, self.porosity)


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

                # Power form factors. A fissile material may shape its own
                # volumetric source through the callables ``radial_profile``
                # f(r, bu) (intended for a TUBRNP-style rim profile later) and/or
                # ``axial_profile`` f(z) (e.g. the chopped cosine). Each callable
                # receives the dof coordinates and the current local burnup. The
                # composite f_r·f_z is normalised once to nodal mean 1, so the
                # shaping redistributes the linear heat rate without changing its
                # integral. Default (no callables): f ≡ 1, the flat source.
                shape = np.ones(len(dofs))
                rprof = mat.get("_radial_profile_func")
                zprof = mat.get("_axial_profile_func")
                if rprof is not None or zprof is not None:
                    coords = self.V_t.tabulate_dof_coordinates()[dofs]
                    bu_vals = (
                        self.burnup.x.array[dofs]
                        if self.burnup is not None
                        else np.zeros(len(dofs))
                    )
                    if rprof is not None:
                        shape = shape * np.asarray(rprof(coords, bu_vals, mat, model=self), dtype=float)
                    if zprof is not None:
                        shape = shape * np.asarray(zprof(coords, bu_vals, mat, model=self), dtype=float)
                    # Nodal-mean normalisation (area-weighted refinement to come
                    # with the TUBRNP profile).
                    mean = shape.mean()
                    if mean > 0:
                        shape = shape / mean

                # Accumulate: if multiple sources are configured on the same
                # material (e.g. fissile + gamma_heating below), they should
                # add — not overwrite.
                self.q_third.x.array[dofs] += q_val * shape
                print(f"  q_third += {q_val:.3e} W/m³ × f(r,bu)·f(z) (fissile, mean f = 1)")
                print(f"  Heat flux = {self.lhr / self.perimeter:.3e} W/m2")

                # Integrated-power diagnostic: the exact FE integral of the
                # fissile source over this material, with the regime weight
                # (2πr in axisymmetric). For a rod this should track LHR·Lz;
                # a radially peaked profile deviates slightly because the mean-1
                # normalisation is nodal, not area-weighted. The form is compiled
                # once and cached (q_third updates in place).
                if not hasattr(self, "_power_forms"):
                    self._power_forms = {}
                if name not in self._power_forms:
                    x_sc = ufl.SpatialCoordinate(self.mesh)
                    w_int = 2.0 * ufl.pi * x_sc[0] if self.regime == "axisymmetric" else 1.0
                    dx_mat = ufl.Measure(
                        "dx", domain=self.mesh,
                        subdomain_data=self.cell_tags,
                        subdomain_id=self.label_map[name],
                    )
                    self._power_forms[name] = dolfinx.fem.form(w_int * self.q_third * dx_mat)
                P_int = dolfinx.fem.assemble_scalar(self._power_forms[name])
                P_int = self.mesh.comm.allreduce(P_int, op=MPI.SUM)
                unit = {"axisymmetric": "W", "3d": "W", "2d": "W/m",
                        "plane_stress": "W/m", "1d": "W/m²"}.get(self.regime, "W")
                print(f"  [INFO] Integrated fissile power in {name}: {P_int:.6e} {unit}")

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
                # Per-material reference surface for the cylindrical/spherical
                # decay correlation. Defaults to the geometry inner_radius
                # (backward-compatible). A material that sits inboard of the
                # geometry reference (e.g. a thermal shield in front of the
                # vessel) sets `gamma_inner_radius` so its K_0 profile is
                # normalised at its own inner surface rather than the vessel's.
                gamma_Ri = float(mat.get("gamma_inner_radius", self.inner_radius))

                def f(x, q_third_0=q_third_0, mu=mu, gamma_Ri=gamma_Ri):
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

                        return q_third_0 * sp.k0(mu * radius) / sp.k0(mu * gamma_Ri)
                    elif self.geometry_type == "sphere":
                        r = np.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)
                        return (
                            q_third_0
                            * (gamma_Ri / r)
                            * np.exp(-mu * (r - gamma_Ri))
                        )

                f_func = dolfinx.fem.Function(self.V_t)
                f_func.interpolate(f)
                # Accumulate (see CODE-P1-10): allows fissile + gamma_heating
                # on the same material to combine rather than overwrite.
                self.q_third.x.array[dofs] += f_func.x.array[dofs]

        self.q_third.x.scatter_forward()

    def update_state(self, dt):
        """Advance each material's own history over a step of ``dt`` seconds.
        Called once per step, *after* the solve.

        Burnup: a fissile material accumulates its local burnup from the deposited
        fission power. The volumetric source ``q_third`` [W/m³] is the energy
        deposited per unit fuel volume per second; dividing by the heavy-metal mass
        density ``ρ·HM_frac`` gives the specific power per unit heavy metal [W/kgU],
        whose time integral is the burnup. The unit conversion W·s/kgU → MWd/kgU
        divides by 86400 s/day × 1e6 W/MW = 8.64e10::

            Δbu = q_third · dt / (ρ · HM_frac · 8.64e10)   [MWd/kgU]

        No feedback is applied here — burnup is recorded. Downstream behaviours
        that consume it (fuel-k(bu), swelling(bu,T), FGR) read this field, so a
        fissile case without such behaviour is unaffected in its solve.
        """
        if self.burnup is None or dt <= 0.0:
            return

        SECONDS_PER_MWD = 8.64e10  # 86400 s/day × 1e6 W/MW

        # Capture burnup at step-start so SCIANTIX gets the (old, new) pair
        bu_old = (self.burnup.x.array[self._sciantix_dofs].copy()
                  if self.sciantix_field is not None else None)

        for name, mat in self.materials.items():
            if not mat.get("fissile", False):
                continue
            rho = mat.get("rho")
            if rho is None:
                print(f"  [update_state] '{name}' is fissile but has no 'rho'; "
                      f"skipping burnup accumulation.")
                continue
            hm = float(mat.get("heavy_metal_fraction", 0.8815))
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_t)
            q = self.q_third.x.array[dofs]
            self.burnup.x.array[dofs] += q * dt / (float(rho) * hm * SECONDS_PER_MWD)

        self.burnup.x.scatter_forward()
        print(f"[update_state] burnup max = {self.burnup.x.array.max():.4e} MWd/kgU")

        # --. SCIANTIX gaseous swelling (opt-in; rides the eigenstrain bus) --..
        if self.sciantix_field is not None:
            self._update_sciantix(dt, bu_old)

    def _init_sciantix_field(self):
        """Build the per-dof SCIANTIX field over the fissile region (opt-in coupling).

        One SCIANTIX integration point per ``V_t`` dof of every fissile material,
        seeded from the case's ``input_initial_conditions.txt``. Requires
        ``libsciantix.so`` (``config.sciantix_lib`` or ``$SCIANTIX_LIB``) plus
        ``input_settings.txt`` and ``input_initial_conditions.txt`` in the run
        directory — the same files a SCIANTIX standalone run needs.
        """
        from z3st.coupling.sciantix.sciantix_binding import SciantixField

        dof_lists = [
            self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_t)
            for name, mat in self.materials.items() if mat.get("fissile", False)
        ]
        self._sciantix_dofs = (
            np.unique(np.concatenate(dof_lists)) if dof_lists else np.array([], dtype=np.int64)
        )
        self.gas_swelling = dolfinx.fem.Function(self.V_t, name="Gaseous swelling")
        self.gas_swelling.x.array[:] = 0.0

        # Per-dof total fission-gas (Xe + Kr) concentration fields [at/m^3], for the
        # Paraview output and the FG plots — one Function per SCIANTIX gas state.
        # Names mirror SciantixField.gas_concentrations() keys.
        self.fg_fields = {
            "produced":       dolfinx.fem.Function(self.V_t, name="FG produced (at_m3)"),
            "in_grain":       dolfinx.fem.Function(self.V_t, name="FG in grain (at_m3)"),
            "grain_boundary": dolfinx.fem.Function(self.V_t, name="FG at grain boundary (at_m3)"),
            "released":       dolfinx.fem.Function(self.V_t, name="FG released (at_m3)"),
        }
        for fn in self.fg_fields.values():
            fn.x.array[:] = 0.0

        n = int(self._sciantix_dofs.size)
        if n == 0:
            print("[sciantix] no fissile dofs; gaseous-swelling field stays zero.")
            return
        print(f"[sciantix] building SCIANTIX field over {n} fissile dofs ...")
        self.sciantix_field = SciantixField(
            n, libpath=self.sciantix_lib, ic_path=self.sciantix_ic
        )
        print("[sciantix] field initialised (library loaded, initial conditions seeded).")

    def _update_sciantix(self, dt, bu_old):
        """Advance every SCIANTIX point by ``dt`` and refresh ``gas_swelling``.

        Local conditions per dof: temperature from ``self.T``; volumetric fission
        rate from the deposited power, ``fission_rate = q_third / E_fission``
        [fiss/m³ s] with ``E_fission`` from ``config.sciantix_energy_per_fission``;
        and the host burnup pair (``bu_old`` captured at step-start, ``bu_new`` the
        just-updated ``self.burnup``) transferred to SCIANTIX (it does not compute
        burnup in a coupling build — Z3ST's RADAR model owns it).
        """
        dofs = self._sciantix_dofs
        T = self.T.x.array[dofs]
        fission_rate = self.q_third.x.array[dofs] / self.sciantix_energy_per_fission
        bu_new = self.burnup.x.array[dofs]
        gs = self.sciantix_field.step(dt, T, fission_rate,
                                      burnup_old=bu_old, burnup_new=bu_new)
        self.gas_swelling.x.array[dofs] = gs
        self.gas_swelling.x.scatter_forward()
        print(f"[update_state] gaseous swelling max = {self.gas_swelling.x.array.max():.4e} (ΔV/V)")

        # Fission-gas (Xe + Kr) concentration fields for output / FG plots.
        conc = self.sciantix_field.gas_concentrations()
        for key, fn in self.fg_fields.items():
            fn.x.array[dofs] = conc[key]
            fn.x.scatter_forward()
        prod = self.fg_fields["produced"].x.array[dofs]
        rel = self.fg_fields["released"].x.array[dofs]
        fgr = float(rel.sum() / prod.sum()) if prod.sum() > 0 else 0.0
        print(f"[update_state] fission gas: produced max = {prod.max():.4e} at/m³, "
              f"fuel-avg FGR = {fgr:.4f}")

    _SNAPSHOT_FIELDS = (
        "T", "u", "D", "H", "burnup", "gas_swelling", "c", "c_n",
        "p", "ep", "p_n", "ep_n",
    )  # dolfinx Functions
    _SNAPSHOT_DICTS = ("eps_cr", "_dgamma0")
    _SNAPSHOT_MATERIAL_KEYS = ("E", "nu", "bulk_modulus", "_lhr_max")
    _SNAPSHOT_MATERIAL_CONSTANTS = ("lmbda", "G")

    def snapshot_state(self):
        """Deep-copy every step-level state field so a time step can be rolled
        back and retried at a smaller dt (adaptive time-stepping, piece B).

        The roster of what is captured lives in the class constants
        ``_SNAPSHOT_FIELDS`` / ``_SNAPSHOT_DICTS`` / ``_SNAPSHOT_MATERIAL_*`` —
        extend those when adding new persistent state. Captured (only those
        present, per active physics):

        - primary fields T, u, D and the crack-driving history H;
        - the burnup accumulator and the cluster pair c / c_n;
        - the plasticity history p, ep, p_n, ep_n;
        - the per-material creep dicts eps_cr and _dgamma0 (predictor,
          updated every iteration, so polluted even by a failed attempt);
        - per-material cracking scalars (``_lhr_max`` is a running max that
          does NOT decrease, so without rollback a bisected sub-step inherits
          the failed attempt's stiffness degradation) and the live lmbda/G
          Constants.

        NOT captured: iteration scratch (``_aitken_R_prev``, ``_aitken_p_R_prev``,
        ``_h_gap_prev``) — solve_staggered resets those at entry.

        IMPORTANT: take the snapshot at the very start of a (sub)step, BEFORE
        parameters()/set_power()/update_state() run, so ``_lhr_max`` and burnup
        are captured at their last-converged values.
        """
        snap = {"fields": {}, "dicts": {}, "materials": {}}
        for name in self._SNAPSHOT_FIELDS:
            fn = getattr(self, name, None)
            if isinstance(fn, dolfinx.fem.Function):
                snap["fields"][name] = fn.x.array.copy()
        for name in self._SNAPSHOT_DICTS:
            d = getattr(self, name, None)
            if isinstance(d, dict):
                snap["dicts"][name] = {
                    k: v.x.array.copy()
                    for k, v in d.items()
                    if isinstance(v, dolfinx.fem.Function)
                }
        for mname, mat in getattr(self, "materials", {}).items():
            msnap = {k: mat[k] for k in self._SNAPSHOT_MATERIAL_KEYS if k in mat}
            for k in self._SNAPSHOT_MATERIAL_CONSTANTS:
                v = mat.get(k)
                if hasattr(v, "value"):  # live dolfinx Constant
                    msnap[k + ".value"] = float(v.value)
            if msnap:
                snap["materials"][mname] = msnap
        # SCIANTIX per-point C-state (variables/diffusion_modes/history) — the
        # gas_swelling dolfinx field alone is not enough to roll back a step.
        if self.sciantix_field is not None:
            snap["sciantix"] = self.sciantix_field.snapshot()
        return snap

    def restore_state(self, snap):
        """Inverse of :meth:`snapshot_state`: write every captured field, dict
        and material scalar back in place, undoing a failed or oversized step so
        it can be retried at a smaller dt. scatter_forward keeps ghost dofs
        consistent after the in-place array overwrite."""
        for name, arr in snap.get("fields", {}).items():
            fn = getattr(self, name)
            fn.x.array[:] = arr
            fn.x.scatter_forward()
        for name, saved in snap.get("dicts", {}).items():
            d = getattr(self, name)
            for k, arr in saved.items():
                d[k].x.array[:] = arr
                d[k].x.scatter_forward()
        for mname, msnap in snap.get("materials", {}).items():
            mat = self.materials[mname]
            for key, val in msnap.items():
                if key.endswith(".value"):
                    target = mat.get(key[: -len(".value")])
                    if hasattr(target, "value"):
                        target.value = val
                else:
                    mat[key] = val
        if "sciantix" in snap and self.sciantix_field is not None:
            self.sciantix_field.restore(snap["sciantix"])

    def solve(self, max_iters=100, dt=0.0):
        print("\n")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print(f"--. spine - solve --..")
        print("--.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---")
        print("\n")

        print(f"Current step = {self.current_step} | dt = {dt:.2e} s")
        print(f"Coupling = {self.coupling}")

        if self.coupling == "staggered":
            # Return the convergence result so the time loop can react to a
            # stalled step: True on convergence, False if it exhausts max_iter.
            return self.solve_staggered(
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
            
            if self.on.get("mechanical", False) and self.creep_active(mat):
                T_field = self.T if self.on.get("thermal", False) else None
                sigma_cr, eps_el_cr = self.creep_output_stress(self.u, mat, T_field)
                self.stress[name] = sigma_cr
                self.energy_density[name] = 0.5 * ufl.inner(sigma_cr, eps_el_cr)

