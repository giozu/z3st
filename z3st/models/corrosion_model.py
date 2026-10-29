# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.4.1 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import math
import os
import sys

import dolfinx
import dolfinx.fem.petsc
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from scipy.special import erf

R_GAS = 8.314462618  # J/mol/K
KB_EV = 8.617333e-5  # eV/K
FARADAY = 96485.33  # C/mol
M_CR = 51.996e-3  # kg/mol
M_NI = 58.693e-3  # kg/mol


# --.. ..- .-.. .-.. --- material functions --.. ..- .-.. .-.. ---
def arrhenius(D0, Q, T, exp=math.exp):
    """D0 exp(-Q / (R T)). ``exp=ufl.exp`` gives a UFL expression."""
    return D0 * exp(-Q / (R_GAS * T))


def hart_diffusivity(card, T, exp=math.exp):
    """Effective Cr diffusivity (m2/s), Hart approximation (Bhave et al. 2023, Eq. 46):
    D_b (1 - d_gb/G) + D_gb d_gb/G. Bulk only when the card has no GB data."""
    D_b = arrhenius(card["D0_bulk"], card["Q_bulk"], T, exp)
    if "gb_width" not in card or "grain_size" not in card:
        return D_b
    f_gb = card["gb_width"] / card["grain_size"]
    D_gb = arrhenius(card["D0_gb"], card["Q_gb"], T, exp)
    return D_b * (1.0 - f_gb) + D_gb * f_gb


def check_diffusivity(D, label):
    """Warn when D (m2/s) is outside the range expected for solid-state diffusion,
    which usually means the card was left in um2/s."""
    if not 1e-24 <= D <= 1e-10:
        print(f"  [WARNING] Cr diffusivity for '{label}' is {D:.3e} m2/s, outside "
              f"[1e-24, 1e-10] m2/s. Check that D0 is in m2/s.")


def mass_fraction_cr(x, M_Cr=M_CR, M_Ni=M_NI):
    """Cr mass fraction of a binary Ni-Cr alloy from the Cr atomic fraction."""
    return x * M_Cr / (x * M_Cr + (1.0 - x) * M_Ni)


def atomic_fraction_cr(w, M_Cr=M_CR, M_Ni=M_NI):
    """Cr atomic fraction of a binary Ni-Cr alloy from the Cr mass fraction."""
    return (w / M_Cr) / (w / M_Cr + (1.0 - w) / M_Ni)


# --.. ..- .-.. .-.. --- salt equilibrium --.. ..- .-.. .-.. ---
def _gibbs_fcc_ni(T):
    """G of Ni on FCC (eV/atom), Bhave et al. 2023, Eq. A.2."""
    return (-5179.159 + 117.854 * T - 22.096 * T * math.log(T) - 4.8407e-3 * T**2) / FARADAY


def _gibbs_cr_in_ni(T):
    """G of Cr on the FCC Ni lattice (eV/atom), Eqs. A.3 and A.8."""
    G_fcc_cr = (-1572.94 + 157.643 * T - 26.908 * T * math.log(T) + 1.89435e-3 * T**2
                - 1.47721e-6 * T**3 + 139250.0 / T)
    H_xs, S_xs = -1.5644e4, -1.5601
    return (G_fcc_cr + H_xs - S_xs * T) / FARADAY


def _gibbs_bcc_cr(T):
    """G of Cr on BCC (eV/atom), Eq. 32."""
    return (-1.53e-11 * T**3 + 1.96e-8 * T**2 - 2.78e-4 * T * math.log(T)
            + 1.63e-3 * T - 9.18e-2 + 1.44 / T)


def redox_surface_fraction(T, E_F, E0_F, E0_Cr=-0.39, E0_Ni=0.473,
                           x_Cr_salt=2.5e-6, x_Ni_salt=3.69e-11):
    """Cr and Ni atomic fractions of the alloy surface in equilibrium with the salt
    (Bhave et al. 2023, Eqs. 30-31, B.2-B.4, B.8-B.9).

    The salt redox state enters as the fluoride-potential shift E_F - E0_F, applied
    to both Cr2+ and Ni2+ as Eq. 30 states it. With x_Ni_salt fixed and A_Ni >> 1
    (Bhave's conditions, E_F above about 3 V) the shift cancels and the result
    depends only on x_Cr_salt / x_Ni_salt; at lower E_F it does not."""
    kT = KB_EV * T
    shift = E_F - E0_F
    G0_v = 1.56 - 3.3 * kT
    G0_cr2 = _gibbs_bcc_cr(T) + 2.0 * (E0_Cr + shift)
    G0_ni2 = _gibbs_fcc_ni(T) + 2.0 * (E0_Ni + shift)
    mu_cr = G0_cr2 + kT + kT * math.log(x_Cr_salt)
    mu_ni = G0_ni2 + kT + kT * math.log(x_Ni_salt)
    A_cr = math.exp((mu_cr - (_gibbs_cr_in_ni(T) - G0_v)) / kT)
    A_ni = math.exp((mu_ni - (_gibbs_fcc_ni(T) - G0_v)) / kT)
    x_cr = A_cr / (1.0 + A_ni + A_cr)
    x_ni = A_ni / (1.0 + A_ni + A_cr)
    if not math.isfinite(x_cr) or x_cr > 0.5:
        raise ValueError(
            f"redox surface Cr fraction is {x_cr:.3e}: the surface is Cr-saturated. "
            f"Check the salt inputs: E_F = {E_F} V, E0_F = {E0_F} V, "
            f"x_Cr_salt = {x_Cr_salt}, x_Ni_salt = {x_Ni_salt}.")
    return x_cr, x_ni


# --.. ..- .-.. .-.. --- analytic references --.. ..- .-.. .-.. ---
def erfc_profile(d, t, D, x0, xs):
    """Cr atomic fraction at depth d below a surface held at xs, semi-infinite alloy."""
    return xs + (x0 - xs) * erf(np.asarray(d) / (2.0 * np.sqrt(D * t)))


def profile_mass_loss(t, D, x0, xs, rho, n=20000):
    """Cr mass loss per area (kg/m2) of the erfc profile, integrated in mass fraction."""
    L = 10.0 * 2.0 * np.sqrt(D * t)
    h = L / n
    d = (np.arange(n) + 0.5) * h
    return rho * np.sum(mass_fraction_cr(x0) - mass_fraction_cr(erfc_profile(d, t, D, x0, xs))) * h


def eq50_mass_loss(w0, ws, rho, D, t):
    """Cr mass loss per area (kg/m2), Bhave et al. 2023 Eq. 50."""
    return 2.0 * (w0 - ws) * rho * math.sqrt(D * t / math.pi)


# --.. ..- .-.. .-.. --- model --.. ..- .-.. .-.. ---
class CorrosionModel:
    """
    Selective Cr depletion of an alloy in contact with a molten salt (reduced model).

    Solves dx/dt = div(D(T) grad x) for the Cr atomic fraction x in the alloy, with
    x = x_surf on the salt-facing boundary and zero flux elsewhere. D(T) is the Hart
    effective diffusivity of the material card. x_surf is given (surface.type: fixed)
    or computed from equilibrium with the salt (surface.type: redox). The salt is not
    meshed. In 1D this reproduces the electrochemical phase-field model of Bhave et
    al. (2023) to within 5 % (their Sec. 4).
    """

    def __init__(self):
        self.corr_cfg = self.input_file.get("corrosion", {})
        self.corr_T = self.corr_cfg.get("T")
        if self.corr_T is not None:
            self.corr_T = float(self.corr_T)
        elif not self.on.get("thermal", False):
            print("  [ERROR] corrosion.T is required when the thermal model is off.")
            sys.exit(1)
        self.corr_threshold = float(self.corr_cfg.get("depletion_threshold", 0.843))
        self.corr_axis = int(self.corr_cfg.get("depth_axis", 0))
        self.corr_solver = self.corr_cfg.get("linear_solver", "direct_mumps")

        surf = self.corr_cfg.get("surface", {})
        surf_type = surf.get("type", "fixed")
        if surf_type == "fixed":
            self.x_surf = float(surf["x_surf"])
        elif surf_type == "redox":
            if self.corr_T is None:
                print("  [ERROR] corrosion.T is required for surface.type: redox.")
                sys.exit(1)
            optional = {k: float(surf[k]) for k in ("E0_Cr", "E0_Ni", "x_Cr_salt", "x_Ni_salt") if k in surf}
            self.x_surf, _ = redox_surface_fraction(
                self.corr_T, float(surf["E_F"]), float(surf["E0_F"]), **optional)
        else:
            print(f"  [ERROR] Unknown corrosion surface type '{surf_type}' (fixed | redox).")
            sys.exit(1)
        if not 0.0 <= self.x_surf <= 1.0:
            print(f"  [ERROR] Corrosion surface Cr fraction {self.x_surf} is outside [0, 1].")
            sys.exit(1)

        print("[CorrosionModel] Initialized")
        print(f"  [CorrosionModel] Surface condition: {surf_type}, x_surf = {self.x_surf:.4e}")

    def set_corrosion_initial_conditions(self):
        """x = x0 of each material's corrosion card. Every material needs a card."""
        for name, mat in self.materials.items():
            card = mat.get("corrosion")
            if card is None:
                print(f"  [ERROR] Material '{name}' has no 'corrosion' card.")
                sys.exit(1)
            # YAML reads 2.852e5 (no exponent sign) as a string
            mat["_corr"] = {k: float(v) for k, v in card.items()}
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_x)
            self.x_cr.x.array[dofs] = mat["_corr"]["x0"]
            self.x_cr_n.x.array[dofs] = mat["_corr"]["x0"]
            if self.corr_T is not None:
                D = hart_diffusivity(mat["_corr"], self.corr_T)
                check_diffusivity(D, name)
                print(f"  [Corrosion] '{name}': x0 = {mat['_corr']['x0']}, D({self.corr_T:.1f} K) = {D:.4e} m2/s")
            else:
                print(f"  [Corrosion] '{name}': x0 = {mat['_corr']['x0']}, D checked at the first solved temperature")
        # With T from the thermal model the initial field is T_ref, not the service
        # temperature, so the D range check waits for the first solved T.
        self._corr_D_checked = self.corr_T is not None
        self.x_cr.x.scatter_forward()
        self.x_cr_n.x.scatter_forward()
        self._x0_dofs = self.x_cr.x.array.copy()

    def set_corrosion_boundary_conditions(self, V_x):
        """Dirichlet x = x_surf on the salt-facing regions of boundary_conditions['corrosion']."""
        self.dirichlet_corrosion = []
        self._corr_surface_ids = []
        surface_dofs = []
        self._x_surf_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.x_surf))

        for label, bc_list in self.boundary_conditions.get("corrosion", {}).items():
            for bc_info in bc_list:
                region_name = bc_info.get("region")
                region_id = self.label_map.get(region_name)
                if region_id is None:
                    print(f"  [ERROR] Region '{region_name}' not found in label_map for corrosion BC.")
                    sys.exit(1)
                if bc_info.get("type") != "Dirichlet":
                    print(f"  [ERROR] Corrosion BC on '{region_name}' must be Dirichlet (salt-facing surface).")
                    sys.exit(1)
                facets = self.facet_tags.find(region_id)
                if self.mesh.comm.allreduce(len(facets), op=MPI.SUM) == 0:
                    print(f"  [ERROR] Corrosion BC region '{region_name}' has no boundary facets "
                          f"(is it a volume label?).")
                    sys.exit(1)
                dofs = dolfinx.fem.locate_dofs_topological(V_x, self.fdim, facets)
                self.dirichlet_corrosion.append(dolfinx.fem.dirichletbc(self._x_surf_const, dofs, V_x))
                surface_dofs.append(dofs)
                self._corr_surface_ids.append(region_id)
                print(f"  [INFO] Corrosion surface on '{label}' at region '{region_name}': x = {self.x_surf:.4e}")

        if not self._corr_surface_ids:
            print("  [ERROR] No salt-facing surface: boundary_conditions has no 'corrosion' entry.")
            sys.exit(1)
        self._corr_surface_dofs = np.unique(np.concatenate(surface_dofs))
        for name, mat in self.materials.items():
            if self.x_surf > mat["_corr"]["x0"]:
                print(f"  [WARNING] x_surf = {self.x_surf:.4e} exceeds x0 = {mat['_corr']['x0']} "
                      f"of '{name}': the salt enriches the alloy in Cr.")

    # --.. ..- .-.. .-.. --- staggered step --.. ..- .-.. .-.. ---
    def _corrosion_step(self, x_new, x_old, dt, T_current):
        """
        Backward-Euler step of the Cr diffusion equation.

        ``x_old`` is the t^n level (``self.x_cr_n``) and is not modified: the step is
        re-solved from it on every staggered iteration, as in ``_cluster_step``.
        One-way coupled: it reads T, never feeds back, and never blocks convergence.
        """
        if dt <= 0.0:
            x_new.x.array[:] = x_old.x.array
            return

        x, v = ufl.TrialFunction(self.V_x), ufl.TestFunction(self.V_x)
        dt_c = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(dt))
        if self.on.get("thermal", False) and T_current is not None:
            T = T_current
            if not self._corr_D_checked:
                T_max = self.mesh.comm.allreduce(float(np.max(T_current.x.array)), op=MPI.MAX)
                for name, mat in self.materials.items():
                    check_diffusivity(hart_diffusivity(mat["_corr"], T_max), name)
                self._corr_D_checked = True
        else:
            T = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.corr_T))

        a, L = 0, 0
        for name, mat in self.materials.items():
            D = hart_diffusivity(mat["_corr"], T, exp=ufl.exp)
            dxm = self.dx_tags[self.label_map[name]]
            a += (x / dt_c) * v * dxm + D * ufl.inner(ufl.grad(x), ufl.grad(v)) * dxm
            L += (x_old / dt_c) * v * dxm

        problem = dolfinx.fem.petsc.LinearProblem(
            a, L, u=x_new, bcs=self.dirichlet_corrosion,
            petsc_options=self.get_solver_options("corrosion", self.corr_solver),
            petsc_options_prefix="corrosion_",
        )
        problem.solve()
        x_new.x.scatter_forward()

    # --.. ..- .-.. .-.. --- diagnostics --.. ..- .-.. .-.. ---
    def _corrosion_temperature(self):
        if self.on.get("thermal", False) and self.T is not None:
            return self.T
        return dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(self.corr_T))

    def corrosion_diagnostics(self):
        """Cr mass loss per salt-facing area (kg/m2), depletion depth (m) and mean
        outward Cr flux through the surface (atomic fraction x m/s)."""
        comm = self.mesh.comm
        T = self._corrosion_temperature()

        def integrate(expr):
            return comm.allreduce(dolfinx.fem.assemble_scalar(dolfinx.fem.form(expr)), op=MPI.SUM)

        ds_s = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.facet_tags,
                           subdomain_id=tuple(self._corr_surface_ids))
        n = ufl.FacetNormal(self.mesh)
        area = integrate(dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(1.0)) * ds_s)

        loss, flux = 0.0, 0.0
        for name, mat in self.materials.items():
            card = mat["_corr"]
            tag = self.label_map[name]
            dx_m = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_tags, subdomain_id=tag)
            loss += card["rho"] * integrate((mass_fraction_cr(card["x0"]) - mass_fraction_cr(self.x_cr)) * dx_m)
            # Surface facets of this material only: the indicator is 1 on its cells
            ind = dolfinx.fem.Function(self.Q)
            cells = self.cell_tags.find(tag)
            ind.x.array[dolfinx.fem.locate_dofs_topological(self.Q, self.mesh.topology.dim, cells)] = 1.0
            ind.x.scatter_forward()
            D = hart_diffusivity(card, T, exp=ufl.exp)
            flux += integrate(-ind * D * ufl.dot(ufl.grad(self.x_cr), n) * ds_s)

        # Depth along depth_axis, from the mean coordinate of the surface dofs
        coords = self.V_x.tabulate_dof_coordinates()[:, self.corr_axis]
        n_owned = self.V_x.dofmap.index_map.size_local
        surf = self._corr_surface_dofs[self._corr_surface_dofs < n_owned]
        X_surf = comm.allreduce(coords[surf].sum(), op=MPI.SUM) / comm.allreduce(surf.size, op=MPI.SUM)
        x = self.x_cr.x.array[:n_owned]
        depleted = x < self.corr_threshold * self._x0_dofs[:n_owned]
        local_depth = float(np.max(np.abs(coords[:n_owned][depleted] - X_surf))) if depleted.any() else 0.0
        depth = comm.allreduce(local_depth, op=MPI.MAX)

        return {"mass_loss": loss / area, "depth": depth, "flux": flux / area}

    def write_corrosion_history(self, t, step, output_dir="output"):
        """Append t, mass loss (kg/m2 and mg/cm2), depth (m) and flux to output/corrosion.tsv."""
        d = self.corrosion_diagnostics()
        if self.mesh.comm.rank != 0:
            return
        path = os.path.join(output_dir, "corrosion.tsv")
        mode = "w" if step == 0 else "a"
        with open(path, mode) as f:
            if step == 0:
                f.write("t_s\tmass_loss_kg_m2\tmass_loss_mg_cm2\tdepth_m\tflux\n")
            f.write(f"{t:.6e}\t{d['mass_loss']:.6e}\t{d['mass_loss'] * 100.0:.6e}\t"
                    f"{d['depth']:.6e}\t{d['flux']:.6e}\n")
        print(f"  [Corrosion] t = {t:.3e} s: mass loss = {d['mass_loss'] * 100.0:.4e} mg/cm2, "
              f"depth = {d['depth'] * 1e6:.3f} um")
