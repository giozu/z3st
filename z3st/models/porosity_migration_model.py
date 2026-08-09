# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

import dolfinx
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from z3st.core.solver import aitken_omega


class PorosityMigrationModel:
    """
    Porosity migration model representing thermal-gradient-driven pore transport.

    Two discretisations of the same transport equation. The CG path (default) uses
    SU/SUPG stabilisation and reproduces Barani et al. (2022). The DG path uses an
    upwind flux with an optional SIPG block, SSP-RK3 in time, a Kuzmin vertex
    limiter for the discrete maximum principle, and a conservative saturation cap;
    it is exercised by cases/verification/fuel/porosity_migration_dg.

    The DG derivation -- weak form, facet-flux consistency, CFL estimate,
    boundedness proof and the mass-conservation statement -- lives in a separate
    write-up that is not part of this repository.
    """

    def __init__(self):
        print("[PorosityMigrationModel] initializer")

        self.porosity_cfg = self.input_file.get("porosity", {})

        # Sens (1972) pore-velocity correlation, J. Nucl. Mater. 43, 293, Eqs. 20-24,
        # as used by Novascone et al. (2018) Eq. 4 and Barani et al. (2022) Eq. 9
        # (the benchmark reproduced here), in SI units:
        #   v = v0 (c1 + c2 T + c3 T^2 + c4 T^3) T^-2.5 exp(-Hs/RT) grad T   [m/s]
        # - (c1..c4) are the UO2 molecular-volume polynomial from the Christensen
        #   density data, Sens Eq. 24: rho_T = 10.97/(c1 + c2 T + c3 T^2 + c4 T^3).
        # - Hs and the vapour-pressure pre-exponential ln(p0) = 33.7085 [dyn/cm^2]
        #   are Sens Table 1 (Ohse data, 1500-2000 K range): the Arrhenius term
        #   exp(-Hs/RT) is the UO2 vapour pressure p = p0 exp(-Hs/RT), Sens Eq. 18.
        # - v0 is the single lumped prefactor (Sens Eq. 21 geometric/kinetic factor
        #   4 p0 Omega_T DeltaH ... / (3 pi sigma^2 P N ...), which Novascone calls
        #   "the product of many other constants"). v0 = 1.303427e8 reproduces
        #   Sens' own Table 2 disc-pore velocity (4.2e-11 m/s at T=2000 K,
        #   gradT=2500 K/cm) to ~8%, and yields the published central-void radius
        #   (~0.2 r/Ro).
        self.v0 = float(self.porosity_cfg.get("v0", 1.303427e8))
        self.c1 = float(self.porosity_cfg.get("c1", 0.988))
        self.c2 = float(self.porosity_cfg.get("c2", 6.395e-6))
        self.c3 = float(self.porosity_cfg.get("c3", 3.543e-9))
        self.c4 = float(self.porosity_cfg.get("c4", 3.0e-12))
        self.H_s = float(self.porosity_cfg.get("Hs", 5.98e5))   # heat of vaporisation [J/mol] (Sens Table 1)

        print("[PorosityMigrationModel] options loaded:")
        for key, val in [
            ("v0", self.v0), ("c1", self.c1), ("c2", self.c2),
            ("c3", self.c3), ("c4", self.c4), ("Hs", self.H_s),
        ]:
            print(f"  {key:<20}: {val}")

    def set_porosity_initial_conditions(self):
        """
        Sets the initial condition for the porosity field using material cards.
        """
        print("\n[PorosityMigration] Setting initial porosity conditions...")
        self.porosity.x.array[:] = 0.0

        for name, mat in self.materials.items():
            dofs = self.mgr.locate_domain_dofs(label=self.label_map[name], V=self.V_p)
            p_init = float(mat.get("initial_porosity", 0.0))
            self.porosity.x.array[dofs] = p_init
            print(f"  [INFO] Subdomain '{name}': set initial porosity to {p_init:.4f} ({len(dofs)} DOFs)")

        self.porosity.x.scatter_forward()
        self.porosity_n.x.array[:] = self.porosity.x.array
        self.porosity_n.x.scatter_forward()

    def update_porosity_dependent_properties(self, T_eval, p_eval):
        """
        Dynamically evaluate and update porosity-dependent material properties.
        This modifies mat["k"] in place so the staggered conduction step sees the
        updated field.
        """
        # k, T and p must be sampled at the same DOFs to combine pointwise. T and
        # k live on V_t (Lagrange-1). With the CG porosity discretisation p is
        # also Lagrange-1, so the arrays align directly. With the DG one p lives
        # on V_p (DG-1) and has a different layout/size, so first interpolate it
        # onto a cached Lagrange-1 helper before the Maxwell-Eucken arithmetic.
        T_vals = T_eval.x.array
        if p_eval.x.array.shape == T_vals.shape:
            p_vals_aligned = p_eval.x.array
        else:
            if getattr(self, "_p_on_Vt", None) is None:
                self._p_on_Vt = dolfinx.fem.Function(self.V_t)
            self._p_on_Vt.interpolate(p_eval)
            p_vals_aligned = self._p_on_Vt.x.array

        for name, mat in self.materials.items():
            if mat.get("thermal_conductivity_model") == "kato_porosity":
                p_vals = p_vals_aligned

                x = float(mat.get("stoichiometry_deviation", 0.025))
                k_He = float(mat.get("helium_conductivity", 0.69))

                # Kato temperature dependent conductivity k(T)
                A = 2.713 * x + 1.595e-2
                B = (2.493 - 2.625 * x) * 1.0e-4
                C = 1.541e11
                D = 1.522e4

                T_clamped = np.clip(T_vals, 298.0, 4000.0)
                k_T = 1.0 / (A + B * T_clamped) + (C / (T_clamped**2.5)) * np.exp(-D / T_clamped)

                # Maxwell-Eucken correction for porosity p
                p_clamped = np.clip(p_vals, 0.0, 1.0)
                num = k_He + 2.0 * k_T - 2.0 * p_clamped * (k_T - k_He)
                den = k_He + 2.0 * k_T + p_clamped * (k_T - k_He)

                k_Tp = k_T * (num / den)

                mat["k"].x.array[:] = k_Tp
                mat["k"].x.scatter_forward()

    def _apply_vertex_limiter(self, p_func):
        """
        Vertex-based (Kuzmin 2010, J. Comput. Appl. Math. 233, 3077) slope
        limiter for P1 DG porosity. Each cell's solution is written as
        mean + slope about the centroid; the slope is scaled by a single factor
        alpha_K in [0, 1], the largest value keeping every vertex value within
        the range of the cell means of the cells sharing that vertex.

        Two guarantees, both provable from the construction above:
          * cell mean (hence ∫p, the void volume) is preserved exactly for any
            alpha_K -> conservative, unlike np.clip;
          * no new extrema (discrete maximum principle): if all cell means are
            in [0, 1] the limited field is in [0, 1] -> bounded without a clamp.

        Operates in place on the owned + ghost DG-1 dofs and scatters forward
        (owners win on the partition boundary; a one-deep ghost layer gives the
        boundary vertices the correct neighbouring cell means).
        """
        mesh = self.mesh
        tdim = mesh.topology.dim
        mesh.topology.create_connectivity(tdim, 0)
        c2v = mesh.topology.connectivity(tdim, 0)
        cv = c2v.array.reshape(-1, tdim + 1)          # (num_cells, nverts) local vertex ids
        num_cells = cv.shape[0]

        dofs = np.asarray(self.V_p.dofmap.list).reshape(num_cells, -1)  # P1 DG: dof k ↔ vertex k
        arr = p_func.x.array
        cell_vals = arr[dofs]                         # (num_cells, ndofs) vertex values
        pbar = cell_vals.mean(axis=1)                 # P1 simplex: centroid value = vertex mean = cell mean

        nverts = int(cv.max()) + 1 if cv.size else 0
        vmax = np.full(nverts, -np.inf)
        vmin = np.full(nverts, np.inf)
        flat_v = cv.ravel()
        flat_pbar = np.repeat(pbar, cv.shape[1])      # cell-major, matches cv.ravel()
        np.maximum.at(vmax, flat_v, flat_pbar)        # vertex <- max of incident cell means
        np.minimum.at(vmin, flat_v, flat_pbar)        # vertex <- min of incident cell means

        diff = cell_vals - pbar[:, None]              # p_i - pbar (the slope contribution)
        hi = vmax[cv] - pbar[:, None]                 # positive headroom to the local max
        lo = vmin[cv] - pbar[:, None]                 # negative headroom to the local min
        tiny = 1e-14
        ratio = np.ones_like(diff)
        pos = diff > tiny
        neg = diff < -tiny
        ratio[pos] = np.minimum(1.0, hi[pos] / diff[pos])
        ratio[neg] = np.minimum(1.0, lo[neg] / diff[neg])
        np.clip(ratio, 0.0, 1.0, out=ratio)
        alpha = ratio.min(axis=1)                     # (num_cells,) one factor per cell

        arr[dofs] = pbar[:, None] + alpha[:, None] * diff
        p_func.x.scatter_forward()

    def _build_saturation_cache(self):
        """Cache for the conservative saturation cap: cell→cell (face-neighbour)
        adjacency (CSR `nbr`/`ptr`), DG-1 cell dofs, cell areas |K|, cell centroids,
        and the per-edge centroid offset `d = c_K − c_neighbour` (aligned with
        `nbr`). Serial / owned-cell scope (the pore-migration case runs serial)."""
        mesh = self.mesh
        V = self.V_p
        tdim = mesh.topology.dim
        mesh.topology.create_connectivity(tdim, tdim - 1)
        mesh.topology.create_connectivity(tdim - 1, tdim)
        c2f = mesh.topology.connectivity(tdim, tdim - 1)
        f2c = mesh.topology.connectivity(tdim - 1, tdim)
        nc = mesh.topology.index_map(tdim).size_local

        dofs = np.asarray(V.dofmap.list).reshape(-1, tdim + 1)[:nc]
        cen = V.tabulate_dof_coordinates()[:, :2][dofs].mean(axis=1)   # (nc, 2) centroids

        ptr = np.zeros(nc + 1, dtype=np.int64)
        flat, dvec = [], []
        for K in range(nc):
            s = sorted({int(c) for f in c2f.links(K) for c in f2c.links(f)
                        if c != K and c < nc})
            ptr[K + 1] = ptr[K] + len(s)
            flat.extend(s)
            for j in s:
                dvec.append(cen[K] - cen[j])
        nbr = np.asarray(flat, dtype=np.int64)
        nbr_d = np.asarray(dvec).reshape(-1, 2) if dvec else np.zeros((0, 2))

        Mv = dolfinx.fem.assemble_vector(dolfinx.fem.form(ufl.TestFunction(V) * ufl.dx))
        Mv.scatter_reverse(dolfinx.la.InsertMode.add)
        area = Mv.array[dofs].sum(axis=1)              # |K| = Σ_{i∈K} ∫φ_i (DG dofs are cell-exclusive)
        self._sat_cache = (nbr, ptr, area, dofs, nbr_d)

    def _apply_saturation_cap(self, p_func, v_vec=None):
        """Conservative saturation cap: enforce cell-mean porosity ≤ 1 (a void
        cannot exceed unit porosity) without discarding mass. The excess mass of an
        over-saturated cell is redistributed to its neighbours in proportion to
        their remaining capacity (1 − p̄), swept until no cell mean exceeds 1.
        Capacity is largest in the emptier (outer) cells, so the excess flows
        outward and the saturated void core grows. Total ∫p is preserved to
        round-off; the subsequent vertex limiter keeps the DG dofs in [0,1].
        Enabled by porosity.saturation_cap."""
        if getattr(self, "_sat_cache", None) is None:
            self._build_saturation_cache()
        nbr, ptr, area, dofs, nbr_d = self._sat_cache
        nc = area.shape[0]
        arr = p_func.x.array
        cell_vals = arr[dofs]
        old_pbar = cell_vals.mean(axis=1)
        mass = old_pbar * area
        tol = 1e-13
        total0 = mass.sum()
        sweeps = int(self.porosity_cfg.get("saturation_sweeps", 200))

        # Per-cell velocity (DG0) to direct the excess upstream (against v, i.e.
        # outward): a saturated cell pushes its excess back toward the cells the
        # pores came from, so the void front advances against the pore flow and
        # the redistribution converges in O(front-depth) sweeps, not diffusively.
        v_cell = None
        if v_vec is not None:
            gdim = self.mesh.geometry.dim
            V0 = dolfinx.fem.functionspace(self.mesh, ("DG", 0, (gdim,)))
            vf = dolfinx.fem.Function(V0)
            vf.interpolate(dolfinx.fem.Expression(v_vec, V0.element.interpolation_points))
            v_cell = vf.x.array.reshape(-1, gdim)[:nc, :2]   # planar (xy) components

        capped = False
        for _ in range(sweeps):
            pbar = mass / area
            over = pbar > 1.0 + 1e-12
            if not np.any(over):
                break
            capped = True
            excess = np.maximum(0.0, pbar - 1.0) * area        # excess mass per cell
            mass = np.minimum(mass, area)                       # cap means at 1
            recv = np.zeros(nc)
            cap = np.maximum(0.0, 1.0 - mass / area) * area     # capacity mass per cell
            for K in np.nonzero(over)[0]:
                a, b = ptr[K], ptr[K + 1]
                js = nbr[a:b]
                if js.size == 0:
                    recv[K] += excess[K]
                    continue
                # 1) upstream weight: w_j = max(0, v_K · (c_K − c_j))  (neighbour j
                #    is upstream when v points from j into K)
                w = None
                if v_cell is not None:
                    w = np.maximum(0.0, nbr_d[a:b] @ v_cell[K])
                    if w.sum() <= tol:
                        w = None
                # 2) else prefer neighbours with capacity; 3) else pass through equally
                if w is None:
                    w = cap[js]
                    if w.sum() <= tol:
                        w = np.ones(js.size)
                recv[js] += excess[K] * (w / w.sum())
            mass = mass + recv
        if not capped:
            return
        pbar = mass / area
        # conservation guard: tiny round-off in the sweep is corrected by a global
        # rescale of the change so Σ mass is exactly preserved.
        drift = mass.sum() - total0
        if abs(drift) > tol * max(total0, 1.0):
            mass -= drift / nc
            pbar = mass / area
        arr[dofs] = cell_vals + (pbar - old_pbar)[:, None]      # uniform shift preserves slope
        p_func.x.scatter_forward()
        self._apply_vertex_limiter(p_func)                       # DMP -> dofs back into [0,1]

    def _ssprk3_advance(self, p_new, p_n, dt, v_vec):
        """
        Explicit 3-stage SSP-RK3 (Shu–Osher) advance of the DG pore-advection
        equation over [t^n, t^n + dt], starting from p_n, with internal
        sub-stepping to the advective CFL and the vertex limiter applied after
        every stage (SSP needs the limiter at every stage, not just at the end).

        Semi-discrete form:  M dp/dt = R(p),  R(p) = -S(p) + g, where S collects
        the upwind facet flux, the -∫ p v·∇w volume term, the outflow boundary
        term and the optional SIPG block, and g is the inflow boundary vector.
        M is row-sum-lumped: for DG-1, M_i = ∫ φ_i dx (partition of unity), which
        is diagonal so each forward-Euler stage is a single sparse assemble + a
        pointwise divide — no linear solve. Mass (Σ_i M_i p_i = ∫p) is conserved
        to round-off because 1ᵀR = -(boundary efflux) and the limiter preserves
        cell means.

        v_vec is passed in (UFL) so the same routine serves a unit test with a
        prescribed closed velocity field.
        """
        import math
        mesh = self.mesh
        V = self.V_p
        v_test = ufl.TestFunction(V)
        n_vec = ufl.FacetNormal(mesh)
        h_e = ufl.CellDiameter(mesh)

        # --- CFL sub-step count from per-cell h/|v|_max ---
        # Use the per-cell maximum |v| (over the cell's DG1 vertices), not the
        # centroid value: the pore velocity peaks where ∇T is steepest, so a
        # centroid (DG0) estimate understates the peak and the explicit step goes
        # unstable (negative p). A safety factor dg_cfl_safety adds margin.
        tdim = mesh.topology.dim
        ncl = mesh.topology.index_map(tdim).size_local
        cells = np.arange(ncl, dtype=np.int32)
        h_cells = mesh.h(tdim, cells)
        V1 = dolfinx.fem.functionspace(mesh, ("DG", 1))
        vmag1 = dolfinx.fem.Function(V1)
        vmag1.interpolate(dolfinx.fem.Expression(
            ufl.sqrt(ufl.dot(v_vec, v_vec) + 1e-30),
            V1.element.interpolation_points))
        d1 = np.asarray(V1.dofmap.list).reshape(-1, tdim + 1)[:ncl]
        vmag_cellmax = vmag1.x.array[d1].max(axis=1)
        ratios = h_cells / (vmag_cellmax + 1e-30)
        local_min = float(ratios.min()) if ratios.size else np.inf
        C_cfl = float(self.porosity_cfg.get("dg_cfl", 1.0 / 3.0))
        safety = float(self.porosity_cfg.get("dg_cfl_safety", 0.8))
        dt_cfl = C_cfl * safety * mesh.comm.allreduce(local_min, op=MPI.MIN)
        if dt_cfl > 0.0 and np.isfinite(dt_cfl):
            N_sub = max(1, int(math.ceil(dt / dt_cfl)))
        else:
            N_sub = 1
        cap = int(self.porosity_cfg.get("dg_max_substeps", 5000))
        if N_sub > cap:
            print(f"  SSP-RK3: CFL wants {N_sub} substeps, capped at {cap} "
                  f"(limiter still bounds; check dt/mesh)")
            N_sub = cap
        dt_sub = dt / N_sub

        # --- residual 1-form R(p_new, w); p_new is the live stage coefficient ---
        vn = ufl.dot(ufl.avg(v_vec), n_vec("+"))
        p_up = ufl.conditional(vn > 0, p_new("+"), p_new("-"))
        vn_b = ufl.dot(v_vec, n_vec)
        rim_inflow = self.porosity_cfg.get("rim_inflow_porosity", None)
        p_in = (dolfinx.fem.Constant(mesh, PETSc.ScalarType(float(rim_inflow)))
                if rim_inflow is not None else p_n)
        R_form = (
            p_new * ufl.dot(v_vec, ufl.grad(v_test)) * ufl.dx
            - vn * p_up * ufl.jump(v_test) * ufl.dS
            - ufl.conditional(vn_b > 0, vn_b * p_new * v_test, 0.0) * ufl.ds
            - ufl.conditional(vn_b < 0, vn_b * p_in * v_test, 0.0) * ufl.ds
        )
        D_p = float(self.porosity_cfg.get("diffusion", 0.0))
        if D_p > 0.0:
            D_const = dolfinx.fem.Constant(mesh, PETSc.ScalarType(D_p))
            gamma = float(self.porosity_cfg.get("sipg_penalty", 10.0))
            h_avg = (h_e("+") + h_e("-")) / 2.0
            R_form -= D_const * ufl.dot(ufl.grad(p_new), ufl.grad(v_test)) * ufl.dx
            R_form += D_const * ufl.dot(ufl.avg(ufl.grad(p_new)), ufl.jump(v_test, n_vec)) * ufl.dS
            R_form += D_const * ufl.dot(ufl.avg(ufl.grad(v_test)), ufl.jump(p_new, n_vec)) * ufl.dS
            R_form -= D_const * (gamma / h_avg) * ufl.jump(p_new) * ufl.jump(v_test) * ufl.dS
        R_compiled = dolfinx.fem.form(R_form)

        # --- lumped mass M_i = ∫ φ_i dx (positive, diagonal) ---
        M_vec = dolfinx.fem.assemble_vector(dolfinx.fem.form(v_test * ufl.dx))
        M_vec.scatter_reverse(dolfinx.la.InsertMode.add)
        M = M_vec.array
        no = V.dofmap.index_map.size_local

        def residual():
            p_new.x.scatter_forward()
            b = dolfinx.fem.assemble_vector(R_compiled)
            b.scatter_reverse(dolfinx.la.InsertMode.add)
            return b.array

        limiter = str(self.porosity_cfg.get("dg_limiter", "vertex")).lower()

        def stage(a0, a1, c0, c1, cR, R):
            # owned forward-Euler update, then the per-stage limiter (SSP needs it
            # every stage). "vertex" preserves cell means -> conservative.
            p_new.x.array[:no] = (c0 * a0[:no] + c1 * a1[:no]
                                  + cR * dt_sub * (R[:no] / M[:no]))
            if limiter == "vertex":
                self._apply_vertex_limiter(p_new)
            elif limiter == "clamp":
                p_new.x.array[:] = np.clip(p_new.x.array, 0.0, 1.0)
                p_new.x.scatter_forward()
            else:
                p_new.x.scatter_forward()

        # start from p_n (each staggered call recomputes from the fixed t^n level)
        p_new.x.array[:] = p_n.x.array
        p_new.x.scatter_forward()
        for _ in range(N_sub):
            p0 = p_new.x.array.copy()
            R = residual(); stage(p0, p0, 1.0, 0.0, 1.0, R)        # u1 = p0 + dt L
            u1 = p_new.x.array.copy()
            R = residual(); stage(p0, u1, 0.75, 0.25, 0.25, R)     # u2
            u2 = p_new.x.array.copy()
            R = residual(); stage(p0, u2, 1.0/3.0, 2.0/3.0, 2.0/3.0, R)  # u3
        return N_sub

    def _porosity_step(self, p_new, p_n, dt, T_current, stag_tol_p, prev_res_p):
        """
        Solve one pore-advection sub-step (Barani et al. 2022, Eq. 1b) with SU
        stabilisation, inside the staggered T-p fixed-point loop.

        p_new : current porosity iterate, updated in place.
        p_n   : porosity at the start of the time step (t^n). Held fixed here and
                used only for the backward-Euler time term — it must not be
                overwritten, otherwise each staggered iteration would advect a
                full dt forward from the previous iterate.
        """
        if dt <= 0.0:
            print("  Porosity step: dt <= 0.0 -> preserving porosity (converged: True)")
            return True, 0.0, 0.0, prev_res_p

        # Previous staggered iterate: reference for the Eq. (2) convergence check
        # and for under-relaxation. Distinct from p_n (the t^n time level).
        p_prev_iter = p_new.x.array.copy()

        u_p, v_test = ufl.TrialFunction(self.V_p), ufl.TestFunction(self.V_p)
        dt_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(dt))

        R = 8.314  # universal gas constant [J/(mol K)]
        c1, c2, c3, c4 = self.c1, self.c2, self.c3, self.c4
        Hs = self.H_s

        # Pore velocity (Novascone 2018 Eq. 4 / Barani 2022 Eq. 9), SI. The lumped
        # prefactor v0 bundles Sens' c0, the heat-of-vaporisation prefactor and the
        # vapour-pressure pre-exponential; the Arrhenius term uses Hs [J/mol] and R
        # [SI]. v points up the thermal gradient (inward, toward the hot centre).
        poly = c1 + c2 * T_current + c3 * T_current**2 + c4 * T_current**3
        mobility = self.v0 * poly * (T_current**(-2.5)) * ufl.exp(-Hs / (R * T_current))
        v_vec = mobility * ufl.grad(T_current)

        # Stabilisation for the advection-dominated pore transport, selected by
        # porosity.stabilisation:
        #  - "su"   (default): streamline-upwind artificial diffusion
        #            K = (h_e/2|v|) v⊗v (Barani Eq. 3-4), added only to the
        #            bilinear form. Robust but inconsistent — it also diffuses the
        #            converged solution, smearing the void front.
        #  - "supg": streamline-upwind Petrov-Galerkin. The test function is
        #            perturbed by tau (v·∇w) and weighted against the full strong
        #            residual p/dt + ∇·(v p) - p_n/dt, so the added diffusion
        #            vanishes at the exact solution (consistent -> sharper front).
        #  - "dg":  upwind discontinuous Galerkin (porosity.discretisation: dg),
        #            the same operator family as cluster dynamics. The upwind
        #            facet flux is the stabilisation (no SU/SUPG, no mass-matrix
        #            surgery); an optional SIPG block adds physical diffusion.
        h_e = ufl.CellDiameter(self.mesh)
        n_vec = ufl.FacetNormal(self.mesh)
        disc = str(self.porosity_cfg.get("discretisation", "cg")).lower()

        dg_explicit = False
        if disc == "dg":
            # Time integrator: "ssprk3" (default,
            # explicit, bounded, sidesteps the staggered stiffness) or "be"
            # (implicit Backward-Euler, the assembled system below).
            integrator = str(self.porosity_cfg.get("dg_integrator", "ssprk3")).lower()
            if integrator == "ssprk3":
                n_sub = self._ssprk3_advance(p_new, p_n, dt, v_vec)
                print(f"  Porosity step: SSP-RK3 explicit, {n_sub} CFL substep(s)")
                dg_explicit = True

        if disc == "dg" and not dg_explicit:
            # Conservative form, integrated by parts cell-wise:
            #   ∫ (p/dt) w dx - ∫ p v·∇w dx + Σ_facets <flux,[w]> = ∫ (p_n/dt) w dx
            a_p = (u_p / dt_const) * v_test * ufl.dx
            a_p -= u_p * ufl.dot(v_vec, ufl.grad(v_test)) * ufl.dx
            L_p = (p_n / dt_const) * v_test * ufl.dx

            # Interior facets: pure upwinding on the (discontinuous, since T is
            # P1) normal velocity vn = avg(v)·n⁺. The upwind trace is taken from
            # the cell the flow comes from; jump(w) = w⁺ - w⁻.
            vn = ufl.dot(ufl.avg(v_vec), n_vec("+"))
            p_up = ufl.conditional(vn > 0, u_p("+"), u_p("-"))
            a_p += vn * p_up * ufl.jump(v_test) * ufl.dS

            # Exterior facets: outflow (v·n>0) carries the interior value (kept in
            # the bilinear form); inflow (v·n<0) carries an upstream value, moved
            # to the RHS. Default inflow value is p_n on the boundary, so the cold
            # rim keeps its fabrication porosity (matching the CG natural-BC
            # intent); porosity.rim_inflow_porosity overrides it with a constant.
            vn_b = ufl.dot(v_vec, n_vec)
            rim_inflow = self.porosity_cfg.get("rim_inflow_porosity", None)
            p_in = (
                dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(float(rim_inflow)))
                if rim_inflow is not None else p_n
            )
            a_p += ufl.conditional(vn_b > 0, vn_b * u_p * v_test, 0.0) * ufl.ds
            L_p -= ufl.conditional(vn_b < 0, vn_b * p_in * v_test, 0.0) * ufl.ds

            # Optional SIPG diffusion (off by default -> pure advection, since the
            # micrometric pore diffusivity is negligible in Barani's model).
            D_p = float(self.porosity_cfg.get("diffusion", 0.0))
            if D_p > 0.0:
                D_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(D_p))
                gamma = float(self.porosity_cfg.get("sipg_penalty", 10.0))
                h_avg = (h_e("+") + h_e("-")) / 2.0
                a_p += D_const * ufl.dot(ufl.grad(u_p), ufl.grad(v_test)) * ufl.dx
                a_p -= D_const * ufl.dot(ufl.avg(ufl.grad(u_p)), ufl.jump(v_test, n_vec)) * ufl.dS
                a_p -= D_const * ufl.dot(ufl.avg(ufl.grad(v_test)), ufl.jump(u_p, n_vec)) * ufl.dS
                a_p += D_const * (gamma / h_avg) * ufl.jump(u_p) * ufl.jump(v_test) * ufl.dS
        elif disc != "dg":
            v_mag = ufl.sqrt(ufl.dot(v_vec, v_vec) + 1e-30)
            stab = str(self.porosity_cfg.get("stabilisation", "su")).lower()

            # Conservative advection integrated by parts:
            #   ∫ ∇·(v p) w dx = -∫ p v·∇w dx + ∫_∂Ω p (v·n) w ds
            # Default: drop the boundary term (homogeneous natural BC, as in Barani
            # Sec. 4.1 for SU). v·n = 0 on the symmetry faces and is negligible at
            # the cold outer rim, so the rim keeps its fabrication porosity without
            # an imposed value. Set porosity.rim_inflow_porosity to force a value.
            a_p = (u_p / dt_const) * v_test * ufl.dx
            a_p -= u_p * ufl.dot(v_vec, ufl.grad(v_test)) * ufl.dx

            L_p = (p_n / dt_const) * v_test * ufl.dx

            if stab == "supg":
                # SUPG parameter for transient advection: tends to dt/2 for small
                # dt (time-step limited) and to the cell transit time h_e/(2|v|)
                # for large dt / fast flow. For P1 the second derivative in
                # div(v u_p) = u_p div(v) + v·grad(u_p) drops, so the strong
                # residual is exact on the element interiors.
                tau = 1.0 / ufl.sqrt((2.0 / dt_const) ** 2 + (2.0 * v_mag / h_e) ** 2)
                supg_w = tau * ufl.dot(v_vec, ufl.grad(v_test))
                a_p += (u_p / dt_const + ufl.div(v_vec * u_p)) * supg_w * ufl.dx
                L_p += (p_n / dt_const) * supg_w * ufl.dx
            else:
                K_term = (h_e / (2.0 * v_mag)) * ufl.dot(v_vec, ufl.grad(u_p)) * ufl.dot(v_vec, ufl.grad(v_test))
                a_p += K_term * ufl.dx

            rim_inflow = self.porosity_cfg.get("rim_inflow_porosity", None)
            if rim_inflow is not None:
                rim_label = self.porosity_cfg.get("rim_label", "outer")
                if rim_label not in self.label_map:
                    raise KeyError(
                        f"porosity.rim_inflow_porosity: rim label '{rim_label}' is "
                        "not in the geometry label map — set porosity.rim_label to "
                        "the rim facet label."
                    )
                ds_outer = self.ds_tags[self.label_map[rim_label]]
                p_inflow_const = dolfinx.fem.Constant(self.mesh, PETSc.ScalarType(float(rim_inflow)))
                # Standard weak inflow/outflow split: the prescribed value enters
                # only where v·n < 0 (inflow); where the rim is locally outflow the
                # unknown leaves through the LHS term.
                v_n = ufl.dot(v_vec, n_vec)
                v_n_in = ufl.conditional(ufl.lt(v_n, 0.0), v_n, 0.0)
                v_n_out = ufl.conditional(ufl.gt(v_n, 0.0), v_n, 0.0)
                a_p += v_n_out * u_p * v_test * ds_outer
                L_p -= p_inflow_const * v_n_in * v_test * ds_outer

        if not dg_explicit:
            petsc_opts = self.get_solver_options(
                solver_type=self.porosity_cfg.get("linear_solver", "direct_mumps"),
                physics="porosity",
                rtol=self.porosity_cfg.get("rtol", 1.0e-8),
            )
            problem_p = dolfinx.fem.petsc.LinearProblem(
                a_p, L_p, bcs=[], u=p_new,
                petsc_options=petsc_opts,
                petsc_options_prefix="porosity_",
            )
            problem_p.solve()

        # Enforce the physical bounds [0, 1]. The CG path can only project
        # pointwise (np.clip), which is non-conservative. The DG path defaults to
        # the vertex-based limiter (porosity.dg_limiter), which is bounded and
        # mass-conserving (it rescales the slope, leaving the cell mean and hence
        # ∫p untouched); "clamp" clips pointwise, "none" disables it.
        # SSP-RK3 already limited every stage, so skip the post-step limiter there.
        if disc == "dg":
            if not dg_explicit:
                limiter = str(self.porosity_cfg.get("dg_limiter", "vertex")).lower()
                if limiter == "vertex":
                    self._apply_vertex_limiter(p_new)
                elif limiter == "clamp":
                    p_new.x.array[:] = np.clip(p_new.x.array, 0.0, 1.0)
                    p_new.x.scatter_forward()
                # "none": leave the raw DG solve untouched (upwind keeps it bounded)
            # Conservative saturation cap (porosity.saturation_cap): enforce p ≤ 1
            # by redistributing the over-saturated excess outward, not clipping it.
            # Applied after the limiter (which already gave positivity/DMP). Runs
            # for both ssprk3 and be on the DG path.
            if bool(self.porosity_cfg.get("saturation_cap", False)):
                self._apply_saturation_cap(p_new, v_vec=v_vec)
        else:
            p_new.x.array[:] = np.clip(p_new.x.array, 0.0, 1.0)
            p_new.x.scatter_forward()

        # Under-relaxation against the previous staggered iterate. Two options:
        #  - fixed factor porosity.relax in (0, 1];
        #  - Aitken Δ² (porosity.aitken: true): the factor is recomputed each
        #    staggered iteration from the last two raw residuals r_k = p_raw - p_k,
        #    ω_{k+1} = -ω_k (r_{k-1}·Δr)/|Δr|², Δr = r_k - r_{k-1}, clamped to a
        #    sane band. p_raw is the clamped solve output, so ω in [0.05, 1] keeps
        #    p_new a convex combination of two values in [0, 1] — no re-clamp.
        #    State (_aitken_p_*) is reset per time step by solve_staggered and is
        #    kept separate from the displacement Aitken (_aitken_*).
        if bool(self.porosity_cfg.get("aitken", False)):
            r_k = p_new.x.array - p_prev_iter
            r_prev = getattr(self, "_aitken_p_R_prev", None)
            omega0 = float(self.porosity_cfg.get("aitken_omega0", 0.5))
            # Unlike the displacement loop, porosity restarts from omega0 on the
            # first staggered iteration of a step rather than carrying omega over.
            usable = r_prev is not None and r_prev.shape == r_k.shape
            omega = aitken_omega(
                r_k, r_prev,
                float(getattr(self, "_aitken_p_omega", omega0)) if usable else omega0,
                self.mesh.comm, self.V_p.dofmap.index_map.size_local, 0.05, 1.0)
            self._aitken_p_R_prev = r_k.copy()
            self._aitken_p_omega = omega
            p_new.x.array[:] = p_prev_iter + omega * r_k
            p_new.x.scatter_forward()
            print(f"  Porosity step: Aitken omega = {omega:.4f}")
        else:
            relax_p = float(self.porosity_cfg.get("relax", 1.0))
            if relax_p < 1.0:
                p_new.x.array[:] = relax_p * p_new.x.array + (1.0 - relax_p) * p_prev_iter
                p_new.x.scatter_forward()

        # Convergence between successive staggered iterates. Two metrics
        # (porosity.conv_metric):
        #  - "max_dof" (default): mixed rel/abs max over DOFs (Barani Eq. 2).
        #  - "integral": relative change of the conserved void volume
        #    Q = ∫p between iterates, |Q_k - Q_{k-1}|/Q < stag_tol_rel. At
        #    saturation the front DOFs flicker (max_dof never converges) while
        #    the total void volume is steady — the physically meaningful measure
        #    of staggered convergence: dQ/dt = -(boundary efflux) holds exactly for
        #    this scheme, so a steady Q means a converged transport step.
        eps_rel = float(self.porosity_cfg.get("stag_tol_rel", 1.0e-6))
        eps_abs = float(self.porosity_cfg.get("stag_tol_abs", 1.0e-8))
        n_owned = self.V_p.dofmap.index_map.size_local
        p_now = p_new.x.array[:n_owned]
        diff = np.abs(p_now - p_prev_iter[:n_owned])
        check_vals = diff - np.abs(p_now) * eps_rel - eps_abs
        local_max = float(np.max(check_vals)) if check_vals.size > 0 else -1.0
        max_val = self.mesh.comm.allreduce(local_max, op=MPI.MAX)

        conv_metric = str(self.porosity_cfg.get("conv_metric", "max_dof")).lower()
        if conv_metric == "integral":
            Q_now = self.mesh.comm.allreduce(
                dolfinx.fem.assemble_scalar(dolfinx.fem.form(p_new * ufl.dx)), op=MPI.SUM)
            p_prev_fn = dolfinx.fem.Function(self.V_p)
            p_prev_fn.x.array[:] = p_prev_iter
            Q_prev = self.mesh.comm.allreduce(
                dolfinx.fem.assemble_scalar(dolfinx.fem.form(p_prev_fn * ufl.dx)), op=MPI.SUM)
            dQ = abs(Q_now - Q_prev) / max(abs(Q_now), 1e-300)
            # Dedicated, physically-scaled tolerance: the void volume is converged
            # when it changes < conv_integral_tol (default 1e-4 = 0.01%) between
            # staggered iterates. Looser than the DOF tolerance on purpose — at
            # saturation the front flickers at the dof level (~1e-2) while the
            # total void volume is steady to ~1e-5.
            q_tol = float(self.porosity_cfg.get("conv_integral_tol", 1.0e-4))
            conv_p = dQ < q_tol
            print(f"  Porosity step: |dQ|/Q = {dQ:.3e} (tol {q_tol:.1e}, converged: {conv_p}) "
                  f"[max_dof residual = {max_val:.3e}]")
        else:
            conv_p = max_val < 0.0
            print(f"  Porosity step: max(diff - rel - abs) = {max_val:.3e} (converged: {conv_p})")

        print(f"  p_new: min={p_new.x.array.min():.4f}, max={p_new.x.array.max():.4f}, mean={p_new.x.array.mean():.4f}")

        return conv_p, 0.0, 0.0, prev_res_p
