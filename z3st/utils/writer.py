# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
OutputWriter: unified VTU / XDMF writer for z3st.

Both backends share the same field set (temperature, displacement, per-material
stress / von Mises / hydrostatic / strain-energy density, global strain,
per-material heat flux, damage + crack-driving-force, cluster density,
cumulative plastic strain).

The class pre-allocates all FE function spaces, Function objects, and
``dolfinx.fem.Expression`` objects in ``__init__``. Per-step ``write(t, step)``
does only interpolation + I/O, so the UFL JIT-compile that used to happen
inside the time loop in ``__main__.py`` now happens once at setup time.

Naming conventions (preserved from the previous code path):
  - VTU, ``n_steps == 1``                : ``output/<basename>.vtu``
  - VTU, ``n_steps > 1``                 : ``output/<basename>_NNNN.vtu``
  - XDMF                                 : single ``output/<basename>.xdmf`` + ``.h5``

A ``vtkhdf`` backend was prototyped against dolfinx 0.10 but reverted: the
0.10 ``dolfinx.io.vtkhdf`` functional API
(``write_mesh / write_point_data / write_cell_data``) has no field-name
parameter, so the many named z3st fields cannot be written into a single
file. Revisit when dolfinx ships a class-based ``VTKHDFFile`` with per-field
naming.
"""

import os

import dolfinx
import numpy as np
import pyvista
import ufl
from dolfinx.io import XDMFFile


# ---------------------------------------------------------------------------
# Derived-field helpers (numpy on host arrays after interpolation)
# ---------------------------------------------------------------------------

def _von_mises(sig):
    """Von Mises equivalent stress for an (..., 3, 3) batch.

    In the axisymmetric regime the hoop strain eps_tt = u_r / r is singular at
    the axis r=0, so the stress carries NaN at axis nodes. Suppress the
    resulting 'invalid value' warning and return a finite field (0 at the
    singular nodes) rather than poisoning the exported VonMises field with NaN.
    """
    s_xx, s_yy, s_zz = sig[..., 0, 0], sig[..., 1, 1], sig[..., 2, 2]
    s_xy, s_yz, s_zx = sig[..., 0, 1], sig[..., 1, 2], sig[..., 2, 0]
    with np.errstate(invalid="ignore"):
        vm = np.sqrt(
            0.5 * ((s_xx - s_yy) ** 2 + (s_yy - s_zz) ** 2 + (s_zz - s_xx) ** 2)
            + 3.0 * (s_xy ** 2 + s_yz ** 2 + s_zx ** 2)
        )
    return np.nan_to_num(vm, nan=0.0, posinf=0.0, neginf=0.0)


def _hydrostatic(sig):
    """Hydrostatic pressure p = tr(sigma) / 3 on an (..., 3, 3) batch. Guarded
    for the same axisymmetric r=0 NaN as :func:`_von_mises`."""
    p = (sig[..., 0, 0] + sig[..., 1, 1] + sig[..., 2, 2]) / 3.0
    return np.nan_to_num(p, nan=0.0, posinf=0.0, neginf=0.0)


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

class OutputWriter:
    """Unified VTU / XDMF writer with one-shot setup.

    Parameters
    ----------
    problem : Spine
        The z3st problem. Must have already been initialised (``initialize_fields``,
        ``set_boundary_conditions``) and have ``get_results()`` populating the
        symbolic stress / strain / energy-density UFL expressions.
    output_format : {"vtu", "xdmf"}
    output_dir : str
        Directory to write into (created if absent).
    filename : str | None
        Base filename. Defaults to ``fields.vtu`` or ``fields.xdmf`` based on
        ``output_format``.
    n_steps : int
        Total number of time steps. Controls whether VTU files carry a step
        suffix (``n_steps > 1``) or use the plain basename (``n_steps == 1``).
    """

    _DEFAULT_EXT = {"vtu": ".vtu", "xdmf": ".xdmf"}

    def __init__(self, problem, output_format="vtu", output_dir="output",
                 filename=None, n_steps=1):
        fmt = output_format.lower()
        if fmt not in ("vtu", "xdmf"):
            raise ValueError(
                f"Unsupported output format '{output_format}'; expected 'vtu' or 'xdmf'. "
                "('vtkhdf' was prototyped but reverted — see module docstring.)"
            )

        # Fallback
        if fmt == "vtu" and problem.mesh.comm.size > 1:
            if problem.mesh.comm.rank == 0:
                print(f"[WARNING] VTU output is serial-only; running on "
                      f"{problem.mesh.comm.size} MPI ranks -> switching to XDMF "
                      "to avoid a corrupted file. Set 'format: xdmf' to silence this.")
            fmt = "xdmf"

        self.problem = problem
        self.format = fmt
        self.output_dir = output_dir
        self.n_steps = max(1, int(n_steps))
        os.makedirs(output_dir, exist_ok=True)

        # Normalise to the format's extension
        ext = self._DEFAULT_EXT[fmt]
        if filename is None:
            filename = "fields"
        filename = os.path.splitext(filename)[0] + ext
        self.filename = filename

        # The UFL expressions pre-compiled below live on problem.strain,
        # problem.stress[name], etc., populated by get_results(). Call it now
        # if the user hasn't.
        if not hasattr(problem, "stress") or not hasattr(problem, "strain"):
            problem.get_results()

        self.mesh = problem.mesh
        self.tdim = problem.tdim

        # Quadrature metadata for plasticity-aware projections.
        self._metadata = {}
        if getattr(problem, "q_degree", None) is not None:
            self._metadata = {
                "quadrature_degree": problem.q_degree,
                "quadrature_scheme": "default",
            }

        # FE spaces (always 3 for tensor padding, even in 2D / axisymmetric).
        dim = 3
        self.V_tensor_cells = dolfinx.fem.functionspace(self.mesh, ("DG", 0, (dim, dim)))
        self.V_vector_cells = dolfinx.fem.functionspace(self.mesh, ("DG", 0, (dim,)))
        self.V_scalar_cells = dolfinx.fem.functionspace(self.mesh, ("DG", 0))
        self.V_tensor_points = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1, (dim, dim)))
        self.V_scalar_points = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))

        # Cluster: DG1 → CG1 interpolation for XDMF visualisation.
        self.V_c_cg = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1))

        # Pre-allocate Function targets and compile Expressions for each
        # symbolic UFL field, gated on the relevant physics being active.
        self._strain_fn_cells = None
        self._strain_fn_points = None
        self._strain_expr_cells = None
        self._strain_expr_points = None
        # Merged output: single domain-wide fields (one Stress, one HeatFlux, ...)
        # instead of one per material. Each material fills only its own cells, so
        # cell values are bit-identical to the old per-material fields; point
        # values can differ only at material interfaces (last writer wins).
        self._stress_fn_cells = None
        self._stress_fn_points = None
        self._psi_fn_cells = None
        self._psi_fn_points = None
        self._heatflux_fn = None
        self._stress_sources = []      # list of (expr_cells, expr_points, cells)
        self._psi_sources = []         # list of (expr_cells, expr_points, cells)
        self._heatflux_sources = []    # list of (expr_cells, cells)
        self._D_cell_fn = None
        self._D_cell_expr = None
        self._H_cell_fn = None
        self._H_cell_expr = None
        self._c_cg_fn = None
        self._p_cg_fn = None         # CG1 view of DG porosity for output (lazy)
        self._p_proj_problem = None  # cumulative plastic strain (lazy)
        self._u_out_fn = None        # CG1 vertex view of u for P2+ output (lazy)

        on = problem.on

        # Strain (global, mechanical only)
        if on.get("mechanical", False) and problem.strain is not None:
            self._strain_fn_cells = dolfinx.fem.Function(self.V_tensor_cells, name="Strain")
            self._strain_fn_points = dolfinx.fem.Function(self.V_tensor_points, name="Strain")
            self._strain_expr_cells = self._make_interp_or_proj(problem.strain, self.V_tensor_cells)
            self._strain_expr_points = self._make_interp_or_proj(problem.strain, self.V_tensor_points)

        # Per-material expressions feed the single merged fields. Allocate each
        # merged field once (if any material contributes), then collect every
        # material's interpolation source together with the cells it must fill.
        cell_tags = getattr(problem, "cell_tags", None) or getattr(problem, "tags", None)

        def _mat_cells(name):
            # Local cell indices of a material's subdomain, for cell-restricted
            # interpolation. None -> fill the whole mesh (single-domain case).
            if cell_tags is None:
                return None
            try:
                return cell_tags.find(problem.label_map[name])
            except Exception:
                return None

        mats = problem.materials
        mech = on.get("mechanical", False)
        therm = on.get("thermal", False) and getattr(problem, "T", None) is not None

        if mech and any(n in problem.stress for n in mats):
            self._stress_fn_cells = dolfinx.fem.Function(self.V_tensor_cells, name="Stress")
            self._stress_fn_points = dolfinx.fem.Function(self.V_tensor_points, name="Stress")
        if mech and any(n in problem.energy_density for n in mats):
            self._psi_fn_cells = dolfinx.fem.Function(self.V_scalar_cells, name="StrainEnergyDensity")
            self._psi_fn_points = dolfinx.fem.Function(self.V_scalar_points, name="StrainEnergyDensity")
        if therm:
            self._heatflux_fn = dolfinx.fem.Function(self.V_vector_cells, name="HeatFlux")

        for name, material in mats.items():
            cells = _mat_cells(name)
            if mech and name in problem.stress:
                self._stress_sources.append((
                    self._make_interp_or_proj(problem.stress[name], self.V_tensor_cells),
                    self._make_interp_or_proj(problem.stress[name], self.V_tensor_points),
                    cells,
                ))
                if name in problem.energy_density:
                    self._psi_sources.append((
                        self._make_interp_or_proj(problem.energy_density[name], self.V_scalar_cells),
                        self._make_interp_or_proj(problem.energy_density[name], self.V_scalar_points),
                        cells,
                    ))
            if therm:
                q_expr = -material["k"] * ufl.grad(problem.T)
                self._heatflux_sources.append((
                    self._make_interp_or_proj(q_expr, self.V_vector_cells),
                    cells,
                ))

        # Damage + crack-driving-force, cell projections
        if on.get("damage", False) and getattr(problem, "D", None) is not None:
            self._D_cell_fn = dolfinx.fem.Function(self.V_scalar_cells, name="Damage")
            self._D_cell_expr = self._make_interp_or_proj(problem.D, self.V_scalar_cells)
            if getattr(problem, "H", None) is not None:
                self._H_cell_fn = dolfinx.fem.Function(
                    self.V_scalar_cells, name="CrackDrivingForce"
                )
                self._H_cell_expr = self._make_interp_or_proj(problem.H, self.V_scalar_cells)

        # Cluster: pre-allocated CG1 Function for XDMF visualisation.
        if on.get("cluster", False) and getattr(problem, "c", None) is not None:
            self._c_cg_fn = dolfinx.fem.Function(self.V_c_cg, name="ClusterDensity")

        # Porosity: with the DG discretisation, point_data needs a CG1 projection
        # (the same DG->CG1 interpolation the cluster field uses). The CG path
        # writes its Lagrange-1 field directly, so no helper is needed there.
        if (
            on.get("porosity", False)
            and getattr(problem, "porosity", None) is not None
            and getattr(problem, "porosity_discretisation", "cg") == "dg"
        ):
            self._p_cg_fn = dolfinx.fem.Function(self.V_c_cg, name="Porosity")

        # XDMF: open file once and write mesh header.
        self._ts_file = None
        if self.format == "xdmf":
            full_path = os.path.join(self.output_dir, self.filename)
            self._ts_file = XDMFFile(self.mesh.comm, full_path, "w")
            self._ts_file.write_mesh(self.mesh)

    # -----------------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------------

    def write(self, t=0.0, step=0):
        """Refresh all interpolated Functions and write the current state."""
        self._refresh()
        if self.format == "xdmf":
            self._write_timeseries(t)
        else:
            self._write_vtu(step)

    def close(self):
        """Close the underlying time-series file, if any. Safe to call repeatedly."""
        if self._ts_file is not None:
            self._ts_file.close()
            self._ts_file = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    # -----------------------------------------------------------------------
    # Internal
    # -----------------------------------------------------------------------

    # ---------- field-update helpers ----------

    def _make_interp_or_proj(self, ufl_expr, V):
        """Build either a ``dolfinx.fem.Expression`` (for interpolation) or a
        cached ``LinearProblem`` (for L2 projection), depending on whether the
        UFL expression can be directly interpolated into ``V``.

        Quadrature-sourced expressions (e.g. plasticity in ``custom`` mode,
        where the stress depends on quadrature internal variables) fail at
        ``Expression`` construction with::

            ValueError: Mismatch of tabulation points and element points.

        For those we fall back to projection. The fallback path costs one
        linear solve per write() call, which is acceptable here because the
        ``LinearProblem`` matrix is assembled once and re-used per step.
        """
        try:
            return dolfinx.fem.Expression(ufl_expr, V.element.interpolation_points)
        except Exception as e:
            print(
                f"[OutputWriter] interpolation not available for field on "
                f"{V.ufl_element()} ({type(e).__name__}); falling back to L2 projection."
            )
            u_tr = ufl.TrialFunction(V)
            v = ufl.TestFunction(V)
            dx = ufl.dx(metadata=self._metadata) if self._metadata else ufl.dx
            a = ufl.inner(u_tr, v) * dx
            L = ufl.inner(ufl_expr, v) * dx
            return dolfinx.fem.petsc.LinearProblem(
                a, L,
                petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
                petsc_options_prefix=f"z3st_proj_{id(ufl_expr):x}_",
            )

    def _u_for_output(self, u):
        """Return displacement as a CG1 vertex field for VTU/XDMF.

        VTU and legacy XDMF expect one value per mesh vertex, but a P2 (or higher)
        mechanical space carries extra dofs (edge midpoints), so ``u.x.array`` is
        longer than ``#vertices × tdim`` and a direct reshape/write fails. For such
        spaces we interpolate ``u`` once per step onto a cached CG1 vector field;
        for P1 ``u`` (already vertex-resolution) it is returned unchanged.
        """
        n_vert = self.mesh.geometry.x.shape[0]
        tdim = self.mesh.topology.dim
        if u.x.array.size == n_vert * tdim:
            return u
        if self._u_out_fn is None:
            V = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1, (tdim,)))
            self._u_out_fn = dolfinx.fem.Function(V, name="Displacement")
        self._u_out_fn.interpolate(u)
        self._u_out_fn.x.scatter_forward()
        return self._u_out_fn

    @staticmethod
    def _refresh_into(target_fn, source, cells=None):
        """Update ``target_fn`` from a pre-prepared source (Expression or cached
        projection LinearProblem). When ``cells`` is given (the local cell
        indices of a material's subdomain), the interpolation fills only those
        cells of the shared field, leaving the rest for the other materials."""
        if isinstance(source, dolfinx.fem.Expression):
            if cells is None:
                target_fn.interpolate(source)
            else:
                target_fn.interpolate(source, cells0=cells)
        else:
            # LinearProblem path (quadrature-sourced, e.g. custom plasticity).
            # Only used by single-material cases, so a whole-mesh copy is exact.
            result = source.solve()
            target_fn.x.array[:] = result.x.array[:]

    def _refresh(self):
        """Refresh every interpolated/projected Function from its source."""
        if self._strain_expr_cells is not None:
            self._refresh_into(self._strain_fn_cells, self._strain_expr_cells)
            self._refresh_into(self._strain_fn_points, self._strain_expr_points)

        for src_c, src_p, cells in self._stress_sources:
            self._refresh_into(self._stress_fn_cells, src_c, cells)
            self._refresh_into(self._stress_fn_points, src_p, cells)

        for src_c, src_p, cells in self._psi_sources:
            self._refresh_into(self._psi_fn_cells, src_c, cells)
            self._refresh_into(self._psi_fn_points, src_p, cells)

        for src_c, cells in self._heatflux_sources:
            self._refresh_into(self._heatflux_fn, src_c, cells)

        if self._D_cell_expr is not None:
            self._refresh_into(self._D_cell_fn, self._D_cell_expr)
        if self._H_cell_expr is not None:
            self._refresh_into(self._H_cell_fn, self._H_cell_expr)

        # Cluster: DG1 → CG1 Function-to-Function interpolation for XDMF/VTU.
        if self._c_cg_fn is not None:
            self._c_cg_fn.interpolate(self.problem.c)

        # Porosity (DG only): same DG1 → CG1 projection for output.
        if self._p_cg_fn is not None:
            self._p_cg_fn.interpolate(self.problem.porosity)

    def _write_timeseries(self, t):
        """Append the current state to the open XDMF file at time ``t``."""
        p = self.problem
        on = p.on
        write = self._ts_file.write_function

        if on.get("thermal", False) and p.T is not None:
            write(p.T, t)
        if on.get("mechanical", False) and p.u is not None:
            write(self._u_for_output(p.u), t)
        if self._strain_fn_cells is not None:
            write(self._strain_fn_cells, t)
        if self._stress_fn_cells is not None:
            write(self._stress_fn_cells, t)
        if self._psi_fn_cells is not None:
            write(self._psi_fn_cells, t)
        if self._heatflux_fn is not None:
            write(self._heatflux_fn, t)
        if on.get("damage", False) and p.D is not None:
            write(p.D, t)
        if self._c_cg_fn is not None:
            write(self._c_cg_fn, t)
        if getattr(p, "burnup", None) is not None:
            write(p.burnup, t)
        # SCIANTIX coupling fields (only present when models.fission_gas is on)
        if getattr(p, "gas_swelling", None) is not None:
            write(p.gas_swelling, t)
        if getattr(p, "fg_fields", None):
            for fn in p.fg_fields.values():
                write(fn, t)
        if on.get("porosity", False) and getattr(p, "porosity", None) is not None:
            write(self._p_cg_fn if self._p_cg_fn is not None else p.porosity, t)

    def _write_vtu(self, step):
        p = self.problem
        on = p.on

        topology, cell_types, geometry = dolfinx.plot.vtk_mesh(self.mesh, self.tdim)
        n_points = geometry.shape[0]
        grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

        # Nodal global fields (T, u)
        if on.get("thermal", False) and p.T is not None:
            grid.point_data["Temperature"] = p.T.x.array
        if on.get("mechanical", False) and p.u is not None:
            u_out = self._u_for_output(p.u)
            grid.point_data["Displacement"] = u_out.x.array.reshape(n_points, self.tdim)

        # Material IDs
        ct = getattr(p, "cell_tags", None) or getattr(p, "tags", None)
        if ct is not None:
            n_cells = self.mesh.topology.index_map(self.tdim).size_local
            mids = np.zeros(n_cells, dtype=np.int32)
            for name, tag in p.label_map.items():
                if "face" in name:
                    continue
                try:
                    mids[ct.find(tag)] = tag
                except RuntimeError:
                    pass
            grid.cell_data["MaterialID"] = mids

        # Strain (global)
        if self._strain_fn_cells is not None:
            grid.cell_data["Strain (cells)"] = self._strain_fn_cells.x.array.reshape(-1, 9)
            grid.point_data["Strain (points)"] = self._strain_fn_points.x.array.reshape(n_points, 9)

        # Stress + derived (von Mises, hydrostatic) — single domain-wide field,
        # each cell/node filled from its own material's stress.
        if self._stress_fn_cells is not None:
            s_cells = self._stress_fn_cells.x.array.reshape(-1, 3, 3)
            s_points = self._stress_fn_points.x.array.reshape(n_points, 3, 3)
            grid.cell_data["Stress (cells)"] = s_cells.reshape(-1, 9)
            grid.point_data["Stress (points)"] = s_points.reshape(n_points, 9)
            grid.cell_data["VonMises (cells)"] = _von_mises(s_cells)
            grid.point_data["VonMises (points)"] = _von_mises(s_points)
            grid.cell_data["Hydrostatic (cells)"] = _hydrostatic(s_cells)
            grid.point_data["Hydrostatic (points)"] = _hydrostatic(s_points)

        # Strain-energy density — single domain-wide field
        if self._psi_fn_cells is not None:
            grid.cell_data["StrainEnergyDensity (cells)"] = self._psi_fn_cells.x.array.copy()
            grid.point_data["StrainEnergyDensity (points)"] = self._psi_fn_points.x.array.copy()

        # Heat flux — single domain-wide field (cells)
        if self._heatflux_fn is not None:
            grid.cell_data["HeatFlux (cells)"] = self._heatflux_fn.x.array.reshape(-1, 3)

        # Cluster density
        if self._c_cg_fn is not None:
            grid.point_data["ClusterDensity"] = self._c_cg_fn.x.array

        # Burnup (nodal fuel state, MWd/kgU)
        if getattr(p, "burnup", None) is not None:
            grid.point_data["Burnup"] = p.burnup.x.array

        # SCIANTIX coupling fields (only present when models.fission_gas is on):
        # total gaseous swelling ΔV/V and the Xe+Kr concentrations (at/m^3).
        if getattr(p, "gas_swelling", None) is not None:
            grid.point_data["Gaseous swelling"] = p.gas_swelling.x.array
        if getattr(p, "fg_fields", None):
            for fn in p.fg_fields.values():
                grid.point_data[fn.name] = fn.x.array

        # Damage + crack-driving-force
        if on.get("damage", False) and p.D is not None:
            grid.point_data["Damage"] = p.D.x.array.copy()
            if self._D_cell_fn is not None:
                grid.cell_data["Damage"] = self._D_cell_fn.x.array.copy()
            if self._H_cell_fn is not None:
                grid.cell_data["CrackDrivingForce"] = self._H_cell_fn.x.array.copy()

        # Porosity (DG is projected to CG1 first; CG writes its nodal field).
        if on.get("porosity", False) and getattr(p, "porosity", None) is not None:
            poro_out = self._p_cg_fn if self._p_cg_fn is not None else p.porosity
            grid.point_data["Porosity"] = poro_out.x.array.copy()

        # Plasticity: cumulative plastic strain p lives on a quadrature space.
        # Projection per step is needed; the LinearProblem is built once.
        if on.get("plasticity", False) and getattr(p, "p", None) is not None:
            if self._p_proj_problem is None:
                u_tr = ufl.TrialFunction(self.V_scalar_cells)
                v = ufl.TestFunction(self.V_scalar_cells)
                dx = ufl.dx(metadata=self._metadata) if self._metadata else ufl.dx
                a = ufl.inner(u_tr, v) * dx
                L = ufl.inner(p.p, v) * dx
                self._p_proj_problem = dolfinx.fem.petsc.LinearProblem(
                    a, L,
                    petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
                    petsc_options_prefix="z3st_proj_p_",
                )
            p_projected = self._p_proj_problem.solve()
            grid.cell_data["CumulativePlasticStrain"] = p_projected.x.array.copy()

        # File naming: single file for n_steps==1, otherwise step-numbered.
        if self.n_steps == 1:
            out_path = os.path.join(self.output_dir, self.filename)
        else:
            base, _ext = os.path.splitext(self.filename)
            out_path = os.path.join(self.output_dir, f"{base}_{step:04d}.vtu")
        grid.save(out_path)
