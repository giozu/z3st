# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST <-> SCIANTIX coupling — Python binding (PROTOTYPE / DRAFT)
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Thin ctypes binding to the SCIANTIX C-linkage coupling entry point.

SCIANTIX already exposes, for the TRANSURANUS coupling, an ``extern "C"`` API in
``coupling/TUSrcCoupling.h``::

    void getSciantixOptions(int options[], double scaling_factors[]);
    void callSciantix(int options[], double history[], double variables[],
                      double scaling_factors[], double diffusion_modes[]);

so no new C/C++ is needed beyond compiling SCIANTIX as a SHARED library (see
README.md). This wrapper loads that library, owns the per-point state arrays,
and exposes a single ``advance(dt, T0, T1, fission_rate, ...)`` step that returns
the engineering outputs Z3ST consumes — gaseous swelling (for the eigenstrain
bus) and burnup / FGR.

The design intentionally mirrors the way SCIANTIX is *meant* to be driven by a
host code: arrays in, arrays out, all state carried in ``history`` (old/new
pairs) and ``variables`` between calls — exactly the per-point, stateful-via-
arrays pattern Z3ST already uses for the burnup and creep DG0 states. One
``SciantixSolver`` instance == one integration point (a radial ring, to start).

STATUS: draft. The index map and array sizes below were verified 2026-06-18
against SCIANTIX v2.2.1 source (``include/MainVariables.h``,
``src/operations/SetVariablesFunctions.C``, ``src/operations/SetVariables.C``,
``src/coupling/TUSrcCoupling.C``): array sizes options[40]/history[20]/
variables[300]/scaling_factors[20]/diffusion_modes[720]; history slots and the
swelling(24/36)/burnup(38) output slots all confirmed; ``history[6]`` is the
time step in seconds (drives ``physics_variable["Time step"]``).

VALIDATED end-to-end 2026-06-18 against SCIANTIX's own standalone regression case
``regression/baker/test_Baker1977__1273K`` (see ``validate_baker.py``): built
SCIANTIX as a shared lib (g++ -fPIC -shared over src/*.C), seeded initial
conditions via ``load_initial_conditions``, and drove the same 1273 K / 1e19
fiss/m3 s / 100x55 h history. All four engineering outputs match the standalone
gold to ~1e-7 relative error (FGR 0.132097, intragranular swelling 3.07e-4,
intergranular swelling 0.0417, burnup 6.719).

Two non-obvious facts the validation surfaced, both handled in
``load_initial_conditions`` / ``_apply_initialization``:
  1. In coupling mode SCIANTIX does NOT read input_initial_conditions.txt (that is
     standalone-only, file_manager/InputReading.C) — the host must seed
     ``variables[]`` first (grain radius, fuel density, U content, ...).
  2. The standalone's one-time ``Initialization()`` (file_manager/Initialization.C),
     also skipped by the coupling entry, sets grain-boundary defaults absent from
     the IC file (``variables[25]``=2e13, ``[35]``=0.5, ``[37]``=1.0) and converts
     U% -> at/m3 using density. Without the grain-boundary defaults the
     intergranular model returns nan and releases nothing.

Burnup ownership (2026-06-19): the production design is a -DCOUPLING_TU build, where
SCIANTIX skips its own Burnup()/EffectiveBurnup()/Densification (Simulation.C:43)
and reads burnup from history[7] (old) / history[8] (new) (SetVariables.C:74). The
HOST owns burnup — Z3ST computes it with its RADAR model and feeds it via
``advance(..., burnup_old=, burnup_new=)``. Verified against the Baker gold with a
-DCOUPLING_TU lib: feeding the gold burnup trajectory reproduces gas swelling/FGR to
~1e-7 and burnup exactly.

Z3ST integration is wired in (SciantixField + spine + materials.sciantix_swelling;
README section 4). Remaining: fresh-fuel only (the diffusion-mode projection for a
pre-irradiated restart is not implemented).
"""

import ctypes
import os

import numpy as np

# --. array sizes (SCIANTIX include/MainVariables.h, v2.2.1; verified 2026-06-18) --..
OPTIONS_SIZE = 40
HISTORY_SIZE = 20
VARIABLES_SIZE = 300
SCALING_FACTORS_SIZE = 20
DIFFUSION_MODES_SIZE = 720    # = N_MODE_BLOCKS(18) * N_DIFFUSION_MODES(40), exact

# --. history[] layout — what the HOST writes each step (old/new pairs) --..
# (src/operations/SetVariablesFunctions.C::initializeHistoryVariable)
H_TEMPERATURE_OLD = 0     # (K)            x scaling_factors[4]
H_TEMPERATURE_NEW = 1
H_FISSION_RATE_OLD = 2    # (fiss/m^3 s)   x scaling_factors[5]
H_FISSION_RATE_NEW = 3
H_HYDRO_STRESS_OLD = 4    # (MPa)
H_HYDRO_STRESS_NEW = 5
H_TIME_STEP = 6           # (s) -> physics_variable["Time step"] (SetVariables.C:47)
# Burnup transfer (COUPLING build only): with -DCOUPLING_TU, SCIANTIX SKIPS its own
# Burnup/EffectiveBurnup/Densification (Simulation.C:43) and READS burnup from these
# two history slots instead — old from history[7], new from history[8]
# (SetVariables.C:74). So the HOST owns burnup (Z3ST RADAR model) and feeds it in.
# In a NON-coupling build these slots are Time(h)/step-number (output-only) instead.
H_BURNUP_OLD = 7          # (MWd/kgUO2)  COUPLING_TU build only
H_BURNUP_NEW = 8          # (MWd/kgUO2)  COUPLING_TU build only
H_STEAM_PRESSURE_OLD = 9  # (atm)
H_STEAM_PRESSURE_NEW = 10

# --. variables[] layout — outputs the host READS --..
# (src/operations/SetVariablesFunctions.C::initializeSciantixVariable; verified
# against the named SciantixVariable("Xe ...", variables[i]) registrations)
V_XE_PRODUCED = 1               # (at/m^3)  "Xe produced"
V_XE_INGRAIN = 2                # (at/m^3)  "Xe in grain"          (intragranular)
V_XE_GRAINBOUNDARY = 5          # (at/m^3)  "Xe at grain boundary"
V_XE_RELEASED = 6               # (at/m^3)  "Xe released"
V_KR_PRODUCED = 7               # (at/m^3)  "Kr produced"
V_KR_INGRAIN = 8                # (at/m^3)  "Kr in grain"
V_KR_GRAINBOUNDARY = 11         # (at/m^3)  "Kr at grain boundary"
V_KR_RELEASED = 12              # (at/m^3)  "Kr released"
V_INTRAGRAN_GAS_SWELLING = 24   # (/)  intragranular gas-bubble swelling, DV/V
V_INTERGRAN_GAS_SWELLING = 36   # (/)  intergranular gas swelling, DV/V  (dominant at high T)
V_BURNUP = 38                   # (MWd/kgUO2)

# Total fission gas (Xe + Kr) in each state -> (produced, in-grain, gb, released)
# variables[] index pairs; summed to a single per-state concentration (at/m^3).
_FG_STATES = {
    "produced":       (V_XE_PRODUCED, V_KR_PRODUCED),
    "in_grain":       (V_XE_INGRAIN, V_KR_INGRAIN),
    "grain_boundary": (V_XE_GRAINBOUNDARY, V_KR_GRAINBOUNDARY),
    "released":       (V_XE_RELEASED, V_KR_RELEASED),
}

# --. initial-conditions slot map: (variables[] start index, count) per data
# line of a standalone input_initial_conditions.txt, IN FILE ORDER. Mirrors
# src/file_manager/InputReading.C:154-232. The coupling entry (callSciantix)
# does NOT read that file, so the host must seed these before the first call.
# Values are passed verbatim (the U line is % of heavy atoms; SCIANTIX converts
# to at/m^3 internally using the fuel density).
_IC_SLOTS = [
    (0, 1),    # grain radius (m)
    (1, 6),    # Xe:  produced, intragranular, intra-solution, intra-bubbles, gb, released
    (7, 6),    # Kr:  same six
    (13, 6),   # He:  same six
    (19, 2),   # intragranular bubbles: concentration (at/m^3), radius (m)
    (38, 1),   # burnup (MWd/kgUO2)
    (39, 1),   # effective burnup (MWd/kgUO2)
    (65, 1),   # irradiation time (h)
    (40, 1),   # fuel density (kg/m^3)
    (41, 5),   # U234 U235 U236 U237 U238 (% of heavy atoms)
    (48, 7),   # Xe133: produced, intragranular, solution, bubbles, decayed, gb, released
    (57, 7),   # Kr85m: same seven
    (66, 1),   # stoichiometry deviation
    (150, 1),  # chromium content (ppm)
]


class SciantixSolver:
    """One SCIANTIX integration point (e.g. one fuel radial ring).

    Holds the persistent state arrays and advances them one step per call.
    """

    def __init__(self, libpath=None, lib=None, options=None, scaling_factors=None):
        """Create one integration point.

        ``lib`` (a pre-loaded ``ctypes.CDLL``) and ``options`` / ``scaling_factors``
        (numpy arrays from another point) can be passed to SHARE them across a whole
        field of points — this avoids re-loading the ``.so`` and re-reading
        ``input_settings.txt`` per point (see :class:`SciantixField`). When
        ``options`` is omitted, the per-point defaults are pulled from SCIANTIX via
        ``getSciantixOptions`` (which reads ``input_settings.txt`` from the CWD).
        """
        if lib is not None:
            self._lib = lib
        else:
            libpath = libpath or os.environ.get("SCIANTIX_LIB", "libsciantix.so")
            self._lib = ctypes.CDLL(libpath)
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_dbl_p = ctypes.POINTER(ctypes.c_double)
        self._lib.callSciantix.argtypes = [c_int_p, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p]
        self._lib.callSciantix.restype = None
        self._lib.getSciantixOptions.argtypes = [c_int_p, c_dbl_p]
        self._lib.getSciantixOptions.restype = None

        # persistent state (contiguous, C-typed, owned by this instance)
        self.history = np.zeros(HISTORY_SIZE, dtype=np.float64)
        self.variables = np.zeros(VARIABLES_SIZE, dtype=np.float64)
        self.diffusion_modes = np.zeros(DIFFUSION_MODES_SIZE, dtype=np.float64)

        if options is not None:
            # share the model selection / scaling factors read once elsewhere
            self.options = np.array(options, dtype=np.int32)
            self.scaling_factors = np.array(scaling_factors, dtype=np.float64)
        else:
            self.options = np.zeros(OPTIONS_SIZE, dtype=np.int32)
            self.scaling_factors = np.ones(SCALING_FACTORS_SIZE, dtype=np.float64)
            # pull SCIANTIX's default options + scaling factors (which models are on)
            self._lib.getSciantixOptions(self._ip(self.options), self._dp(self.scaling_factors))
        self._initialised = False

    # ctypes pointer helpers
    @staticmethod
    def _ip(a):
        return a.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    @staticmethod
    def _dp(a):
        return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    def load_initial_conditions(self, path="input_initial_conditions.txt"):
        """Seed the persistent ``variables[]`` state from a SCIANTIX standalone
        ``input_initial_conditions.txt``.

        REQUIRED before the first ``advance``: the coupling entry
        (``callSciantix``) does not read this file — that is standalone-only
        (``file_manager/InputReading.C``) — so without seeding, fuel density and
        heavy-metal content are zero and burnup -> inf, gas -> nan. This mirrors
        ``InputReading.C:154-232``: data lines (blank/``#`` lines skipped) are
        consumed in file order and mapped to ``variables[]`` via ``_IC_SLOTS``.
        Returns self so it can be chained after the constructor.
        """
        data_lines = []
        with open(path) as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                data_lines.append([float(tok) for tok in s.replace(",", " ").split()])
        if len(data_lines) < len(_IC_SLOTS):
            raise ValueError(
                f"{path}: {len(data_lines)} data lines, expected >= {len(_IC_SLOTS)}"
            )
        for (start, count), vals in zip(_IC_SLOTS, data_lines):
            if len(vals) < count:
                raise ValueError(
                    f"{path}: data line for variables[{start}:{start+count}] has "
                    f"{len(vals)} values, expected {count}"
                )
            self.variables[start:start + count] = vals[:count]
        self._apply_initialization()
        return self

    def _apply_initialization(self):
        """Apply the derived initial state that the standalone does in its one-time
        ``Initialization()`` (``file_manager/Initialization.C``) — the coupling
        entry skips it entirely. Without the grain-boundary defaults the
        intergranular model divides by zero (-> nan swelling, zero release).
        Run once, after the raw ``variables[]`` are seeded from file.
        """
        v = self.variables
        # grain-boundary defaults (NOT in input_initial_conditions.txt; hardcoded)
        v[25] = 2.0e13   # intergranular bubble concentration (bub/m2)
        v[35] = 0.5      # intergranular saturation fractional coverage
        v[37] = 1.0      # intergranular fractional intactness
        # heavy-metal content: % of heavy atoms -> at/m^3 (uses fuel density v[40];
        # 6.022e24 = N_A per kg, 0.8815 = U mass fraction in UO2)
        density = v[40]
        for idx, mass in ((41, 234.04095), (42, 235.04393), (43, 236.04557),
                          (44, 237.04873), (45, 238.05079)):
            v[idx] *= density * 6.022e24 * 0.8815 / mass
        v[64] = 1.0      # intragranular similarity ratio
        # fabrication / residual porosity (theoretical density 10960 kg/m^3)
        v[70] = v[71] = 1.0 - density / 10960.0
        v[73] = 0.75 * v[71]
        # NB: Initialization() also projects the initial intra-granular gas onto the
        # diffusion modes. For fresh fuel (all initial gas = 0) that projection is
        # identically zero, so the zero-initialised diffusion_modes are already
        # correct. A pre-irradiated restart (non-zero v[2..4], v[8..10], v[14..16])
        # would additionally need that spectral projection — not yet implemented.
        if any(self.variables[i] != 0.0 for i in (2, 3, 4, 8, 9, 10, 14, 15, 16)):
            print("[WARNING] non-zero initial intra-granular gas: diffusion-mode "
                  "projection not implemented; restart state will be approximate.")

    def advance(self, dt, T_old, T_new, fission_rate, hydro_stress=0.0,
                steam_pressure=0.0, burnup_old=None, burnup_new=None):
        """Advance this point by dt (s) with the given local conditions and
        return a dict of the engineering outputs.

        Parameters
        ----------
        dt            : time-step length (s)
        T_old, T_new  : temperature at the start/end of the step (K)
        fission_rate  : volumetric fission rate (fiss / m^3 s); held over the step
        hydro_stress  : hydrostatic stress (MPa), optional
        steam_pressure: steam pressure (atm), optional
        burnup_old,
        burnup_new    : host-computed burnup (MWd/kgUO2) at the start/end of the
                        step. REQUIRED for a -DCOUPLING_TU build, where SCIANTIX does
                        not compute burnup and reads it from history[7]/[8] (the
                        intended design: Z3ST's RADAR model owns burnup). Leave None
                        for a non-coupling build (SCIANTIX integrates burnup itself).
        """
        h = self.history
        h[H_TEMPERATURE_OLD], h[H_TEMPERATURE_NEW] = T_old, T_new
        h[H_FISSION_RATE_OLD] = h[H_FISSION_RATE_NEW] = fission_rate
        h[H_HYDRO_STRESS_OLD] = h[H_HYDRO_STRESS_NEW] = hydro_stress
        h[H_STEAM_PRESSURE_OLD] = h[H_STEAM_PRESSURE_NEW] = steam_pressure
        h[H_TIME_STEP] = dt
        if burnup_new is not None:
            h[H_BURNUP_OLD] = burnup_new if burnup_old is None else burnup_old
            h[H_BURNUP_NEW] = burnup_new

        self._lib.callSciantix(
            self._ip(self.options), self._dp(self.history), self._dp(self.variables),
            self._dp(self.scaling_factors), self._dp(self.diffusion_modes),
        )
        self._initialised = True
        return self.results()

    def results(self):
        v = self.variables
        xe_prod = v[V_XE_PRODUCED]
        intrag = v[V_INTRAGRAN_GAS_SWELLING]
        interg = v[V_INTERGRAN_GAS_SWELLING]
        return {
            "gaseous_swelling": float(intrag + interg),   # total DV/V -> eigenstrain bus
            "intragranular_gas_swelling": float(intrag),
            "intergranular_gas_swelling": float(interg),
            "burnup_MWd_kgUO2": float(v[V_BURNUP]),
            "fission_gas_release": float(v[V_XE_RELEASED] / xe_prod) if xe_prod > 0 else 0.0,
        }


# --.. ..- .-.. .-.. --- field driver (one SCIANTIX point per dof) --.. ..- .-.. .-.. ---
class SciantixField:
    """A field of SCIANTIX integration points, one per finite-element dof.

    This is what Z3ST drives: ``spine.initialize_fields`` builds one of these over
    the fissile dofs, ``spine.update_state(dt)`` calls :meth:`step` with the local
    temperature and fission-rate arrays each step, and the returned gaseous-swelling
    array is written into a dolfinx Function that the ``sciantix_swelling``
    eigenstrain callable reads (the eigenstrain bus). The library is loaded once and
    the model selection (``options``/``scaling_factors``) is read once, then shared
    across all points; only the per-point state arrays (history/variables/
    diffusion_modes) are independent.

    The SCIANTIX C++ ``Simulation`` is a singleton scratch object re-initialised
    from the passed arrays on every call, so points are independent as long as each
    owns its own arrays — which they do.
    """

    def __init__(self, n_points, libpath=None, ic_path="input_initial_conditions.txt"):
        libpath = libpath or os.environ.get("SCIANTIX_LIB", "libsciantix.so")
        self._lib = ctypes.CDLL(libpath)
        # one template reads the models + initial conditions; the rest clone them
        template = SciantixSolver(lib=self._lib)
        if ic_path is not None:
            template.load_initial_conditions(ic_path)
        self.points = [template]
        for _ in range(1, n_points):
            pt = SciantixSolver(lib=self._lib, options=template.options,
                                scaling_factors=template.scaling_factors)
            pt.variables[:] = template.variables          # same seeded initial state
            pt.diffusion_modes[:] = template.diffusion_modes
            self.points.append(pt)
        self.n_points = n_points

    def step(self, dt, T, fission_rate, hydro_stress=None,
             burnup_old=None, burnup_new=None):
        """Advance every point by ``dt`` (s) and return the per-point total gaseous
        swelling ΔV/V as a numpy array (length ``n_points``).

        ``T`` and ``fission_rate`` are per-point arrays (K, fiss/m^3 s);
        ``hydro_stress`` is an optional per-point array (MPa). ``burnup_old`` /
        ``burnup_new`` are per-point host burnup arrays (MWd/kgUO2) — pass them for a
        -DCOUPLING_TU build so SCIANTIX consumes Z3ST's RADAR burnup instead of
        computing its own (see :meth:`SciantixSolver.advance`).
        """
        out = np.empty(self.n_points, dtype=np.float64)
        for i, pt in enumerate(self.points):
            hs = 0.0 if hydro_stress is None else float(hydro_stress[i])
            bo = None if burnup_old is None else float(burnup_old[i])
            bn = None if burnup_new is None else float(burnup_new[i])
            res = pt.advance(dt, float(T[i]), float(T[i]), float(fission_rate[i]),
                             hydro_stress=hs, burnup_old=bo, burnup_new=bn)
            out[i] = res["gaseous_swelling"]
        return out

    def gas_concentrations(self):
        """Per-dof total fission-gas (Xe + Kr) concentrations in each state.

        Reads the *current* ``variables[]`` of every point (call after :meth:`step`)
        and returns a dict of numpy arrays, each length ``n_points`` (units at/m^3):

            ``produced``        gas created so far (Xe + Kr produced)
            ``in_grain``        retained intragranular (in grain)
            ``grain_boundary``  accumulated on grain faces (intergranular)
            ``released``        released to the free volume (FGR numerator)

        These are diagnostic/output fields (the eigenstrain only needs the swelling
        from :meth:`step`); Z3ST writes them to the Paraview output and the FG plots.
        """
        out = {key: np.empty(self.n_points, dtype=np.float64) for key in _FG_STATES}
        for i, pt in enumerate(self.points):
            v = pt.variables
            for key, (i_xe, i_kr) in _FG_STATES.items():
                out[key][i] = v[i_xe] + v[i_kr]
        return out

    # --. checkpoint / restore — for adaptive time-stepping rollback (CODE-FEATURE-4
    #     shares this invariant with the adaptive controller; see punch_list) --..
    def snapshot(self):
        return [(pt.variables.copy(), pt.diffusion_modes.copy(), pt.history.copy())
                for pt in self.points]

    def restore(self, snap):
        for pt, (var, modes, hist) in zip(self.points, snap):
            pt.variables[:] = var
            pt.diffusion_modes[:] = modes
            pt.history[:] = hist
