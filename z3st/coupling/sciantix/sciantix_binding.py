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

STATUS: draft. The index map and array sizes below are read from
SCIANTIX ``include/MainVariables.h`` and ``src/operations/SetVariablesFunctions.C``
(v2.2.1); the indices marked (confirm) should be checked by the SCIANTIX author
against the current source before production use. Not yet run end-to-end here —
needs the shared library (README.md).
"""

import ctypes
import os

import numpy as np

# --. array sizes (SCIANTIX include/MainVariables.h, v2.2.1) --..
OPTIONS_SIZE = 40
HISTORY_SIZE = 20
VARIABLES_SIZE = 300
DIFFUSION_MODES_SIZE = 8000   # (confirm) SCIANTIX spectral-mode state; generous, zero-init

# --. history[] layout — what the HOST writes each step (old/new pairs) --..
# (src/operations/SetVariablesFunctions.C::initializeHistoryVariable)
H_TEMPERATURE_OLD = 0     # (K)            x scaling_factors[4]
H_TEMPERATURE_NEW = 1
H_FISSION_RATE_OLD = 2    # (fiss/m^3 s)   x scaling_factors[5]
H_FISSION_RATE_NEW = 3
H_HYDRO_STRESS_OLD = 4    # (MPa)
H_HYDRO_STRESS_NEW = 5
H_TIME_STEP = 6          # (confirm) time-step length used by the integrator
H_STEAM_PRESSURE_OLD = 9  # (atm)
H_STEAM_PRESSURE_NEW = 10

# --. variables[] layout — outputs the host READS --..
# (src/operations/SetVariablesFunctions.C::initializeSciantixVariable)
V_XE_PRODUCED = 1               # (at/m^3)
V_XE_RELEASED = 6               # (at/m^3)
V_INTRAGRAN_GAS_SWELLING = 24   # (/)  intragranular gas-bubble swelling, DV/V
V_INTERGRAN_GAS_SWELLING = 36   # (/)  intergranular gas swelling, DV/V  (dominant at high T)
V_BURNUP = 38                   # (MWd/kgUO2)


class SciantixSolver:
    """One SCIANTIX integration point (e.g. one fuel radial ring).

    Holds the persistent state arrays and advances them one step per call.
    """

    def __init__(self, libpath=None):
        libpath = libpath or os.environ.get("SCIANTIX_LIB", "libsciantix.so")
        self._lib = ctypes.CDLL(libpath)
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_dbl_p = ctypes.POINTER(ctypes.c_double)
        self._lib.callSciantix.argtypes = [c_int_p, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p]
        self._lib.callSciantix.restype = None
        self._lib.getSciantixOptions.argtypes = [c_int_p, c_dbl_p]
        self._lib.getSciantixOptions.restype = None

        # persistent state (contiguous, C-typed, owned by this instance)
        self.options = np.zeros(OPTIONS_SIZE, dtype=np.int32)
        self.history = np.zeros(HISTORY_SIZE, dtype=np.float64)
        self.variables = np.zeros(VARIABLES_SIZE, dtype=np.float64)
        self.scaling_factors = np.ones(OPTIONS_SIZE, dtype=np.float64)
        self.diffusion_modes = np.zeros(DIFFUSION_MODES_SIZE, dtype=np.float64)

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

    def advance(self, dt, T_old, T_new, fission_rate, hydro_stress=0.0, steam_pressure=0.0):
        """Advance this point by dt (s) with the given local conditions and
        return a dict of the engineering outputs.

        Parameters
        ----------
        dt            : time-step length (s)
        T_old, T_new  : temperature at the start/end of the step (K)
        fission_rate  : volumetric fission rate (fiss / m^3 s); held over the step
        hydro_stress  : hydrostatic stress (MPa), optional
        steam_pressure: steam pressure (atm), optional
        """
        h = self.history
        h[H_TEMPERATURE_OLD], h[H_TEMPERATURE_NEW] = T_old, T_new
        h[H_FISSION_RATE_OLD] = h[H_FISSION_RATE_NEW] = fission_rate
        h[H_HYDRO_STRESS_OLD] = h[H_HYDRO_STRESS_NEW] = hydro_stress
        h[H_STEAM_PRESSURE_OLD] = h[H_STEAM_PRESSURE_NEW] = steam_pressure
        h[H_TIME_STEP] = dt

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
