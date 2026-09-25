#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for verification/cohesive/bar_1D.

Streams the structural response to ``output/response.csv``: one row per step
with the prescribed displacement, the axial stress and the maximum of the phase
field. ``non-regression.py`` reads this CSV and compares the whole curve with
the closed-form solution of the paper.

The bar is in equilibrium under a uniform axial stress, so the stress is taken
as its volume average, which is exact here and insensitive to which element
localizes.
"""

import os

import dolfinx
import numpy as np
import ufl
from mpi4py import MPI

_CSV = os.path.join(os.path.dirname(__file__), "output", "response.csv")
_HEADER = "step,U_t_m,sigma_xx_Pa,alpha_max,E_el_J,E_frac_J\n"

_run_started = False


def per_step(problem, step, t):
    global _run_started

    material = problem.materials["bar"]
    u, p, q = problem.split_state(problem.w)

    # Volume-averaged axial stress (uniform along the bar at equilibrium).
    sigma = problem.sigma_cohesive(u, p, q, material)
    dx = problem.dx_tags[problem.label_map["bar"]]
    comm = problem.mesh.comm
    sigma_int = comm.allreduce(
        dolfinx.fem.assemble_scalar(dolfinx.fem.form(sigma[0, 0] * dx)), op=MPI.SUM
    )
    length = comm.allreduce(
        dolfinx.fem.assemble_scalar(dolfinx.fem.form(dolfinx.fem.Constant(
            problem.mesh, dolfinx.default_scalar_type(1.0)) * dx)), op=MPI.SUM
    )
    sigma_avg = sigma_int / length

    # Prescribed displacement of the loaded end, as the solver applied it.
    U_t = 0.0
    for bc_list in problem.dirichlet_mechanical.values():
        for bc in bc_list:
            if isinstance(bc, dict) and isinstance(bc.get("raw"), list):
                U_t = float(np.atleast_1d(bc["const"].value)[0])

    alpha_max = comm.allreduce(float(problem.D.x.array.max()), op=MPI.MAX)
    E_el, E_frac = problem.compute_cohesive_energy_balance(problem.w, problem.D)

    if comm.rank != 0:
        return

    os.makedirs(os.path.dirname(_CSV), exist_ok=True)
    mode = "w" if not _run_started else "a"
    with open(_CSV, mode) as f:
        if not _run_started:
            f.write(_HEADER)
        f.write(f"{step:d},{U_t:.10e},{sigma_avg:.10e},{alpha_max:.10e},"
                f"{E_el:.10e},{E_frac:.10e}\n")
    _run_started = True
