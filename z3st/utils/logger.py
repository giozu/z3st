# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""Framework-wide logger.

It writes to **stdout**, not stderr, so its output joins the same stream the
`print` diagnostics use: `__main__` wraps stdout in a markdown filter, the case
Allrun redirects stdout to `log_z3st.md`, and CI dumps the tail of that file when
a case fails.

The formatter emits `[LEVEL] message`, which the markdown filter already
recognises and renders as `**[LEVEL]** message`.
"""

import logging
import sys

from mpi4py import MPI


class _LateBoundStdout:
    """Resolves the target stream at write time, not at import time.

    ``__main__`` wraps ``sys.stdout`` early, so forwarding on every write keeps the
    handler independent of import order. That wrapper carries the MPI rank gate and
    lets ``[WARNING]`` and ``[ERROR]`` lines through from every rank.
    """

    def write(self, message):
        return sys.stdout.write(message)

    def flush(self):
        sys.stdout.flush()


class _RankFilter(logging.Filter):
    """Chatter from rank 0 only; warnings and errors from every rank.

    Under MPI every rank runs the same code, so an unguarded INFO line appears
    once per process.

    The only gate when z3st is imported as a library, without going through
    ``python -m z3st``.
    """

    def __init__(self):
        super().__init__()
        self.rank = MPI.COMM_WORLD.rank

    def filter(self, record):
        return self.rank == 0 or record.levelno >= logging.WARNING


log = logging.getLogger("z3st")
log.setLevel(logging.INFO)

if not log.handlers:
    _handler = logging.StreamHandler(_LateBoundStdout())
    _handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    _handler.addFilter(_RankFilter())
    log.addHandler(_handler)

# Own handler only, no propagation to the root logger: no duplicated lines, and
# z3st does not configure logging for third-party libraries.
log.propagate = False
