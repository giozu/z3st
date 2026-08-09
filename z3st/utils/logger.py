# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.2.0 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---

"""Framework-wide logger.

It writes to **stdout**, not stderr, so its output joins the same stream the
`print` diagnostics use: `__main__` wraps stdout in a markdown filter, the case
Allrun redirects stdout to `log_z3st.md`, and CI dumps the tail of that file when
a case fails. On stderr — where a plain `logging.basicConfig` puts it — the mesh
diagnostics never reached the log file at all, so a case that died on a mesh
problem left no record of which mesh it had loaded.

The formatter emits `[LEVEL] message`, which the markdown filter already
recognises and renders as `**[LEVEL]** message`.
"""

import logging
import sys

from mpi4py import MPI


class _LateBoundStdout:
    """Resolves the target stream at write time, not at import time.

    ``__main__`` wraps ``sys.stdout`` early; binding a handler to whatever object
    happens to be installed when this module is first imported would make behaviour
    depend on import order. Forwarding on every write removes that coupling.

    That wrapper also carries the MPI rank gate, and it lets ``[WARNING]`` and
    ``[ERROR]`` lines through from every rank -- which is why the formatter below
    emits the level in brackets.
    """

    def write(self, message):
        return sys.stdout.write(message)

    def flush(self):
        sys.stdout.flush()


class _RankFilter(logging.Filter):
    """Chatter from rank 0 only; warnings and errors from every rank.

    Under MPI every rank runs the same code, so an unguarded INFO line appears
    once per process. Warnings and errors are let through everywhere because a
    failure on one rank is exactly what you need to see.

    Kept even though ``__main__``'s stdout wrapper gates the same way: this is the
    only gate when z3st is imported as a library, without going through
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

# Own handler only: propagating would also hit the root logger and duplicate every
# line. This also keeps z3st from configuring logging for third-party libraries,
# which an earlier logging.basicConfig here did as a side effect.
log.propagate = False
