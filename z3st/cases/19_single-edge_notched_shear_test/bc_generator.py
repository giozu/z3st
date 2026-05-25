"""Generate the boundary-conditions YAML.

A top-edge shear-displacement ramp paired with a fixed bottom edge and a
unit-damage Dirichlet on the crack. The ramp
is defined as a sequence of segments ``(u_start, u_end, delta_u)``, so
the resolution can be finer around the bifurcation point and coarser
in the elastic ramp-up and post-peak softening tail.

For a uniform single-stage ramp, pass a one-element ``SEGMENTS`` list.

Re-run this script after editing the segments; the simulation's
``input.yaml::n_steps`` must match the step count printed below.
"""

import yaml

class IndentedList(list):
    """Marker subclass: signals that this sequence should render in YAML
    flow style (e.g. ``[ux, uy]``) rather than expanded block style."""
    pass

def _indented_list_presenter(dumper, data):
    return dumper.represent_sequence(
        "tag:yaml.org,2002:seq", data, flow_style=True
    )

yaml.add_representer(IndentedList, _indented_list_presenter)

# Tolerance for floating-point comparisons at segment boundaries (m).
# Tighter than the round-to-12-decimals applied below, looser than 0.
_FP_TOL = 1.0e-12

def make_ramp_segments(segments: list) -> list:
    """Build the per-step displacement list from a multi-stage ramp.

    Each segment ``(u_start, u_end, delta_u)`` generates increments of
    size ``delta_u`` from ``u_start`` up to and including ``u_end``. The
    first segment contributes its ``u_start`` as step 0; every subsequent
    segment skips its first point to avoid duplicating the previous
    segment's endpoint. Segments must be contiguous (``u_end`` of one
    equals ``u_start`` of the next, within ``_FP_TOL``) and each segment's
    span must be an integer multiple of its ``delta_u``.

    Parameters
    ----------
    segments : list of (float, float, float)
        ``(u_start_m, u_end_m, delta_u_m)`` for each segment. Pass a
        single-element list for a uniform ramp.

    Returns
    -------
    list of IndentedList
        Per-step ``[u_x, u_y]`` vectors with ``u_y == 0.0``, ready to be
        embedded under ``mechanical.steel[ymax].displacement`` in the
        boundary-conditions YAML.

    Raises
    ------
    ValueError
        On any of: empty segment list, non-positive number of increments
        in a segment, segment span not an integer multiple of its
        ``delta_u``, or non-contiguous segments.

    Examples
    --------
    Three-stage Ambati-style ramp (fine around the bifurcation):

    >>> segs = [
    ...     (0.0,      8.0e-6,  1.0e-7),
    ...     (8.0e-6,  15.0e-6,  1.0e-8),
    ...     (15.0e-6, 30.0e-6,  1.0e-7),
    ... ]
    >>> len(make_ramp_segments(segs))
    931
    """
    if not segments:
        raise ValueError("segments must be a non-empty list of (u_start, u_end, delta_u).")

    # Validate continuity between consecutive segments.
    for i in range(1, len(segments)):
        prev_end = segments[i - 1][1]
        cur_start = segments[i][0]
        if abs(prev_end - cur_start) > _FP_TOL:
            raise ValueError(
                f"segments[{i}].u_start = {cur_start} != "
                f"segments[{i-1}].u_end = {prev_end}; segments must be contiguous."
            )

    # Validate each segment internally and emit points.
    u_vals = []
    for i, (u0, u1, du) in enumerate(segments):
        if u1 <= u0:
            raise ValueError(
                f"segments[{i}]: u_end ({u1}) must be > u_start ({u0})."
            )
        if du <= 0.0:
            raise ValueError(f"segments[{i}]: delta_u must be positive (got {du}).")

        n_inc_float = (u1 - u0) / du
        n_inc = int(round(n_inc_float))
        if n_inc <= 0:
            raise ValueError(
                f"segments[{i}]: delta_u ({du}) >= segment span ({u1 - u0}); "
                f"need at least one positive increment."
            )
        if abs(n_inc * du - (u1 - u0)) > _FP_TOL:
            raise ValueError(
                f"segments[{i}]: (u_end - u_start)/delta_u = {n_inc_float} "
                f"is not an integer (off by {abs(n_inc - n_inc_float):.3g}). "
                f"Each segment span must be an integer multiple of its delta_u."
            )

        # Skip k=0 on every segment except the first, to avoid
        # duplicating the previous segment's endpoint.
        k_start = 0 if i == 0 else 1
        for k in range(k_start, n_inc + 1):
            u_vals.append(round(u0 + k * du, 12))

    return [IndentedList([u, 0.0]) for u in u_vals]


def make_ramp(u0: float, u1: float, delta_u: float) -> list:
    """Single-segment convenience wrapper around `make_ramp_segments`."""
    return make_ramp_segments([(u0, u1, delta_u)])


def write_bc_yaml(steps: list, filename: str = "boundary_conditions.yaml") -> None:
    """Write the SENS boundary-conditions YAML to ``filename``."""
    zero_disp = IndentedList([0.0, 0.0])

    bc_data = {
        "mechanical": {
            "steel": [
                {"type": "Dirichlet", "region": "ymin", "displacement": zero_disp},
                {"type": "Dirichlet", "region": "ymax", "displacement": steps},
            ],
        },
        "damage": {
            "steel": [
                {"type": "Dirichlet", "region": "crack", "value": 1.0},
            ],
        },
    }

    with open(filename, "w") as f:
        yaml.dump(bc_data, f, sort_keys=False, default_flow_style=False, indent=2)

    print(f"'{filename}' generated.")


def _print_segment_summary(segments: list) -> None:
    """Pretty-print the ramp segments to stdout."""
    print()
    print("Ramp summary (multi-stage):")
    print(f"  {'segment':>7s}   {'u_start':>10s}   {'u_end':>10s}   "
          f"{'delta_u':>12s}   {'increments':>10s}")
    total = 1  # the first point (u_start of segment 0)
    for i, (u0, u1, du) in enumerate(segments):
        n_inc = int(round((u1 - u0) / du))
        total += n_inc
        print(f"  {i:>7d}   {u0*1e6:8.3f} um   {u1*1e6:8.3f} um   "
              f"{du*1e9:8.2f} nm   {n_inc:>10d}")
    print(f"  {'':>7s}   {'':>10s}   {'':>10s}   "
          f"{'TOTAL':>12s}   {total:>10d}  (= n_steps)")
    print()


if __name__ == "__main__":
    # ─── Ramp configuration ──────────────────────────────────────────
    # Each segment: (u_start_m, u_end_m, delta_u_m).
    #
    # Default: three-stage Ambati-style ramp.
    #   coarse pre-bifurcation     (0      ->  8 um)   delta_u = 100 nm
    #   FINE bifurcation window    (8      -> 15 um)   delta_u =  10 nm   <-- converged kink-angle window
    #   coarse softening tail      (15     -> 30 um)   delta_u = 100 nm
    #
    # For a uniform single-stage ramp, use one segment instead, e.g.
    #     SEGMENTS = [(0.0, 30.0e-6, 1.0e-8)]    # 3001 steps, ~3-5 h wall-time
    #     SEGMENTS = [(0.0, 30.0e-6, 1.0e-7)]    # 301 steps, loose (Ambati's "larger" delta_u)
    # ─────────────────────────────────────────────────────────────────
    SEGMENTS = [
        (0.0,      5.0e-6,   1.0e-6),
        (5.0e-6,   8.0e-6,   2.0e-7),
        (8.0e-6,  15.0e-6,   1.0e-8),
        (15.0e-6, 30.0e-6,   2.0e-7),
    ]

    _print_segment_summary(SEGMENTS)

    steps = make_ramp_segments(SEGMENTS)
    print(f"steps generated = {len(steps)}  ->  set n_steps in input.yaml to {len(steps)}")

    write_bc_yaml(steps)
