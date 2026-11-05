# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# beam1d: Euler-Bernoulli model of the hexagonal wrapper
# Version: 0.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
A beam1d case directory, laid out as a z3st case:

    input.yaml                 geometry_path, boundary_conditions_path,
                               materials: {name: card}, models, section, element
    geometry.yaml              D, t, H, orientation, nz, stations {name: z (m)}
    boundary_conditions.yaml   thermal: {name: [T entry]}  (beam1d.thermal)
                               mechanical: {name: [entries]}

Stations play the role of the gmsh labels: every station is a node of the
mesh and a region for the boundary conditions; bottom (z = 0) and top (z = H)
always exist. Mechanical entries, components in the global axes:

    Dirichlet  region, displacement: [u_x, u_y, u_z], rotation: [u_x', u_y']
               (null leaves a component free)
    Neumann    region, force: [F_x, F_y, F_z] (N), moment: [m_x, m_y] (N m),
               m_a conjugate to u_a'
    Contact    region, direction: x | y, gap (m), side: +1 | -1
               (unilateral restraint at side * gap, beam1d.solver)
    Body       gravity: [g_x, g_y, g_z] (m/s^2), load rho A_w g per unit length
"""

import os
from types import SimpleNamespace

import yaml

from beam1d.section import HexSection

def load_yaml(case_dir, name):
    with open(os.path.join(case_dir, name)) as f:
        return yaml.safe_load(f)


def load_case(case_dir="."):
    """The case files of case_dir, the section and the stations."""
    case_dir = os.path.abspath(case_dir)
    inp = load_yaml(case_dir, "input.yaml")
    geo = load_yaml(case_dir, inp.get("geometry_path", "geometry.yaml"))
    bcs = load_yaml(case_dir, inp.get("boundary_conditions_path", "boundary_conditions.yaml"))
    mat_name, mat_path = next(iter(inp["materials"].items()))
    mat = load_yaml(case_dir, mat_path)
    H = float(geo["H"])
    stations = {"bottom": 0.0, "top": H}
    stations.update({k: float(v) for k, v in (geo.get("stations") or {}).items()})
    sec = HexSection(float(geo["D"]), float(geo["t"]), geo.get("orientation", "flat"),
                     inp.get("section", {}).get("n_gauss", 4))
    return SimpleNamespace(dir=case_dir, inp=inp, geo=geo, bcs=bcs, mat_name=mat_name,
                           mat_path=os.path.normpath(os.path.join(case_dir, mat_path)), mat=mat,
                           H=H, stations=stations, section=sec)


def mechanical(case, kind):
    """The mechanical entries of one type (Dirichlet, Neumann, Contact, Body)."""
    return [bc for bc in (case.bcs.get("mechanical", {}).get(case.mat_name) or []) if bc.get("type") == kind]


def station(case, region):
    if region not in case.stations:
        raise KeyError(f"region '{region}' is not a station of geometry.yaml ({sorted(case.stations)})")
    return case.stations[region]
