"""Export one static sensitivity case; never launch transport or depletion."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy

import numpy as np
import openmc


def prepare(case, high_statistics=False):
    root = Path(__file__).resolve().parent
    ns = runpy.run_path(str(root / "source_model.py"))
    geometry = ns["geometry"]
    materials = ns["materials_file"]
    settings = ns["settings_file"]
    tallies = ns["tallies_file"]
    target = ns["fuel_101_universe_B1"]
    fuel = next(c for c in target.cells.values() if c.name == "Fuel 101")
    clad = next(c for c in target.cells.values() if c.name == "Al Cladding")
    assert fuel.fill is ns["Fuel101_single_B1"]
    assert ns["B1_pin"].fill is target
    r_fuel, r_outer = ns["r_fuel_101"], ns["r_clad_101"]
    bottom, top = ns["bot_active_101"], ns["top_active_101"]
    active = +bottom & -top
    original_clad_region = clad.region

    initial_radius = r_fuel.r
    radius = {"A": 1.791, "B": 1.791, "C": 1.7975, "D": 1.804}[case]
    density = 6.3 * (initial_radius / radius)**2
    if case in ("C", "D"):
        r_fuel = openmc.ZCylinder(r=radius, name="B1 expanded active fuel outer")
        fuel.region = -r_fuel & active
        fuel.fill.set_density("g/cm3", density)
        fuel.fill.volume = np.pi * radius**2 * (top.z0 - bottom.z0)

    if case in ("B", "C", "D"):
        inner = openmc.ZCylinder(r=1.804, name="B1 physical clad inner")
        # D uses one shared surface for exact fuel-clad contact.
        if case == "D":
            inner = r_fuel
        clad.region = (original_clad_region & ~active) | (+inner & -r_outer & active)
    if case in ("B", "C"):
        helium = openmc.Material(name="B1 negligible-density He4 gap")
        helium.add_nuclide("He4", 1.0)
        helium.set_density("g/cm3", 1.0e-10)
        helium.temperature = ns["T_ref_struct"]
        materials.append(helium)
        # Change only the target's active fuel-gap-clad segment. Preserve
        # all axial reflector, poison-disk and end-cladding geometry.
        gap = openmc.Cell(name="B1 active He gap", fill=helium,
                          region=+r_fuel & -inner & active)
        target.add_cell(gap)

    # Mesh tally subdivides measurement only; fuel composition/cell is unchanged.
    edges = r_fuel.r * np.sqrt(np.arange(6) / 5.0)
    origin = list(ns["B1_pin"].translation)
    mesh = openmc.CylindricalMesh(
        r_grid=edges, z_grid=[bottom.z0 - origin[2], top.z0 - origin[2]],
        origin=origin, name="B1 five equal-area measurement rings")
    integrated = openmc.Tally(name="B1 integrated fission")
    target_filter = openmc.CellFilter(fuel)
    integrated.filters = [target_filter]
    integrated.scores = ["fission"]
    integrated.estimator = "tracklength"
    radial = openmc.Tally(name="B1 radial fission five equal-area rings")
    radial.filters = [target_filter, openmc.MeshFilter(mesh)]
    radial.scores = ["fission"]
    radial.estimator = "tracklength"
    tallies.extend([integrated, radial])

    if high_statistics:
        settings.particles = 5000
        settings.inactive = 20
        settings.batches = 220
        settings.seed = {"A": 101, "B": 202, "C": 303, "D": 404}[case]
    out = root / ("high_statistics" if high_statistics else ".") / f"case_{case}"
    out.mkdir(parents=True, exist_ok=True)
    geometry.export_to_xml(out / "geometry.xml")
    materials.export_to_xml(out / "materials.xml")
    settings.export_to_xml(out / "settings.xml")
    tallies.export_to_xml(out / "tallies.xml")
    # Python-side point-location checks: no OpenMC executable/library init.
    x0, y0, z0 = origin
    z = (bottom.z0 + top.z0) / 2 + z0
    for r, expected in [(1.70, fuel.fill), (1.85, ns["CladdingAl"])]:
        assert geometry.find((x0 + r, y0, z))[-1].fill is expected
    if case in ("B", "C"):
        assert geometry.find((x0 + (radius + 1.804)/2, y0, z))[-1].fill is helium
    assert geometry.find((x0 + radius - 1e-6, y0, z))[-1].fill is fuel.fill
    expected = ns["CladdingAl"] if case in ("A", "D") else helium
    assert geometry.find((x0 + radius + 1e-6, y0, z))[-1].fill is expected
    assert np.isclose(density * fuel.fill.volume,
                      6.3 * np.pi * initial_radius**2 * (top.z0-bottom.z0), rtol=1e-14)
    # The ordinary FE101 universe remains intact in case B.
    assert ns["clad_101_cell"].region is original_clad_region
    xs = Path(ns["xs"])
    assert xs.is_file()
    manifest = {
        "case": case, "power_W_for_postprocessing": 250000.0,
        "fuel_outer_cm": r_fuel.r, "clad_inner_cm": 1.791 if case == "A" else 1.804,
        "clad_outer_cm": r_outer.r, "active_z_cm": [bottom.z0, top.z0],
        "B1_origin_cm": origin, "target_fuel_cell_id": fuel.id,
        "radial_edges_cm": edges.tolist(),
        "ring_volumes_cm3": (np.pi * np.diff(edges**2) * (top.z0-bottom.z0)).tolist(),
        "fuel_density_g_cm3": density, "fuel_mass_g": density * fuel.fill.volume,
        "normalized_radial_edges": (edges/radius).tolist(),
        "gap_radial_um": 0.0 if case == "A" else (1.804-radius)*1e4,
        "gap": None if case in ("A", "D") else {"nuclide": "He4", "density_g_cm3": 1e-10,
                                              "temperature_K": helium.temperature},
        "cross_sections": str(xs), "xs_xml_sha256": hashlib.sha256(xs.read_bytes()).hexdigest(),
        "batches": settings.batches, "inactive": settings.inactive,
        "particles": settings.particles,
        "seed": settings.seed,
        "transport_executed": False,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Exported {out}; no transport/depletion executed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["A", "B", "C", "D"], required=True)
    parser.add_argument("--high-statistics", action="store_true")
    args = parser.parse_args()
    prepare(args.case, args.high_statistics)
