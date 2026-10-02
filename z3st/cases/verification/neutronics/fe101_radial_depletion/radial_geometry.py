"""Notebook preparation cell: five rings in B1, physical gap, no transport."""
import json

Fuel101.depletable = False
Fuel103.depletable = False
R_FUEL_CM = 1.791
R_CLAD_INNER_CM = 1.804
R_CLAD_OUTER_CM = 1.880
ring_edges_cm = R_FUEL_CM * np.sqrt(np.arange(6) / 5.0)
active_height_cm = top_active_101.z0 - bot_active_101.z0
original_fuel_volume_cm3 = np.pi * R_FUEL_CM**2 * active_height_cm

# Clone the target universe only; shared ordinary FE101 regions are untouched.
fuel_101_universe_B1 = fuel_101_universe.clone(
    clone_materials=False, clone_regions=False)
fuel_101_universe_B1.name = "FE101 B1 physical geometry B - five radial rings"
old_target_cell = next(c for c in fuel_101_universe_B1.cells.values()
                       if c.name == "Fuel 101")
fuel_101_universe_B1.remove_cell(old_target_cell)
active_region = +bot_active_101 & -top_active_101
ring_surfaces = [None] + [openmc.ZCylinder(r=float(r)) for r in ring_edges_cm[1:-1]] + [r_fuel_101]
ring_materials, ring_cells = [], []
for i in range(5):
    material = Fuel101.clone()
    material.name = f"Fuel101_B1_ring_{i+1}"
    material.depletable = True
    material.volume = float(np.pi * (ring_edges_cm[i+1]**2-ring_edges_cm[i]**2)
                            * active_height_cm)
    region = -ring_surfaces[i+1] & active_region
    if i > 0:
        region = region & +ring_surfaces[i]
    cell = openmc.Cell(name=f"FE101 B1 active ring {i+1}", fill=material, region=region)
    ring_materials.append(material)
    ring_cells.append(cell)
    fuel_101_universe_B1.add_cell(cell)
    materials_file.append(material)

target_clad = next(c for c in fuel_101_universe_B1.cells.values()
                   if c.name == "Al Cladding")
physical_clad_inner = openmc.ZCylinder(r=R_CLAD_INNER_CM)
# Preserve all axial/end regions; only active fuel-gap-clad changes.
target_clad.region = target_clad.region & (
    +physical_clad_inner | -bot_active_101 | +top_active_101)
gap_helium = openmc.Material(name="FE101 B1 negligible-density He4 gap")
gap_helium.add_nuclide("He4", 1.0)
gap_helium.set_density("g/cm3", 1.0e-10)
gap_helium.temperature = T_ref_struct
materials_file.append(gap_helium)
target_gap = openmc.Cell(name="FE101 B1 active 130 um He gap", fill=gap_helium,
                        region=+r_fuel_101 & -physical_clad_inner & active_region)
fuel_101_universe_B1.add_cell(target_gap)
B1_pin.fill = fuel_101_universe_B1

# Pure Python geometry/inventory audit. No OpenMC executable or lib.init.
geometry.determine_paths(instances_only=True)
depletable = [m for m in materials_file if m.depletable]
assert len(depletable) == 5
assert {m.id for m in depletable} == {m.id for m in ring_materials}
assert {m.id for m in geometry.get_all_materials().values() if m.depletable} == {m.id for m in ring_materials}
assert Fuel101.num_instances == 60 and not Fuel101.depletable
assert np.isclose(sum(m.volume for m in ring_materials), original_fuel_volume_cm3, rtol=1e-14)

def element_of(nuclide):
    return "".join(c for c in nuclide if c.isalpha())

reference_densities = Fuel101.get_nuclide_atom_densities()
reference_mass_g = Fuel101.get_mass(volume=original_fuel_volume_cm3)
origin = np.asarray(B1_pin.translation)
audit_rings = []
for i, (material, cell) in enumerate(zip(ring_materials, ring_cells)):
    assert material.num_instances == cell.num_instances == 1
    assert material.nuclides == Fuel101.nuclides
    assert material.density == 6.3 and material.density_units == "g/cm3"
    assert material.temperature == Fuel101.temperature
    assert material._sab == Fuel101._sab
    densities = material.get_nuclide_atom_densities()
    assert densities == reference_densities
    # OpenMC atomic densities are atom/b-cm, volume is cm3.
    initial_atoms = {n: float(d * 1e24 * material.volume) for n, d in densities.items()}
    metal_atoms = sum(v for n, v in initial_atoms.items() if element_of(n) in ("U", "Zr"))
    u_mass_kg = sum(material.get_mass(n) for n in densities if element_of(n) == "U") / 1000
    h_zr = (sum(d for n, d in densities.items() if element_of(n) == "H") /
            sum(d for n, d in densities.items() if element_of(n) == "Zr"))
    for phi in (0.0, 0.7, 2.0, 4.0):
        for z in (bot_active_101.z0+1e-6, (bot_active_101.z0+top_active_101.z0)/2,
                  top_active_101.z0-1e-6):
            radius = np.sqrt((ring_edges_cm[i]**2+ring_edges_cm[i+1]**2)/2)
            point = origin + [radius*np.cos(phi), radius*np.sin(phi), z]
            assert geometry.find(tuple(point))[-1].fill is material
    audit_rings.append({
        "ring": i+1, "material_id": material.id, "cell_id": cell.id,
        "instances": material.num_instances,
        "r_in_cm": float(ring_edges_cm[i]), "r_out_cm": float(ring_edges_cm[i+1]),
        "r_in_over_R": float(ring_edges_cm[i]/R_FUEL_CM),
        "r_out_over_R": float(ring_edges_cm[i+1]/R_FUEL_CM),
        "volume_cm3": material.volume, "density_g_cm3": material.density,
        "mass_g": material.get_mass(), "initial_U_mass_kg": u_mass_kg,
        "initial_U_Zr_atoms": metal_atoms, "initial_atoms_by_nuclide": initial_atoms,
        "H_Zr_atom_ratio": h_zr, "temperature_K": material.temperature,
    })
assert np.isclose(sum(m.get_mass() for m in ring_materials), reference_mass_g, rtol=1e-14)
for nuc, density in reference_densities.items():
    actual = sum(r["initial_atoms_by_nuclide"][nuc] for r in audit_rings)
    assert np.isclose(actual, density*1e24*original_fuel_volume_cm3, rtol=1e-14)
for radius, expected in [(R_FUEL_CM-1e-6, ring_materials[-1]),
                         (R_FUEL_CM+1e-6, gap_helium),
                         (R_CLAD_INNER_CM-1e-6, gap_helium),
                         (R_CLAD_INNER_CM+1e-6, CladdingAl)]:
    point = origin + [radius, 0, (bot_active_101.z0+top_active_101.z0)/2]
    assert geometry.find(tuple(point))[-1].fill is expected
# Ordinary FE101 cell still has the original fuel/clad regions.
assert fuel_101_cell.fill is Fuel101
assert clad_101_cell.fill is CladdingAl

audit = {
    "geometry_baseline": "physical B; legacy A is not the radial baseline",
    "gap_radial_um": (R_CLAD_INNER_CM-R_FUEL_CM)*1e4,
    "gap_nuclide": "He4", "gap_density_g_cm3": 1e-10,
    "active_z_cm": [bot_active_101.z0, top_active_101.z0],
    "B1_origin_cm": origin.tolist(), "original_volume_cm3": original_fuel_volume_cm3,
    "reconstructed_volume_cm3": sum(m.volume for m in ring_materials),
    "original_mass_g": reference_mass_g,
    "reconstructed_mass_g": sum(m.get_mass() for m in ring_materials),
    "ordinary_FE101_instances": Fuel101.num_instances,
    "depletable_materials": len(depletable), "rings": audit_rings,
    "audit_method": "Python instance counts, analytic volumes and point locations; no transport",
}
print("PASS: five distinct materials, one instance each; 60 ordinary FE101 unchanged")
print("PASS: volumes, mass and every initial nuclide inventory reconstruct the original fuel")
print("PASS: physical 130 um He gap; unchanged axial boundaries")
print(pd.DataFrame(audit_rings)[["ring", "volume_cm3", "mass_g", "initial_U_Zr_atoms",
                               "H_Zr_atom_ratio", "instances"]].to_string(index=False))
