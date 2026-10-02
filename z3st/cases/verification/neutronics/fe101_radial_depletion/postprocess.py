"""Prepared postprocessing; runs only on explicit invocation after depletion."""
def postprocess_radial(run_directory, preparation_audit):
    import json
    from pathlib import Path
    import numpy as np
    import pandas as pd
    import openmc
    import openmc.deplete

    run_directory = Path(run_directory)
    metadata = json.loads(Path(preparation_audit).read_text())
    results_path = run_directory / "depletion_results_FE101_B1_1965_5RINGS_3500h_ED.h5"
    results = openmc.deplete.Results(results_path)
    times = np.asarray(results.get_times(time_units="s"), dtype=float)
    assert len(times) == 21 and np.all(np.diff(times) > 0)
    assert np.isclose(times[-1], 3500*3600)
    ring_metadata = metadata["rings"]
    assert set(results[0].index_mat) == {str(r["material_id"]) for r in ring_metadata}
    # Read the same chain used by the run; include fast-fission contributors,
    # e.g. U238, and actinides formed during irradiation, not only U235/Pu239.
    chain = openmc.deplete.Chain.from_xml(metadata["chain_file"])
    fission_nuclides = [n.name for n in chain.nuclides
                       if any(reaction.type == "fission" for reaction in n.reactions)]
    with openmc.StatePoint(run_directory / "openmc_simulation_n0.h5") as sp:
        score = sp.get_tally(name="FE101 B1 ring energy deposition").scores[0]
        assert score == "heating-local"
    power_by_ring = np.zeros((len(times), 5))
    tally_fission_by_ring = np.zeros_like(power_by_ring)
    # Predictor saves BOS statepoints n0..n19 and the final transport n20.
    for k in range(len(times)):
        path = run_directory / f"openmc_simulation_n{k}.h5"
        with openmc.StatePoint(path) as sp:
            global_heat = float(sp.get_tally(name="factor-for-normalization").mean.ravel()[0])
            if global_heat <= 0:
                raise ValueError(f"Non-positive global heating in {path}")
            ring_heat = sp.get_tally(name="FE101 B1 ring energy deposition")
            ring_fission = sp.get_tally(name="FE101 B1 ring total fission")
            heat_ids = list(ring_heat.filters[0].bins)
            fission_ids = list(ring_fission.filters[0].bins)
            source_per_s = 250000.0 / (global_heat * openmc.data.JOULE_PER_EV)
            for i, ring in enumerate(ring_metadata):
                mid = ring["material_id"]
                power_by_ring[k, i] = 250000.0 * ring_heat.mean[heat_ids.index(mid), 0, 0] / global_heat
                tally_fission_by_ring[k, i] = ring_fission.mean[fission_ids.index(mid), 0, 0] * source_per_s

    rows, nuclide_rate_rows, zest_rows = [], [], []
    for i, ring in enumerate(ring_metadata):
        mid = str(ring["material_id"])
        total_rate = np.zeros(len(times))
        for nuc in fission_nuclides:
            if nuc not in results[0].rates.index_nuc:
                if nuc in results[0].index_nuc:
                    _, atoms = results.get_atoms(mid, nuc)
                    if np.any(atoms > 0):
                        raise ValueError(f"Fissioning nuclide {nuc} has inventory but no saved rate")
                continue
            t, rates = results.get_reaction_rate(mid, nuc, "fission")
            assert np.allclose(t, times)
            total_rate += rates
            for k, rate in enumerate(rates):
                nuclide_rate_rows.append({"time_s": times[k], "ring": i+1,
                                          "nuclide": nuc, "fission_rate_per_s": float(rate)})
        inventory = {}
        for nuc in ("U235", "U238", "Pu239", "Cs137"):
            t, atoms = results.get_atoms(mid, nuc, nuc_units="atoms")
            assert np.allclose(t, times)
            inventory[nuc] = atoms
        for nuc in ("U235", "U238"):
            assert np.isclose(inventory[nuc][0], ring["initial_atoms_by_nuclide"][nuc], rtol=1e-10)
        # First-order left-endpoint integration, matching Predictor's BOS
        # transport sampling. Not an exact Bateman-integrated fission count.
        # Final_step=True supplies the genuine final rate but no extra interval.
        cumulative = np.r_[0.0, np.cumsum(total_rate[:-1] * np.diff(times))]
        energy_j = np.r_[0.0, np.cumsum(power_by_ring[:-1, i] * np.diff(times))]
        fima = cumulative / ring["initial_U_Zr_atoms"]
        burnup = energy_j / (8.64e10 * ring["initial_U_mass_kg"])
        for k, t in enumerate(times):
            row = {
                "time_s": t, "time_h": t/3600, "ring": i+1, "material_id": int(mid),
                "r_in_cm": ring["r_in_cm"], "r_out_cm": ring["r_out_cm"],
                "r_in_over_R": ring["r_in_over_R"], "r_out_over_R": ring["r_out_over_R"],
                "volume_cm3": ring["volume_cm3"],
                "initial_U_Zr_atoms": ring["initial_U_Zr_atoms"],
                "initial_U_mass_kg": ring["initial_U_mass_kg"],
                **{f"{n}_atoms": float(a[k]) for n, a in inventory.items()},
                "fission_rate_per_s": float(total_rate[k]),
                "tally_fission_rate_per_s": float(tally_fission_by_ring[k, i]),
                "cumulative_fissions": float(cumulative[k]),
                "FIMA": float(fima[k]), "FIMA_percent": float(100*fima[k]),
                "deposited_power_W": float(power_by_ring[k, i]),
                "cumulative_deposited_energy_J": float(energy_j[k]),
                "BU_MWd_kgU": float(burnup[k]),
            }
            rows.append(row)
            zest_rows.append({key: row[key] for key in ("time_s", "time_h", "ring", "r_in_cm", "r_out_cm",
                             "r_in_over_R", "r_out_over_R", "FIMA")})
    output = run_directory / "postprocessing"
    output.mkdir(exist_ok=True)
    pd.DataFrame(rows).sort_values(["time_s", "ring"]).to_csv(output/"ring_history.csv", index=False)
    pd.DataFrame(nuclide_rate_rows).to_csv(output/"fission_rates_by_nuclide.csv", index=False)
    pd.DataFrame(zest_rows).sort_values(["time_s", "ring"]).to_csv(output/"zest_FIMA.csv", index=False)
    payload = {
        "geometry_baseline": "physical B", "primary_ZEST_variable": "FIMA (fraction)",
        "fima_definition": "cumulative fissions / fixed initial U+Zr atoms; H excluded",
        "integration": "first-order BOS/left-endpoint; time-discretization approximation",
        "burnup_definition": "energy deposited in ring / initial uranium mass; MWd/kgU",
        "energy_normalization": "250 kW whole-reactor heating-local; no fixed energy/fission",
        "uncertainties": "Monte Carlo/depletion/time-integration errors not propagated into cumulative outputs",
        "rings": ring_metadata, "history": sorted(rows, key=lambda r: (r["time_s"], r["ring"])),
    }
    (output/"ring_history.json").write_text(json.dumps(payload, indent=2)+"\n")
    print(f"Saved {len(rows)} ring/time rows to {output}")
    return pd.DataFrame(rows)
