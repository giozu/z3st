"""Read an existing static statepoint; never run OpenMC."""
import argparse
import json
from pathlib import Path

import numpy as np
import openmc


def read_results(case_dir):
    manifest = json.loads((case_dir / "manifest.json").read_text())
    path = case_dir / f"statepoint.{manifest['batches']}.h5"
    with openmc.StatePoint(path) as sp:
        total = sp.get_tally(name="B1 integrated fission")
        radial = sp.get_tally(name="B1 radial fission five equal-area rings")
        heating = sp.get_tally(name="factor-for-normalization")
        f = float(total.mean.ravel()[0])
        rings = radial.mean.ravel()
        h = float(heating.mean.ravel()[0])
        if h <= 0 or f <= 0:
            raise ValueError("Positive heating and fission tallies required")
        source_rate = manifest["power_W_for_postprocessing"] / (h * openmc.data.JOULE_PER_EV)
        volumes = np.asarray(manifest["ring_volumes_cm3"])
        relative = (rings / volumes) / (f / volumes.sum())
        result = {
            "case": manifest["case"], "keff": sp.keff.nominal_value,
            "keff_std_dev": sp.keff.std_dev,
            "fission_per_source": f,
            "fission_per_source_std_dev": float(total.std_dev.ravel()[0]),
            "source_rate_per_s": source_rate,
            "fission_per_s_at_250kW": f * source_rate,
            "ring_fission_per_source": rings.tolist(),
            "ring_fission_per_source_std_dev": radial.std_dev.ravel().tolist(),
            "ring_fission_per_s_at_250kW": (rings * source_rate).tolist(),
            "relative_to_FE_mean": relative.tolist(),
            "ring_sum_over_integrated": float(rings.sum() / f),
            "uncertainty_note": "Normalized-rate/profile uncertainties require tally covariance; raw tally standard deviations retained.",
        }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    read_results(parser.parse_args().case_dir)
