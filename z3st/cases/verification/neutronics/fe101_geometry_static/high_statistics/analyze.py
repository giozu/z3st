"""Augment saved reader JSON and compare static runs; no transport/depletion."""
import json
import math
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parent


def ratio_error(numerator, denominator, sn, sd):
    ratio = numerator / denominator
    a, b = sn / numerator, sd / denominator
    return ratio * math.hypot(a, b), ratio * abs(a-b), ratio * (a+b)


def augment(case):
    directory = ROOT / f"case_{case}"
    result = json.loads((directory / "results.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    with h5py.File(directory / "statepoint.220.h5") as state:
        tally = state['tallies/tally 3']
        n = int(tally['n_realizations'][()])
        sums = tally['results'][:].reshape(-1, 2)[0]
        mean = sums[0] / n
        std = math.sqrt(max(0, (sums[1]/n-mean**2)/(n-1)))
        assert int(state['seed'][()]) == manifest['seed']
        assert n == 200
    result['heating_eV_per_source'] = float(mean)
    result['heating_eV_per_source_std_dev'] = std
    F, sF = result['fission_per_source'], result['fission_per_source_std_dev']
    rate = result['fission_per_s_at_250kW']
    errors = ratio_error(F, mean, sF, std)
    result['fission_rate_std_dev_covariance_zero'] = errors[0] * rate / (F/mean)
    result['fission_rate_std_dev_covariance_bounds'] = [e * rate/(F/mean) for e in errors[1:]]
    bins = []
    for i, volume in enumerate(manifest['ring_volumes_cm3']):
        f, s = result['ring_fission_per_source'][i], result['ring_fission_per_source_std_dev'][i]
        error, lower, upper = ratio_error(f, F, s, sF)
        bins.append({
            'ring': i+1,
            'r_in_over_R': manifest['normalized_radial_edges'][i],
            'r_out_over_R': manifest['normalized_radial_edges'][i+1],
            'r_in_cm': manifest['radial_edges_cm'][i],
            'r_out_cm': manifest['radial_edges_cm'][i+1],
            'volume_cm3': volume, 'fission_per_source': f,
            'fission_per_source_std_dev': s,
            'density_per_source_cm3': f / volume,
            'density_per_source_cm3_std_dev': s / volume,
            'profile': 5*f/F,
            'profile_std_dev_covariance_zero': 5*error,
            'profile_std_dev_covariance_bounds': [5*lower, 5*upper],
        })
    result['rings'] = bins
    result['seed'] = manifest['seed']
    result['exit_code'] = int((directory/'exit_code.txt').read_text())
    result['uncertainty_note'] = (
        'Raw tally errors are OpenMC standard errors. Ratio errors use first-order '
        'propagation; covariance-zero estimates are approximate. Covariance bounds '
        'use |Cov(X,Y)|<=sigma_X*sigma_Y; they are not confidence intervals. '
        'Distinct seeds permit independent-run propagation; batch correlations '
        'and source convergence are not independently assessed.')
    (directory/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def percent_comparison(x, b, sx, sb):
    ratio = x / b
    sigma = 100 * ratio * math.hypot(sx/x, sb/b)
    delta = 100 * (ratio-1)
    return {'delta_percent': delta, 'sigma_percentage_points': sigma,
            'abs_delta_over_sigma': abs(delta)/sigma}


def main():
    results = {case: augment(case) for case in 'ABCD'}
    baseline = results['B']
    comparisons = {}
    for case in 'ACD':
        x = results[case]
        dk = x['keff']-baseline['keff']
        sk = math.hypot(x['keff_std_dev'], baseline['keff_std_dev'])
        comparison = {'delta_keff': dk, 'delta_keff_pcm': dk*1e5,
                      'sigma_keff': sk, 'sigma_keff_pcm': sk*1e5,
                      'abs_delta_keff_over_sigma': abs(dk)/sk,
                      'fission_tally': percent_comparison(
                          x['fission_per_source'], baseline['fission_per_source'],
                          x['fission_per_source_std_dev'], baseline['fission_per_source_std_dev']),
                      'fission_rate_covariance_zero': percent_comparison(
                          x['fission_per_s_at_250kW'], baseline['fission_per_s_at_250kW'],
                          x['fission_rate_std_dev_covariance_zero'], baseline['fission_rate_std_dev_covariance_zero'])}
        comparison['rings'] = []
        for xr, br in zip(x['rings'], baseline['rings']):
            entry = {'ring': xr['ring'],
                     'raw_fission': percent_comparison(
                         xr['fission_per_source'], br['fission_per_source'],
                         xr['fission_per_source_std_dev'], br['fission_per_source_std_dev']),
                     'density': percent_comparison(
                         xr['density_per_source_cm3'], br['density_per_source_cm3'],
                         xr['density_per_source_cm3_std_dev'], br['density_per_source_cm3_std_dev']),
                     'profile_covariance_zero': percent_comparison(
                         xr['profile'], br['profile'], xr['profile_std_dev_covariance_zero'],
                         br['profile_std_dev_covariance_zero'])}
            entry['profile_sigma_percentage_points_covariance_bounds'] = [
                percent_comparison(xr['profile'], br['profile'],
                                   xr['profile_std_dev_covariance_bounds'][j],
                                   br['profile_std_dev_covariance_bounds'][j])['sigma_percentage_points']
                for j in range(2)]
            comparison['rings'].append(entry)
        comparison['fission_rate_sigma_percentage_points_covariance_bounds'] = [
            percent_comparison(x['fission_per_s_at_250kW'], baseline['fission_per_s_at_250kW'],
                               x['fission_rate_std_dev_covariance_bounds'][j],
                               baseline['fission_rate_std_dev_covariance_bounds'][j])['sigma_percentage_points']
            for j in range(2)]
        comparisons[f'{case}-B'] = comparison
    (ROOT/'comparison.json').write_text(json.dumps(comparisons, indent=2)+'\n')
    print(json.dumps(comparisons, indent=2))


if __name__ == '__main__':
    main()
