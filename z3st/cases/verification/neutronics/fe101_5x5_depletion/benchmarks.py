"""Future post-run comparisons; never invoked during preparation."""
def aggregate_history(history, keys):
    import numpy as np
    columns=['volume_cm3','initial_U_Zr_atoms','initial_U_mass_kg','cumulative_fissions',
             'cumulative_deposited_energy_J','fission_rate_per_s','deposited_power_W',
             'U235_atoms','U238_atoms','Pu239_atoms','Cs137_atoms']
    grouped=history.groupby(keys,as_index=False)[columns].sum()
    grouped['FIMA']=grouped.cumulative_fissions/grouped.initial_U_Zr_atoms
    grouped['FIMA_percent']=100*grouped.FIMA
    grouped['BU_MWd_kgU']=grouped.cumulative_deposited_energy_J/(8.64e10*grouped.initial_U_mass_kg)
    weighted=history.assign(weighted_FIMA=history.FIMA*history.initial_U_Zr_atoms,
                            weighted_BU=history.BU_MWd_kgU*history.initial_U_mass_kg)
    independent=weighted.groupby(keys,as_index=False)[['weighted_FIMA','weighted_BU']].sum()
    assert np.allclose(grouped.FIMA,independent.weighted_FIMA/grouped.initial_U_Zr_atoms,rtol=1e-12,atol=1e-18)
    assert np.allclose(grouped.BU_MWd_kgU,independent.weighted_BU/grouped.initial_U_mass_kg,rtol=1e-12,atol=1e-15)
    return grouped


def compare_benchmarks(run_directory, preparation_audit):
    import json
    from pathlib import Path
    import numpy as np
    import pandas as pd
    run=Path(run_directory)
    meta=json.loads(Path(preparation_audit).read_text())
    out=run/'postprocessing'
    history=pd.read_csv(out/'domain_history.csv')
    assert len(history)==525
    rings=aggregate_history(history,['time_h','radial_index'])
    whole=aggregate_history(history,['time_h'])
    axial=aggregate_history(history,['time_h','axial_index'])
    rings.to_csv(out/'radial_aggregated_history.csv',index=False)
    axial.to_csv(out/'axial_aggregated_history.csv',index=False)
    whole.to_csv(out/'FE_mean_history.csv',index=False)
    reference_path=Path(preparation_audit).parent/meta['radial_benchmark']
    reference=pd.DataFrame(json.loads(reference_path.read_text())['history']).rename(columns={'ring':'radial_index'})
    compare=rings.merge(reference,on=['time_h','radial_index'],suffixes=('_2D','_radial'),validate='one_to_one')
    assert len(compare)==105
    for field in ['initial_U_Zr_atoms','initial_U_mass_kg','volume_cm3']:
        assert np.allclose(compare[field+'_2D'],compare[field+'_radial'],rtol=1e-12)
    for field in ['FIMA','BU_MWd_kgU','fission_rate_per_s','deposited_power_W','U235_atoms','U238_atoms','Pu239_atoms','Cs137_atoms']:
        denominator=compare[field+'_radial']
        compare[field+'_relative_difference']=np.divide(compare[field+'_2D']-denominator,denominator,
            out=np.full(len(compare),np.nan),where=denominator.to_numpy()!=0)
    compare.to_csv(out/'benchmark_radial.csv',index=False)
    static_path=Path(preparation_audit).parent/meta['static_benchmark']
    static=json.loads(static_path.read_text())
    static_f=np.asarray(static['matrix_normalized_to_FE_mean'])
    static_power=np.asarray(static['matrix_power_normalized_to_FE_mean'])
    assert static_f.shape==static_power.shape==(5,5)
    initial=history[history.time_h==0].sort_values(['radial_index','axial_index']).copy()
    initial['fission_density']=initial.fission_rate_per_s/initial.volume_cm3
    initial['power_density']=initial.deposited_power_W/initial.volume_cm3
    # FIMA is zero initially; compare its derivative rate/N0, not 0/0.
    initial['initial_FIMA_slope_per_s']=initial.fission_rate_per_s/initial.initial_U_Zr_atoms
    for field in ['fission_density','power_density','initial_FIMA_slope_per_s']:
        initial[field+'_normalized']=initial[field]/np.average(initial[field],weights=initial.volume_cm3)
    initial['static_fission_normalized']=static_f.ravel()
    initial['static_power_normalized']=static_power.ravel()
    initial['fission_pattern_relative_difference']=initial.fission_density_normalized/initial.static_fission_normalized-1
    initial['power_pattern_relative_difference']=initial.power_density_normalized/initial.static_power_normalized-1
    initial['FIMA_slope_pattern_relative_difference']=initial.initial_FIMA_slope_per_s_normalized/initial.static_fission_normalized-1
    initial.to_csv(out/'benchmark_static_initial.csv',index=False)
    summary={'weighted_FE_reconstruction_passed':True,'radial_comparison_rows':len(compare),
             'static_comparison_rows':len(initial),'initial_FIMA_comparison':'Use dFIMA/dt=initial fission rate / fixed initial U+Zr atoms; FIMA(0)=0',
             'interpretation':'Report differences, do not require identity: 2D depletion and Monte Carlo sampling can differ from radial-only results.',
             'uncertainties':'No propagated inter-bin covariance or cumulative uncertainty; differences are diagnostic, not a formal significance test.'}
    (out/'benchmark_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
