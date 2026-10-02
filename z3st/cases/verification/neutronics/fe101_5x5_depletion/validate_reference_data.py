"""Validate compact committed benchmark data; no OpenMC or run outputs needed."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
NEUTRONICS=ROOT.parent

def validate_references():
    radial=NEUTRONICS/'fe101_radial_depletion/reference'
    history=pd.DataFrame(json.loads((radial/'radial_history.json').read_text())['history'])
    summary=json.loads((radial/'radial_summary.json').read_text())
    initial=json.loads((radial/'initial_inventory.json').read_text())
    assert len(history)==105 and not history.isna().any().any()
    assert not history.duplicated(['time_h','ring']).any()
    assert set(history.ring)==set(range(1,6))
    assert np.array_equal(np.sort(history.time_h.unique()),np.arange(21)*175.)
    assert np.allclose(history.FIMA,history.cumulative_fissions/history.initial_U_Zr_atoms,rtol=1e-12,atol=1e-18)
    assert np.allclose(history.BU_MWd_kgU,history.cumulative_deposited_energy_J/(8.64e10*history.initial_U_mass_kg),rtol=1e-12,atol=1e-15)
    for i,ring in enumerate(initial['rings'],1):
        data=history[history.ring==i].sort_values('time_h')
        assert np.allclose(data.volume_cm3,ring['volume_cm3'],rtol=1e-14)
        assert np.allclose(data.initial_U_Zr_atoms,ring['initial_U_Zr_atoms'],rtol=1e-14)
        assert np.allclose(data.initial_U_mass_kg,ring['initial_U_mass_kg'],rtol=1e-14)
        for nuclide in ['U235','U238']:
            assert np.isclose(data.iloc[0][nuclide+'_atoms'],ring['initial_atoms_by_nuclide'][nuclide],rtol=1e-14)
    final=history[history.time_h==3500]
    assert np.isclose(np.average(final.FIMA,weights=final.initial_U_Zr_atoms),summary['FIMA_mean'],rtol=1e-12)
    assert np.isclose(np.average(final.BU_MWd_kgU,weights=final.initial_U_mass_kg),summary['BU_mean_MWd_kgU'],rtol=1e-12)
    static=json.loads((NEUTRONICS/'fe101_axial_static/reference/static_5x5.json').read_text())
    for name in ['matrix_normalized_to_FE_mean','matrix_power_normalized_to_FE_mean']:
        matrix=np.asarray(static[name])
        assert matrix.shape==(5,5) and np.isfinite(matrix).all() and (matrix>0).all()
        assert np.isclose(matrix.mean(),1.,rtol=1e-14)
    axial=np.asarray([r['fission/fission_mean'] for r in static['axial_rows']])
    assert np.allclose(np.asarray(static['matrix_normalized_to_FE_mean']).mean(axis=0),axial,rtol=1e-14)
    geometry=json.loads((NEUTRONICS/'fe101_geometry_static/reference/geometry_results.json').read_text())
    assert len(geometry['cases'])==6
    for case in geometry['cases']:
        assert np.isclose(np.mean(case['relative_to_FE_mean']),1.,rtol=1e-14)
        assert np.isclose(case['ring_sum_over_integrated'],1.,rtol=1e-14)
    print('PASS: compact geometry/static/radial references, indexing, initial inventories and weighted FE FIMA/BU')

if __name__=='__main__':
    validate_references()
