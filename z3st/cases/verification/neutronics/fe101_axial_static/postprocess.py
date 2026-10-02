"""Check conservation and report static axial and radial/axial tally profiles."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import openmc
ROOT=Path(__file__).resolve().parent
m=json.loads((ROOT/'manifest.json').read_text())
v=np.asarray(m['bin_volumes_cm3'])
r=np.asarray(m['radial_edges_cm'])
z=np.asarray(m['axial_edges_cm'])
with openmc.StatePoint(ROOT/'run_preliminary/statepoint.30.h5') as sp:
    integrated=sp.get_tally(name='B1 static integrated').mean.reshape(5,2).sum(axis=0)
    axial_t=sp.get_tally(name='B1 static axial')
    grid_t=sp.get_tally(name='B1 static radial axial')
    axial=axial_t.mean.reshape(5,5,2).sum(axis=0)
    # Mesh bins are (r,phi,z); explicitly map them rather than assuming flattening order.
    bins=list(grid_t.filters[1].bins)
    raw=grid_t.mean.reshape(5,25,2).sum(axis=0)
    raw_std=grid_t.std_dev.reshape(5,25,2)
    grid=np.zeros((5,5,2))
    for b,values in zip(bins,raw):
        ir,iphi,iz=b
        assert iphi==1
        grid[ir-1,iz-1]=values
    global_heat=float(sp.get_tally(name='factor-for-normalization').mean.ravel()[0])
    source_per_s=m['power_W_for_density']/(global_heat*openmc.data.JOULE_PER_EV)
    k={'mean':float(sp.keff.nominal_value),'std_dev':float(sp.keff.std_dev)}
ax_error=(axial.sum(axis=0)-integrated)/integrated
grid_error=(grid.sum(axis=(0,1))-integrated)/integrated
assert np.allclose(axial.sum(axis=0),integrated,rtol=1e-10,atol=0)
assert np.allclose(grid.sum(axis=(0,1)),integrated,rtol=1e-10,atol=0)
assert np.allclose(grid.sum(axis=0),axial,rtol=1e-10,atol=0)
ax_v=v.sum(axis=0)
fd=axial[:,0]*source_per_s/ax_v
pdensity=axial[:,1]*source_per_s*openmc.data.JOULE_PER_EV/ax_v
fd_mean=integrated[0]*source_per_s/v.sum()
pd_mean=integrated[1]*source_per_s*openmc.data.JOULE_PER_EV/v.sum()
norm=grid[:,:,0]*source_per_s/v/fd_mean
local_norm=norm/norm.mean(axis=0)[None,:]
ratio=norm[-1]/norm[0]
rows=pd.DataFrame({'axial_zone':np.arange(1,6),'z_low_cm':z[:-1],'z_high_cm':z[1:],'volume_cm3':ax_v,'fission_tally_per_source':axial[:,0],'heating_local_eV_per_source':axial[:,1],'fission_density_per_s_cm3':fd,'power_density_W_cm3':pdensity,'fission/fission_mean':fd/fd_mean,'power/power_mean':pdensity/pd_mean,'outer_ring/inner_ring':ratio})
rows.to_csv(ROOT/'axial_profile.csv',index=False)
pd.DataFrame(norm,index=np.arange(1,6),columns=[f'z{j}' for j in range(1,6)]).rename_axis('radial_ring').to_csv(ROOT/'radial_axial_normalized.csv')
pd.DataFrame(local_norm,index=np.arange(1,6),columns=[f'z{j}' for j in range(1,6)]).rename_axis('radial_ring').to_csv(ROOT/'radial_profile_at_z.csv')
end_fd=(fd[0]+fd[-1])/2
end_pd=(pdensity[0]+pdensity[-1])/2
summary={'execution':json.loads((ROOT/'execution_status.json').read_text()),'keff':k,'power_normalization_W':m['power_W_for_density'],'FE_volume_cm3':float(v.sum()),'zone_volume_cm3':float(ax_v[0]),'bin_volume_cm3':float(v[0,0]),'FE_fission_rate_per_s':float(integrated[0]*source_per_s),'FE_deposited_power_W':float(integrated[1]*source_per_s*openmc.data.JOULE_PER_EV),'FE_mean_fission_density_per_s_cm3':float(fd_mean),'FE_mean_power_density_W_cm3':float(pd_mean),'axial_sum_relative_errors_fission_heating':ax_error.tolist(),'grid_sum_relative_errors_fission_heating':grid_error.tolist(),'axial_fission_max_min_ratio':float(fd.max()/fd.min()),'axial_power_max_min_ratio':float(pdensity.max()/pdensity.min()),'center_above_mean_end_fission_percent':float(100*(fd[2]/end_fd-1)),'center_above_mean_end_power_percent':float(100*(pdensity[2]/end_pd-1)),'center_above_each_end_fission_percent':(100*(fd[2]/fd[[0,4]]-1)).tolist(),'outer_inner_ratio_range':float(ratio.max()-ratio.min()),'outer_inner_ratio_relative_span_percent':float(100*(ratio.max()-ratio.min())/ratio.mean()),'axial_rows':rows.to_dict(orient='records'),'matrix_normalized_to_FE_mean':norm.tolist(),'radial_profile_normalized_at_each_z':local_norm.tolist(),'mesh_bin_order':[list(map(int,b)) for b in bins], 'uncertainty_note':'No covariance between tally bins or global normalization propagated; normalized-profile significance is not a formal hypothesis test.'}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
