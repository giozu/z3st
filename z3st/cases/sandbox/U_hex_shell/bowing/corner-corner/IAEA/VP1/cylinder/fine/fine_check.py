#!/usr/bin/env python3
"""Fine tube (nt 3, 1 mm layers at the core entry/exit): perimeter-mean hoop
stress (membrane: mean over the perimeter and the three cells through the wall)
against the thin-shell closed form, cell-averaged on each layer."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "..", "..")))
import hex_post as hp

mat = hp.load_material(os.path.join(HERE, "15_15Ti_local.yaml"))
m = hp._import("vp1_temperature", HERE)
d = hp.load_z3st(os.path.join(HERE, "output", "vp1_tube_fine_fields_3d.xdmf"),
                 temperature="vp1_temperature.T_vp1", case_dir=HERE)
Xn = d["X"][d["cells"]]
Xc = Xn.mean(1)
phi = np.arctan2(Xc[:, 1], Xc[:, 0])
tdir = np.stack([-np.sin(phi), np.cos(phi), 0 * phi], 1)
ss = np.einsum("ci,cij,cj->c", tdir, d["S"], tdir)
zz = d["S"][:, 2, 2]
zlo, zhi = Xn[:, :, 2].min(1), Xn[:, :, 2].max(1)
zc = np.round(0.5 * (zlo + zhi), 7)
zu, inv = np.unique(zc, return_inverse=True)
mean_ss = np.bincount(inv, weights=ss) / np.bincount(inv)
lo = np.bincount(inv, weights=zlo) / np.bincount(inv)
hi = np.bincount(inv, weights=zhi) / np.bincount(inv)
E, nu, al = mat["E"], mat["nu"], mat["alpha"]
Rm, t = 0.0716197244, 0.003
beta = (3 * (1 - nu**2)) ** 0.25 / np.sqrt(Rm * t)
s0 = (0.5 * (m.T1_TOP + m.T2_TOP) - m.T_INLET) / (m.Z_CORE_TOP - m.Z_CORE_BOT)
F = lambda u: np.exp(-np.abs(u)) * np.sin(u)
cf = np.zeros_like(zu)
for zk, sk in ((m.Z_CORE_BOT, s0), (m.Z_CORE_TOP, -s0)):
    cf += E * al * sk / (4 * beta) * (F(beta * (hi - zk)) - F(beta * (lo - zk))) / (beta * (hi - lo))
near = np.abs(zu - 1.5) < 0.1
k = np.argmax(np.abs(mean_ss) * near)
ext = lambda v: (float(zu[near][np.abs(v[near]) > 0.1 * np.abs(v[near]).max()].min()),
                 float(zu[near][np.abs(v[near]) > 0.1 * np.abs(v[near]).max()].max()))
above = (zu > 2.8) & (zu < 3.8)
sel_above = above[inv]
out = {"peak_closed_MPa": E * al * s0 / (4 * beta) / 1e6,
       "fe_peak_mean_ss_MPa": float(mean_ss[k] / 1e6), "closed_same_cell_MPa": float(cf[k] / 1e6),
       "z_peak_m": float(zu[k]), "cell_dz_mm": float((hi[k] - lo[k]) * 1e3),
       "rel_err": float(mean_ss[k] / cf[k] - 1),
       "extent_10pct_fe_m": ext(mean_ss), "extent_10pct_closed_m": ext(cf),
       "profile_entry": [(round(float(zu[i]), 4), round(float(mean_ss[i] / 1e6), 4), round(float(cf[i] / 1e6), 4))
                         for i in np.where(near & (np.abs(zu - 1.5) < 0.04))[0][::4]],
       "floor_above_core_MPa": {"max_abs_ss": float(np.abs(ss[sel_above]).max() / 1e6),
                                "max_abs_zz": float(np.abs(zz[sel_above]).max() / 1e6),
                                "max_vm": float(hp.von_mises(d["S"][sel_above]).max() / 1e6)}}
print(json.dumps(out, indent=1))
with open(os.path.join(HERE, "fine_check.json"), "w") as fh:
    json.dump(out, fh, indent=1)
