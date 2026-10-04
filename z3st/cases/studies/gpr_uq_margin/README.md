# gpr_uq_margin

Propagation of the GPR conductivity posterior through a coupled solve.

The same high-rated MA-MOX pin (45 kW/m, surface held at 650 K) is solved at
`xi = -2 .. +2` posterior standard deviations of the log-residual correction,
`k = k_Magni * exp(mean + xi * sigma)`. `run.py` reports the centre temperature
and the margin to the MOX solidus (Adamson et al., JNM 130 (1985) 349).

The card uses `mode: affine`. The width of the spread is set by the posterior of
the synthetic fit built by `studies/magni_gpr_conductivity/make_synthetic_gpr.py`
(called by `Allrun`). No measured MA-MOX data enter. The case shows the
propagation of a posterior, not the uncertainty of a real correlation.

Not in the suite (no gold). Outputs: `output/uq_margin.png`, `output/uq_margin.npz`.

Run: `./Allrun`
