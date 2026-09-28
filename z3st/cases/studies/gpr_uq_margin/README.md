# gpr_uq_margin

Propagation of the GPR conductivity posterior through a coupled solve.

The same high-rated MA-MOX pin (45 kW/m, surface held at 650 K) is solved at
`xi = -2 .. +2` posterior standard deviations of the log-residual correction,
`k = k_Magni * exp(mean + xi * sigma)`. `run.py` reports the centre temperature
and the margin to the MOX solidus (Adamson et al., JNM 130 (1985) 349).

The width of the spread is set by the posterior of the *synthetic* fit shipped
with `studies/magni_gpr_conductivity`, not by measured MA-MOX data: what is
demonstrated is the propagation, not the uncertainty of any real correlation.

Run: `./Allrun`
