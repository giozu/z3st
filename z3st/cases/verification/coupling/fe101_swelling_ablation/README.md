# FE101 swelling ablation

Numerical / model-effect isolation, not physical validation of the provisional U-ZrH swelling correlation.

Reference ON: ../fe101_native_power/runs/full_native (read only). OFF imports the identical native FIMA and power histories; only the material eigenstrain callable is omitted. All other material parameters, boundary conditions, mesh, time horizon and solver settings are unchanged. No source/solver modifications.

Use the 21 saved times, 0 to 3500 h. One MPI process, two OpenMP threads, BLAS/NumExpr one thread. Resources and immutable-file manifests accompany the results. Global-average gap/contact limitation is retained.

The imported provider's Swelling_eigenstrain_native field is a potential coefficient, not the applied constitutive strain in OFF. OFF VTU outputs explicitly rename it as potential and include a separate zero applied-swelling field. FIMA is never replaced with zero.

Reproduction: /home/simone/miniconda3/envs/z3st/bin/python run_ablation.py; then postprocess.py. The runner refuses existing outputs.
