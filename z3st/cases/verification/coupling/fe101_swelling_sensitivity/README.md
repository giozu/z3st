# FE101 swelling intensity sensitivity

Deterministic parametric sensitivity of the provisional swelling law. No probabilistic UQ or physical correlation validation.

Case-local wrapper sensitivity_law.scaled_native_swelling returns s times the unchanged nominal U-ZrH eigenstrain. The factor is opt-in in isolated fuel cards only: 0, 0.5, 1, 1.5, 2. Native FIMA and power remain unchanged. All five solves run at the 21 saved times 0–3500 h with MPI=1, OMP=2 and preflight resource gate. Existing references and all production Python sources are protected by SHA256 manifests.

Run with PYTHONDONTWRITEBYTECODE=1 using the z3st conda Python: run_sensitivity.py then postprocess.py. Refuses pre-existing output directories. Case-local PYTHONPATH resolves the wrapper without any production-code changes.

The provider field Swelling_eigenstrain_native remains a potential nominal coefficient; new output explicitly distinguishes potential nominal strain from Swelling_eigenstrain_applied.
