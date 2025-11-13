Examples
========

Simple steady-state case
------------------------

1. Define the input in ``input.yaml``:
   .. code-block:: yaml

      lhr:
      - 0
      time:
      - 0
      n_steps: 1

2. Run the simulation:
   .. code-block:: bash

        python z3st.py

3. Visualize results:
   .. code-block:: bash

        paraview output/fields.vtu
