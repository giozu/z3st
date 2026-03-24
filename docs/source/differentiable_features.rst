Differentiable Features
=======================

Z3ST exploits the "Differentiable-Native" architecture of UFL to provide Automatic Differentiation (AD) capabilities.

Automatic Differentiation via UFL
--------------------------------
Material laws can be defined in Python as UFL expressions, which are automatically differentiable. 
This is achieved using the ``ufl.derivative()`` operator.

Applications
------------
- **Inverse Modeling**: Identify material parameters from experimental data.
- **Sensitivity Analysis**: Compute the impact of parameter uncertainties on results.
- **Advanced Optimization**: Integrate with gradient-based optimizers for material discovery.
