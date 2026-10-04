"""Case-local deterministic scale; nominal production correlation unchanged."""
import numpy as np
import dolfinx
import ufl
from z3st.materials.fuel_swelling import uzrh_fission_product_swelling

def applied_swelling_scalar(material,model):
    s=float(material['swelling_sensitivity_factor'])
    if not np.isfinite(s) or not 0<=s<=2:raise ValueError('Sensitivity factor outside verified interval')
    if material.get('fima_source')!='native_openmc':raise ValueError('Native FIMA required')
    space=model.Q
    fn=dolfinx.fem.Function(space,name='Swelling_eigenstrain_applied')
    # Evaluate the actual configured constitutive scalar, not only s*FIMA metadata.
    scalar=scaled_native_swelling(model.T,material,model,3)[0,0]
    if s==0:
        assert float(scalar)==0.0
        fn.x.array[:]=0.0
        fn.x.scatter_forward()
        return fn
    expr=dolfinx.fem.Expression(scalar,space.element.interpolation_points,comm=model.mesh.comm)
    fn.interpolate(expr);fn.x.scatter_forward()
    return fn

def scaled_native_swelling(T,material,model=None,dim=3):
    s=float(material['swelling_sensitivity_factor'])
    if not np.isfinite(s) or not 0<=s<=2:raise ValueError('Sensitivity factor outside [0,2]')
    if material.get('fima_source')!='native_openmc':raise ValueError('Native FIMA required')
    return s*uzrh_fission_product_swelling(T,material,model,dim)
