cladding = {
    "name": "T91",

    "E": 170e9, # Young's modulus, Pa
    "nu": 0.3,  # Poisson's ratio
    
    "thermal conductivity": {
        "value": 28.7,
        "units": "W/m-K",
        "reference": "handouts",   
    },
    
    "alpha": 12e-6,  # Thermal expansion coefficient, 1/K, @20°C
    
    "rho": 7750,  # Density, kg/m3
    
    "geometry": {
        "radius_ci": 1.5,  # Inner radius, m
        "radius_co": 2.0,  # Outer radius, m
        "height": 2.0,     # Height, m
    },
}
