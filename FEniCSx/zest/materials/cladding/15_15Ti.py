cladding = {
    "name": "15-15 Ti",

    "E": 170e9, # Young's modulus, Pa
    "nu": 0.3,  # Poisson's ratio
    
    "thermal conductivity": {
        "value": 13.95 + 0.01163 * (T-273.15),
        "units": "W/m-K",
        "reference": "Homework 2024-2025, NDT",   
    },
    
    "alpha": 12e-6,  # Thermal expansion coefficient, 1/K, @20°C
    
    "rho": 7750,  # Density, kg/m3
    
    "geometry": {
        "radius_ci": 1.5,  # Inner radius, m
        "radius_co": 2.0,  # Outer radius, m
        "height": 2.0,     # Height, m
    },
}
