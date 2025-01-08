cladding = {
    "name": "zircaloy-4",

    "E": 99.3e9,  # Young's modulus, Pa
    "nu": 0.37,  # Poisson's ratio
    
    "thermal conductivity": {
        "value": 21.5,
        "units": "W/m-K",
        "reference": "https://www.azom.com/article.aspx?ArticleID=7644",   
    },
    
    "melting temperature" : 2123,  # K

    "alpha": 6e-6,  # Thermal expansion coefficient, 1/K, @20°C

    "rho": 6560,  # Density, kg/m3
        
    "geometry": {
        "radius_ci": 1.5,  # Inner radius, m
        "radius_co": 2.0,  # Outer radius, m
        "height": 2.0,     # Height, m
    },
}
