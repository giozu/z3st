# High-temperature phase fraction y_s(T, H) for Zircaloy (OPERA-HPC D2.1, Eq. 11)

import ufl


def y_beta(T, material):
    """Volume fraction of beta-phase y_s in Zircaloy during alpha -> beta transition.
    T (temperature)         [K]
    H (hydrogen content)    [wppm]
    """
    H = float(material.get("hydrogen_content", 0.0))
    T_mid = 1159.0 - 0.096 * H
    delta_T = 44.0 - 0.026 * H
    
    # return 0.5 * (1.0 + ufl.tanh((T - T_mid) / delta_T))
    return 0.58