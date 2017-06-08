"""
Docstring.
"""
from numpy import ones

#==========================================================================
def properties_for_water(rc):
    """
    Returns physical properties for water for given resolution 'rc'

    For 60 deg from:
    http://www.engineeringtoolbox.com/water-properties-d_1508.html
    -------------------------------------------------------------------------
    """
    # Create and fill matrices for all properties
    rho = ones(rc) * 983.0       # density              [kg/m^3]
    mu = ones(rc) * 0.466E-3     # viscosity            [Pa s]
    cp = ones(rc) * 4185         # thermal capacity     [J/kg/K]
    kappa = ones(rc) * 0.654     # thermal conductivity [W/m/K]

    return rho, mu, cp, kappa  # end of function
