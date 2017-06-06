"""
Standard Python modules.
"""
from numpy import ones

#==========================================================================
def properties_for_air(rc):
    """

    Returns physical properties of air for given resolution 'rc'

    For 60 deg from:
    http://www.engineeringtoolbox.com/air-properties-d_156.html
    -------------------------------------------------------------------------
    """

    # Create and fill matrice for all properties
    rho = ones(rc) * 1.067       # density              [kg/m^3]
    mu = ones(rc) * 20.17E-06    # viscosity            [Pa s]
    cp = ones(rc) * 1009         # thermal capacity     [J/kg/K]
    kappa = ones(rc) * 0.0285    # thermal conductivity [W/m/K]

    return rho, mu, cp, kappa  # end of function
