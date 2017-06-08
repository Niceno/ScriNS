"""
Docstring.
"""
from numpy import empty
from scrins.operators.dif import dif

#==========================================================================
def cartesian_grid(xn, yn, zn):
    """
    Spawns a 3D Cartesian grid from three arrays with node coordinates.
    """
    # Compute cell resolutions
    nx = len(xn)-1
    ny = len(yn)-1
    nz = len(zn)-1

    # Create matrices for cell dimensions ...
    dx = empty((nx, ny, nz))
    dy = empty((nx, ny, nz))
    dz = empty((nx, ny, nz))

    # ... and fill them up!
    dx[:] = dif(xn).reshape(nx, 1, 1)
    dy[:] = dif(yn).reshape(1, ny, 1)
    dz[:] = dif(zn).reshape(1, 1, nz)

    # Compute resolutions for cell-centered and all collocated variables
    rc = nx, ny, nz
    ru = nx-1, ny, nz
    rv = nx, ny-1, nz
    rw = nx, ny, nz-1

    return nx, ny, nz, \
           dx, dy, dz, \
           rc, ru, rv, rw  # end of function
