"""
Docstring.
"""
from numpy import reshape
from scipy.sparse.linalg import bicgstab
from scipy import prod

from scrins.discretization.adj_n_bnds import adj_n_bnds
from scrins.discretization.advection import advection
from scrins.discretization.create_matrix import create_matrix
from scrins.constants.solver import TOL

# ==========================================================================
def calc_t(t, uvwf, rho_cap, kappa, dt, dxyz, obst):
# --------------------------------------------------------------------------
    """
    Docstring.
    """

    # Unpack tuple(s)
    dx, dy, dz = dxyz

    # Fetch the resolution
    rc = t.val.shape

    # Discretize the diffusive part
    A_t, b_t = create_matrix(t, rho_cap/dt, kappa, dxyz, obst, 'n')

    # The advective fluxes
    c_t = advection(rho_cap, t, uvwf, dxyz, dt, 'minmod')

    # Innertial term for enthalpy
    i_t = t.old * rho_cap * dx * dy* dz / dt

    # The entire source term
    f_t = b_t - c_t + i_t

    # Solve for temperature
    res0 = bicgstab(A_t, reshape(f_t, prod(rc)), tol=TOL)
    t.val[:] = reshape(res0[0], rc)

    adj_n_bnds(t)

    return  # end of function
