"""
Docstring.
"""
from numpy import reshape, prod, zeros
from scipy.sparse.linalg import bicgstab

from scrins.discretization.adj_n_bnds     import adj_n_bnds
from scrins.discretization.create_matrix  import create_matrix
from scrins.discretization.vol_balance    import vol_balance
from scrins.discretization.obst_zero_val  import obst_zero_val
from scrins.constants.solver import TOL

# ==========================================================================
def calc_p(p, uvwf, rho, dt, dxyz, obst):
    """
    Calculating the Pressure.
    """
    # -------------------------------------------------------------------------
    # Fetch the resolution
    rc = p.val.shape
    # Create linear system
    A_p, b_p = create_matrix(p, zeros(rc), dt/rho, dxyz, obst, 'n')

    # Compute the source for the pressure.  Important: don't send "obst"
    # as a parameter here, because you don't want to take it into
    # account at this stage.  After velocity corrections, you should.
    b_p = vol_balance(uvwf, dxyz, zeros(rc))

    print('Maximum volume error before correction: %12.5e' % abs(b_p).max())
    print('Volume imbalance before correction    : %12.5e' % b_p.sum())

    # Solve for pressure
    res = bicgstab(A_p, reshape(b_p, prod(rc)), tol=TOL)
    p.val[:] = reshape(res[0], rc)

    print("res[1] = ", res[1])

    # Anchor it to values around zero (the absolute value of pressure
    # correction can get really volatile.  Although it is in prinicple not
    # important for incompressible flows, it is ugly for post-processing.
    p.val[:] = p.val[:] - p.val.mean()

    # Set to zero in obstacle (it can get strange
    # values during the iterative solution procedure)
    if obst.any() != 0:
        p.val[:] = obst_zero_val(p.pos, p.val, obst)

    # Finally adjust the boundary values
    p = adj_n_bnds(p)

    return  # end of function
