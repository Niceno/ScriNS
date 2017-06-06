"""
Docstring.
"""
from collections import namedtuple

from numpy import array, reshape, zeros
from scipy.sparse import spdiags

from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.obst_mod_matrix import obst_mod_matrix
from scrins.operators.cat import cat
from scrins.operators.avg import avg


# ==========================================================================
def create_matrix(phi, inn, mu, dxyz, obst, obc):

    """
    pos   - position of variable (C - central,
                                  X - staggered in x direction,
                                  Y - staggered in y direction,
                                  Z - staggered in z direction)
    inn         - innertial term
    mu         - viscous coefficient
    dx, dy, dz - cell size in x, y and z directions
    obc        - obstacles's boundary condition, ('n' - Neumann,
                                                  'd' - Dirichlet)
    -------------------------------------------------------------------------
    """
    # Unpack tuples
    dx, dy, dz = dxyz

    res = phi.val.shape

    # -------------------------------
    # Create right hand side vector
    # -------------------------------
    b = zeros(res)

    # ------------------------------------
    # Create default matrix coefficients
    # ------------------------------------
    coefficients = namedtuple('matrix_diagonal', 'W E S N B T P')
    c = coefficients(zeros(res), zeros(res), zeros(res), \
                     zeros(res), zeros(res), zeros(res), \
                     zeros(res))

    d = phi.pos

    # Handle central coefficient due to innertia
    c.P[:] = avg(d, inn) * avg(d, dx*dy*dz)

    # Pre-compute geometrical quantities
    sx = dy * dz
    sy = dx * dz
    sz = dx * dy

    if d != X:
        c.W[:] = cat(X, (                                  \
        avg(d, mu[:1, :, :])* avg(d, sx[:1, :, :]) / avg(d, (dx[:1, :, :])/2.0), \
        avg(d, avg(X, mu))* avg(d, avg(X, sx)) / avg(d, avg(X, dx))))

        c.E[:] = cat(X, (                                                    \
        avg(d, avg(X, mu))  * avg(d, avg(X, sx))  / avg(d, avg(X, dx)),       \
        avg(d, mu[-1:, :, :])* avg(d, sx[-1:, :, :])/ avg(d, (dx[-1:, :, :])/2.0)))

    if d != Y:
        c.S[:] = cat(Y, (                                                     \
        avg(d, mu[:, :1, :]) * avg(d, sy[:, :1, :]) / avg(d, (dy[:, :1, :])/2.0),  \
        avg(d, avg(Y, mu))  * avg(d, avg(Y, sy))  / avg(d, avg(Y, dy))))

        c.N[:] = cat(Y, (                                                     \
        avg(d, avg(Y, mu))  * avg(d, avg(Y, sy))  / avg(d, avg(Y, dy)),       \
        avg(d, mu[:, -1:, :]) * avg(d, sy[:, -1:, :]) / avg(d, (dy[:, -1:, :])/2.0)))

    if d != Z:
        c.B[:] = cat(Z, (                                                     \
        avg(d, mu[:, :, :1]) * avg(d, sz[:, :, :1]) / avg(d, (dz[:, :, :1])/2.0),  \
        avg(d, avg(Z, mu))  * avg(d, avg(Z, sz))  / avg(d, avg(Z, dz))))

        c.T[:] = cat(Z, (                                                     \
        avg(d, avg(Z, mu))  * avg(d, avg(Z, sz))  / avg(d, avg(Z, dz)),       \
        avg(d, mu[:, :, -1:]) * avg(d, sz[:, :, -1:]) / avg(d, (dz[:, :, -1:])/2.0)))

     # ---------------------------------
     # Correct for staggered variables
     # ---------------------------------
    if d == X:
        c.W[:] = mu[0:-1, :, :] * sx[0:-1, :, :] / dx[0:-1, :, :]
        c.E[:] = mu[1:, :, :] * sx[1:, :, :] / dx[1:, :, :]
    elif d == Y:
        c.S[:] = mu[:, 0:-1, :] * sy[:, 0:-1, :] / dy[:, 0:-1, :]
        c.N[:] = mu[:, 1:, :] * sy[:, 1:, :] / dy[:, 1:, :]
    elif d == Z:
        c.B[:] = mu[:, :, 0:-1] * sz[:, :, 0:-1] / dz[:, :, 0:-1]
        c.T[:] = mu[:, :, 1:  ] * sz[:, :, 1:  ] / dz[:, :, 1:  ]

    # -----------------------------------------------------------------------
    # Zero them (correct them) for vanishing derivative boundary condition.
    # -----------------------------------------------------------------------

    # The values defined here will be false (numerical value 0)
    # wherever there is or Neumann boundary condition.
    c.W[:1, :, :] *= (phi.bnd[W].typ[:] == DIRICHLET)
    c.E[-1:, :, :] *= (phi.bnd[E].typ[:] == DIRICHLET)
    c.S[:, :1, :] *= (phi.bnd[S].typ[:] == DIRICHLET)
    c.N[:, -1:, :] *= (phi.bnd[N].typ[:] == DIRICHLET)
    c.B[:, :, :1] *= (phi.bnd[B].typ[:] == DIRICHLET)
    c.T[:, :, -1:] *= (phi.bnd[T].typ[:] == DIRICHLET)

     # --------------------------------------------
     # Fill the source terms with boundary values
     # --------------------------------------------
    b[:1, :, :] += c.W[ :1, :, :] * phi.bnd[W].val[:1, :, :]
    b[-1:, :, :] += c.E[-1:, :, :] * phi.bnd[E].val[:1, :, :]
    b[:, :1, :] += c.S[  :, :1, :] * phi.bnd[S].val[:, :1, :]
    b[:, -1:, :] += c.N[  :, -1:, :] * phi.bnd[N].val[:, :1, :]
    b[:, :, :1] += c.B[  :, :, :1] * phi.bnd[B].val[:, :, :1]
    b[:, :, -1:] += c.T[  :, :, -1:] * phi.bnd[T].val[:, :, :1]

     # ---------------------------------------
     # Correct system matrices for obstacles
     # ---------------------------------------
    if obst.any() != 0:
        c = obst_mod_matrix(phi, c, obst, obc)

     # -----------------------------------------------
     # Add all neighbours to the central matrix,
     # and zero the coefficients towards boundaries
     # -----------------------------------------------
    c.P[:] += c.W[:] + c.E[:] + c.S[:] + c.N[:] + c.B[:] + c.T[:]

    c.W[:1, :, :] = 0.0
    c.E[-1:, :, :] = 0.0
    c.S[:, :1, :] = 0.0
    c.N[:, -1:, :] = 0.0
    c.B[:, :, :1] = 0.0
    c.T[:, :, -1:] = 0.0

    # ----------------------
    # Create sparse matrix
    # ----------------------
    nx, ny, nz = res
    n = nx * ny * nz

    data = array([reshape(c.P, n),                     \
                  reshape(-c.W, n), reshape(-c.E, n),   \
                  reshape(-c.S, n), reshape(-c.N, n),   \
                  reshape(-c.B, n), reshape(-c.T, n)])

    diag = array([0, +ny*nz, -ny*nz, +nz, -nz, +1, -1])

    A = spdiags(data, diag, n, n)

    return A, b  # end of function
