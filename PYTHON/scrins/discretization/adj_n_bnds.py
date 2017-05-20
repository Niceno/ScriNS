"""
Adjust Neumann boundary conditions.
"""

from scipy import logical_not as lnot
from scrins.constants.boundary_conditions import NEUMANN
from scrins.constants.compass import W, E, S, N, B, T


def adj_n_bnds(phi):
    """Copies last domain cell values to Neumann boundary condition values.
    """
    # These arrays will hold values true (0) in cells with boundary ...
    # ... condition of Neumann type, and false (0) otherwise

    if_w_n = phi.bnd[W].typ[0, :, :] == NEUMANN  # 1 if west is Neumann
    if_e_n = phi.bnd[E].typ[0, :, :] == NEUMANN  # 1 if east is Neumann
    if_s_n = phi.bnd[S].typ[:, 0, :] == NEUMANN  # 1 if south is Neumann
    if_n_n = phi.bnd[N].typ[:, 0, :] == NEUMANN  # 1 if north is Neumann
    if_b_n = phi.bnd[B].typ[:, :, 0] == NEUMANN  # 1 if bottom is Neumann
    if_t_n = phi.bnd[T].typ[:, :, 0] == NEUMANN  # 1 if top is Neumann

    # In what follows, a linear combination of true (0] and false (0)
    # will copy the values of variable phi to the boundaries.
    phi.bnd[W].val[0, :, :] = (phi.bnd[W].val[0, :, :] * lnot(if_w_n)
                               + phi.val[0, :, :] * if_w_n)
    phi.bnd[E].val[0, :, :] = (phi.bnd[E].val[0, :, :] * lnot(if_e_n)
                               + phi.val[-1, :, :] * if_e_n)

    phi.bnd[S].val[:, 0, :] = (phi.bnd[S].val[:, 0, :] * lnot(if_s_n)
                               + phi.val[:, 0, :] * if_s_n)
    phi.bnd[N].val[:, 0, :] = (phi.bnd[N].val[:, 0, :] * lnot(if_n_n)
                               + phi.val[:, -1, :] * if_n_n)

    phi.bnd[B].val[:, :, 0] = (phi.bnd[B].val[:, :, 0] * lnot(if_b_n)
                               + phi.val[:, :, 0] * if_b_n)
    phi.bnd[T].val[:, :, 0] = (phi.bnd[T].val[:, :, 0] * lnot(if_t_n)
                               + phi.val[:, :, -1] * if_t_n)
