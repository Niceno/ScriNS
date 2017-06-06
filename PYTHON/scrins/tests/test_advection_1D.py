#!/usr/bin/python

"""
Program to test implementation of advection schemes in the code, using
one-dimensional transport of a step function in X, Y or Z
direction, either in positive or negative sense.

The coordinate direction is specified with the local variable "TEST",
which can be either X, Y or Z.

Sense is specified with the variable "FLOW", which can assume values 'p' for
positive, and 'n' for negative sense.
"""


import random

from matplotlib import pyplot as plt
from numpy import zeros

from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.adj_n_bnds import adj_n_bnds
from scrins.discretization.calc_t import calc_t
from scrins.discretization.cartesian_grid import cartesian_grid
from scrins.discretization.nodes import nodes
from scrins.discretization.create_unknown import create_unknown
from scrins.display.print_time_step import print_time_step
from scrins.operators.avg import avg


def main(show_plot=True):
    """
    1D test
    :param show_plot: show plot or not
    :return: None
    """

    #==========================================================================
    #
    # Define problem
    #
    #==========================================================================


    SENS = random.choice([-1, +1])
    TEST = random.choice([X, Y, Z])
    MIN = 10
    MAX = 20

    uvw_bulk = SENS, 0.0, 0.0
    xn = nodes(0, 1, 200)
    yn = nodes(0, 0.1, 4)
    zn = nodes(0, 0.1, 4)
    if TEST == Y:
        uvw_bulk = 0.0, SENS, 0.0
        xn = nodes(0, 0.1, 4)
        yn = nodes(0, 1, 200)
        zn = nodes(0, 0.1, 4)
    elif TEST == Z:
        uvw_bulk = 0.0, 0.0, SENS
        xn = nodes(0, 0.1, 4)
        yn = nodes(0, 0.1, 4)
        zn = nodes(0, 1, 200)

    # Set domain dimensions and grid resolution

    # Cell dimensions
    nx, ny, nz,\
    dx, dy, dz, \
    rc, ru, rv, rw = cartesian_grid(xn, yn, zn)

    # Set physical properties
    rho = zeros(rc)
    mu = zeros(rc)
    kappa = zeros(rc)
    cap = zeros(rc)
    rho  [:, :, :] = 1.0
    mu   [:, :, :] = 0.0
    kappa[:, :, :] = 0.0
    cap  [:, :, :] = 1.0

    # Time-stepping parameters
    dt = 0.002   # time step
    ndt = 350    # number of time steps

    # Create unknowns; names, positions and sizes
    uf = create_unknown('face-u-vel', X, ru, DIRICHLET)
    vf = create_unknown('face-v-vel', Y, rv, DIRICHLET)
    wf = create_unknown('face-w-vel', Z, rw, DIRICHLET)
    t = create_unknown('temperature', C, rc, NEUMANN)

    # Make a tuple with all velocity components
    uvwf = uf, vf, wf

    # Specify inital and boundary conditions
    t.val[:] = MIN
    if TEST == X:
        t.val[nx//2-25:nx//2+25, :, :] = MAX
    elif TEST == Y:
        t.val[:, ny//2-25:ny//2+25, :] = MAX
    elif TEST == Z:
        t.val[:, :, nz//2-25:nz//2+25] = MAX

    adj_n_bnds(t)

    # This is from the advection-1D case
    for i in (X, Y, Z):
        uvwf[i].val[:, :, :] = uvw_bulk[i]
        for j in (W, E, S, N, B, T):
            uvwf[i].bnd[j].val[:, :, :] = uvw_bulk[i]

    obst = zeros(rc)

    #==========================================================================
    #
    # Solution algorithm
    #
    #==========================================================================

    #-----------
    #
    # Time loop
    #
    #-----------
    for ts in range(1, 40):

        print_time_step(ts)

        #------------------
        # Store old values
        #------------------
        t.old[:] = t.val[:]

        #------------------------
        # Temperature (enthalpy)
        #------------------------
        calc_t(t, uvwf, (rho*cap), kappa, dt, (dx, dy, dz), obst)

        #==========================================================================
        #
        # Visualisation
        #
        #==========================================================================
    st = 'ro'
    xc = avg(xn)
    vc = t.val[:, ny//2, nz//2]
    if TEST == Y:
        st = 'go'
        xc = avg(yn)
        vc = t.val[nx//2, :, nz//2]
    elif TEST == Z:
        st = 'bo'
        xc = avg(zn)
        vc = t.val[nx//2, ny//2, :]

    if show_plot:
        plt.plot(xc, vc, st)
        plt.show()

if __name__ == '__main__':
    main()
