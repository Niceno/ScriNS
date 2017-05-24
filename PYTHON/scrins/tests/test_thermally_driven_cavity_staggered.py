#!/usr/bin/python

"""
This program solves thermally driven cavity at Ra = 1.0e6,
in dimensional and non-dimensional forms.

Equations in dimensional form:

D(rho u)/Dt = nabla(mu (nabla u)^T) - nabla p + g
D(rho cp T)/Dt = nabla(lambda (nabla T)^T)

Equations in non-dimensional form, for natural convection problems

DU/Dt = nabla(1/sqrt(Gr) (nabla U)^T) - nabla P + theta
D theta/Dt = nabla(1/(Pr*sqrt(Gr)) (nabla theta)^T)

--------------------------------------------------------------------------
For thermally driven cavity, with properties of air at 60 deg:

nu   =  1.89035E-05;
beta =  0.003;
dT   = 17.126;
L    =  0.1;
Pr   = 0.709;

characteristic non-dimensional numbers are:
Gr = 1.4105E+06
Ra = 1.0000E+06
"""

from numpy import sqrt, zeros

from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN, OUTLET
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.display.plot_isolines import plot_isolines
from scrins.discretization.adj_n_bnds import adj_n_bnds
from scrins.discretization.cartesian_grid import cartesian_grid
from scrins.discretization.nodes import nodes
from scrins.discretization.create_unknown import create_unknown
from scrins.discretization.cfl_max import cfl_max
from scrins.discretization.calc_p import calc_p
from scrins.discretization.calc_t import calc_t
from scrins.discretization.calc_uvw import calc_uvw
from scrins.discretization.corr_uvw import corr_uvw
from scrins.discretization.vol_balance import vol_balance
from scrins.display.print_time_step import print_time_step
from scrins.operators.avg import avg
from scrins.operators.par import par


def main(show_plot=True):


    # ==========================================================================
    #
    # Define problem
    #
    # ==========================================================================

    xn = nodes(0, 1, 64, 1.0 / 256, 1.0 / 256)
    yn = nodes(0, 1, 64, 1.0 / 256, 1.0 / 256)
    zn = nodes(0, 0.1, 5)

    # Cell dimensions
    nx, ny, nz, dx, dy, dz, rc, ru, rv, rw = cartesian_grid(xn, yn, zn)

    # Set physical properties
    grashof = 1.4105E+06
    prandtl = 0.7058
    rho = zeros(rc)
    mu = zeros(rc)
    kappa = zeros(rc)
    cap = zeros(rc)
    rho[:, :, :] = 1.0
    mu[:, :, :] = 1.0 / sqrt(grashof)
    kappa[:, :, :] = 1.0 / (prandtl * sqrt(grashof))
    cap[:, :, :] = 1.0

    # Time-stepping parameters
    dt = 0.02  # time step
    ndt = 1000  # number of time steps

    # Create unknowns; names, positions and sizes
    uf = create_unknown('face-u-vel', X, ru, DIRICHLET)
    vf = create_unknown('face-v-vel', Y, rv, DIRICHLET)
    wf = create_unknown('face-w-vel', Z, rw, DIRICHLET)
    t = create_unknown('temperature', C, rc, NEUMANN)
    p = create_unknown('pressure', C, rc, NEUMANN)
    p_tot = zeros(rc)

    # This is a new test
    t.bnd[W].typ[:] = DIRICHLET
    t.bnd[W].val[:] = -0.5

    t.bnd[E].typ[:] = DIRICHLET
    t.bnd[E].val[:] = +0.5

    for j in (B, T):
        uf.bnd[j].typ[:] = NEUMANN
        vf.bnd[j].typ[:] = NEUMANN
        wf.bnd[j].typ[:] = NEUMANN

    obst = zeros(rc)

    # ==========================================================================
    #
    # Solution algorithm
    #
    # ==========================================================================

    # -----------
    #
    # Time loop
    #
    # -----------
    for ts in range(1, ndt + 1):

        print_time_step(ts)

        # ------------------
        # Store old values
        # ------------------
        t.old[:] = t.val[:]
        uf.old[:] = uf.val[:]
        vf.old[:] = vf.val[:]
        wf.old[:] = wf.val[:]

        # ------------------------
        # Temperature (enthalpy)
        # ------------------------
        calc_t(t, (uf, vf, wf), (rho * cap), kappa, dt, (dx, dy, dz), obst)

        # -----------------------
        # Momentum conservation
        # -----------------------
        ef = zeros(ru), avg(Y, t.val), zeros(rw)

        calc_uvw((uf, vf, wf), (uf, vf, wf), rho, mu, \
                 p_tot, ef, dt, (dx, dy, dz), obst)

        # ----------
        # Pressure
        # ----------
        calc_p(p, (uf, vf, wf), rho, dt, (dx, dy, dz), obst)

        p_tot = p_tot + p.val

        # ---------------------
        # Velocity correction
        # ---------------------
        corr_uvw((uf, vf, wf), p, rho, dt, (dx, dy, dz), obst)

        # Compute volume balance for checking
        err = vol_balance((uf, vf, wf), (dx, dy, dz), obst)
        print('Maximum volume error after correction: %12.5e' % abs(err).max())

        # Check the CFL number too
        cfl = cfl_max((uf, vf, wf), dt, (dx, dy, dz))
        print('Maximum CFL number: %12.5e' % cfl)

        # ==========================================================================
        #
        # Visualisation
        #
        # ==========================================================================
        if show_plot:
            if ts % 20 == 0:
                plot_isolines(t.val, (uf, vf, wf), (xn, yn, zn), Z)

if __name__ == '__main__':
    main()
