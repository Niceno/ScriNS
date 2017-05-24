#!/usr/bin/python

"""
Demonstrates the variation of projection algorythm which computes total
pressure, as a sum of all pressure corrections.

The total pressure being built up in this way counter-balances the
gravity term in momentum equations.

It seems that such an approach is important for bouyancy dominated flows.

Gravity term is under-relaxed here, but it works even without it.
"""

from numpy import outer, zeros

from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN, OUTLET
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.constants.gravitational_constant import G
from scrins.display.plot_isolines import plot_isolines
from scrins.discretization.adj_n_bnds import adj_n_bnds
from scrins.discretization.cartesian_grid import cartesian_grid
from scrins.discretization.nodes import nodes
from scrins.discretization.create_unknown import create_unknown
from scrins.discretization.cfl_max import cfl_max
from scrins.discretization.calc_p import calc_p
from scrins.discretization.calc_uvw import calc_uvw
from scrins.discretization.corr_uvw import corr_uvw
from scrins.discretization.vol_balance import vol_balance
from scrins.display.print_time_step import print_time_step
from scrins.operators.avg import avg
from scrins.operators.par import par

def main(show_plot=True):

    #==========================================================================
    #
    # Define problem
    #
    #==========================================================================

    # Node coordinates
    xn = nodes(0, 10,   80)
    yn = nodes(0,  1,   20)
    zn = nodes(0,  0.5,  5)

    # Cell coordinates
    xc = avg(xn)
    yc = avg(yn)
    zc = avg(zn)

    # Cell dimensions
    nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

    # Set physical properties
    rho   = zeros(rc)
    mu    = zeros(rc)
    kappa = zeros(rc)
    cap   = zeros(rc)
    rho   [:] =    1.15  # density              [kg/m^3]
    mu    [:] =    0.1   # viscosity            [Pa s]

    # Time-stepping parameters
    dt  =   0.1  # time step
    ndt = 200    # number of time steps

    # Create unknowns; names, positions and sizes
    uf = create_unknown('face-u-vel',  X, ru, DIRICHLET)
    vf = create_unknown('face-v-vel',  Y, rv, DIRICHLET)
    wf = create_unknown('face-w-vel',  Z, rw, DIRICHLET)
    p  = create_unknown('pressure',    C, rc, NEUMANN)
    p_tot = zeros(rc)

    # Specify boundary conditions
    uf.bnd[W].typ[:1,:,:] = DIRICHLET
    for k in range(0,nz):
      uf.bnd[W].val[:1,:,k]  = par(0.1, yn)

    uf.bnd[E].typ[:1,:,:] = OUTLET

    for j in (B,T):
      uf.bnd[j].typ[:] = NEUMANN
      vf.bnd[j].typ[:] = NEUMANN
      wf.bnd[j].typ[:] = NEUMANN

    adj_n_bnds(p)

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
    for ts in range(1,ndt+1):

      print_time_step(ts)

      #------------------
      # Store old values
      #------------------
      uf.old[:] = uf.val[:]
      vf.old[:] = vf.val[:]
      wf.old[:] = wf.val[:]

      #-----------------------
      # Momentum conservation
      #-----------------------
      g_v = -G * avg(Y, rho) * min(ts/100,1)

      ef = zeros(ru), g_v, zeros(rw)

      calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu,  \
               p_tot, ef, dt, (dx,dy,dz), obst)

      #----------
      # Pressure
      #----------
      calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obst)

      p_tot = p_tot + p.val

      #---------------------
      # Velocity correction
      #---------------------
      corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), obst)

      # Compute volume balance for checking
      err = vol_balance((uf,vf,wf), (dx,dy,dz), obst)
      print('Maximum volume error after correction: %12.5e' % abs(err).max())

      # Check the CFL number too
      cfl = cfl_max((uf,vf,wf), dt, (dx,dy,dz))
      print('Maximum CFL number: %12.5e' % cfl)

    #==========================================================================
    #
    # Visualisation
    #
    #==========================================================================

      if show_plot:
          if ts % 20 == 0:
            plot_isolines(p_tot, (uf,vf,wf), (xn,yn,zn), Z)
            plot_isolines(p.val, (uf,vf,wf), (xn,yn,zn), Z)

if __name__ == '__main__':
    main()
