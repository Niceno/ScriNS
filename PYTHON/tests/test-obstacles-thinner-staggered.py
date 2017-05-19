#!/usr/bin/python

"""
                                                      o ... scalars
                         (n)                          - ... u velocities
                                                      | ... v velocities
      +-------+-------+-------+-------+-------+
      |       |       |       |       |       |
      |   o   -   o   -   o   -   o   -   o   | j=ny
      |       |       |       |       |       |
      +---|---+---|---+---|---+---|---+---|---+     j=nym
      |       |       |       |       |       |
      |   o   -   o   -   o   -   o   -   o   | ...
      |       |       |       |       |       |
 (w)  +---|---+---|---+---|---+---|---+---|---+    j=2        (e)
      |       |       |       |       |       |
      |   o   -   o   -   o   -   o   -   o   | j=2
      |       |       |       |       |       |
      +---|---+---|---+---|---+---|---+---|---+    j=1 (v-velocity)
      |       |       |       |       |       |
      |   o   -   o   -   o   -   o   -   o   | j=1   (scalar cell)
      |       |       |       |       |       |
      +-------+-------+-------+-------+-------+
 y       i=1     i=2     ...     ...     i=nx      (scalar cells)
^            i=1      i=2    ...    i=nxm      (u-velocity cells)
|
+---> x                  (s)
"""

from math import floor

from numpy import outer, zeros

from scrins.physical_models.properties_for_air import properties_for_air
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
from scrins.discretization.calc_uvw import calc_uvw
from scrins.discretization.corr_uvw import corr_uvw
from scrins.discretization.vol_balance import vol_balance
from scrins.display.print_time_step import print_time_step
from scrins.operators.avg import avg
from scrins.operators.par import par


#==========================================================================
#
# Define problem
#
#==========================================================================

# Node coordinates
xn = nodes(0, 1.25,  256)
yn = nodes(0, 0.125,  32)
zn = nodes(0, 0.125,  32)

# Cell coordinates
xc = avg(xn)
yc = avg(yn)
zc = avg(zn)

# Cell dimensions
nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

# Set physical properties
rho, mu, cap, kappa = properties_for_air(rc)
     
# Time-stepping parameters
dt  =    0.005  # time step
ndt = 2000      # number of time steps

# Create unknowns; names, positions and sizes
uf = create_unknown('face-u-vel',  X, ru, DIRICHLET)
vf = create_unknown('face-v-vel',  Y, rv, DIRICHLET)
wf = create_unknown('face-w-vel',  Z, rw, DIRICHLET)
p  = create_unknown('pressure',    C, rc, NEUMANN)

# Specify boundary conditions
uf.bnd[W].typ[:1,:,:] = DIRICHLET 
uf.bnd[W].val[:1,:,:]  = 0.1 * outer( par(1.0, yn), par(1.0, zn) ) 

uf.bnd[E].typ[:1,:,:] = OUTLET 

for j in (B,T):
  uf.bnd[j].typ[:] = NEUMANN     
  vf.bnd[j].typ[:] = NEUMANN     
  wf.bnd[j].typ[:] = NEUMANN     

adj_n_bnds(p)

# Create obstacles
obst = zeros(rc)

class key:
  ip = -1
  im = -1
  jp = -1
  jm = -1
  kp = -1
  km = -1
  
block = (key(), key(), key(), key());

th = 5;
block[0].im =  3*nx/16            # i minus
block[0].ip =  block[0].im + th   # i plus
block[0].jm =  0                  # j minus
block[0].jp =  3*ny/4             # j plus
block[0].km =  0                  # k minus
block[0].kp =  3*ny/4             # k plus

block[1].im =  5*nx/16            # i minus
block[1].ip =  block[1].im + th   # i plus
block[1].jm =  ny/4               # j minus
block[1].jp =  ny                 # j plus
block[1].km =  ny/4               # k minus
block[1].kp =  ny                 # k plus
  
block[2].im =  7*nx/16            # i minus
block[2].ip =  block[2].im + th   # i plus
block[2].jm =  0                  # j minus
block[2].jp =  3*ny/4             # j plus
block[2].km =  0                  # k minus
block[2].kp =  3*ny/4             # k plus

block[3].im =  9*nx/16            # i minus
block[3].ip =  block[3].im + th   # i plus
block[3].jm =  ny/4               # j minus
block[3].jp =  ny                 # j plus
block[3].km =  ny/4               # k minus
block[3].kp =  ny                 # k plus

for o in range(0, 4):
  for i in range(floor(block[o].im), floor(block[o].ip)):  
    for j in range(floor(block[o].jm), floor(block[o].jp)):
      for k in range(floor(block[o].km), floor(block[o].kp)):
        obst[i,j,k] = 1

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
  ef = zeros(ru), zeros(rv), zeros(rw)
    
  calc_uvw((uf,vf,wf), (uf,vf,wf), rho, mu,  \
           zeros(rc), ef, dt, (dx,dy,dz), obst)
  
  #----------
  # Pressure
  #----------
  calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obst)
  
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

  if ts % 20 == 0:
    plot_isolines(p.val, (uf,vf,wf), (xn,yn,zn), Y)
    plot_isolines(p.val, (uf,vf,wf), (xn,yn,zn), Z)    
