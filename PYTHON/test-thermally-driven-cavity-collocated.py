#!/usr/bin/python

# Standard Python modules
from standard import *

# ScriNS modules
from Constants.all      import *
from Operators.all      import *
from Display.all        import *
from Discretization.all import *
from PhysicalModels.all import *

#==========================================================================
#
# Define problem
#
#==========================================================================

xn = nodes(0, 1,  64, 1.0/256, 1.0/256)
yn = nodes(0, 1,  64, 1.0/256, 1.0/256)
zn = nodes(0, 0.1, 5)

# Cell dimensions
nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

# Set physical properties
grashof = 1.4105E+06
prandtl = 0.7058
rho   = zeros(rc)
mu    = zeros(rc)
kappa = zeros(rc)
cap   = zeros(rc)
rho  [:,:,:] = 1.0 
mu   [:,:,:] = 1.0 / sqrt(grashof)
kappa[:,:,:] = 1.0 / (prandtl * sqrt(grashof)) 
cap  [:,:,:] = 1.0 
     
# Time-stepping parameters
dt     =    0.02       # time step
ndt    = 1000          # number of time steps

# Create unknowns; names, positions and sizes
uc = create_unknown('cell-u-vel',  C, rc, DIRICHLET)
vc = create_unknown('cell-v-vel',  C, rc, DIRICHLET)
wc = create_unknown('cell-w-vel',  C, rc, DIRICHLET)
uf = create_unknown('face-u-vel',  X, ru, DIRICHLET)
vf = create_unknown('face-v-vel',  Y, rv, DIRICHLET)
wf = create_unknown('face-w-vel',  Z, rw, DIRICHLET)
t  = create_unknown('temperature', C, rc, NEUMANN)
p  = create_unknown('pressure',    C, rc, NEUMANN)
p_tot = zeros(rc)

# This is a new test
t.bnd[W].typ[:] = DIRICHLET
t.bnd[W].val[:] = -0.5

t.bnd[E].typ[:] = DIRICHLET
t.bnd[E].val[:] = +0.5
 
for j in (B,T):
  uc.bnd[j].typ[:] = NEUMANN     
  vc.bnd[j].typ[:] = NEUMANN     
  wc.bnd[j].typ[:] = NEUMANN     

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
  t.old[:]  = t.val[:]
  uc.old[:] = uc.val[:]
  vc.old[:] = vc.val[:]
  wc.old[:] = wc.val[:]
  
  #------------------------
  # Temperature (enthalpy)
  #------------------------
  calc_t(t, (uf,vf,wf), (rho*cap), kappa, dt, (dx,dy,dz), obst)

  #-----------------------
  # Momentum conservation
  #-----------------------
  ef = zeros(rc), t.val, zeros(rc)
  
  calc_uvw((uc,vc,wc), (uf,vf,wf), rho, mu, 
           p_tot, ef, dt, (dx,dy,dz), obst)
  
  #----------
  # Pressure
  #----------
  calc_p(p, (uf,vf,wf), rho, dt, (dx,dy,dz), obst)
  
  p_tot = p_tot + p.val
  
  #---------------------
  # Velocity correction
  #---------------------
  corr_uvw((uc,vc,wc), p, rho, dt, (dx,dy,dz), obst)
  corr_uvw((uf,vf,wf), p, rho, dt, (dx,dy,dz), obst)
 
  # Compute volume balance for checking 
  err = vol_balance((uf,vf,wf), (dx,dy,dz), obst)
  print('Maximum volume error after correction: %12.5e' % abs(err).max())

  # Check the CFL number too 
  cfl = cfl_max((uc,vc,wc), dt, (dx,dy,dz))
  print('Maximum CFL number: %12.5e' % cfl)

#==========================================================================
#
# Visualisation
#
#==========================================================================

  if ts % 20 == 0:
    plot_isolines(t.val, (uc,vc,wc), (xn,yn,zn), Z)
  