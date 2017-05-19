#==========================================================================
# This scripts tests inlet and outlet conditons at various places in
# computational domain.
#
# Computational domain is a simple box, and Reynolds number is rather 
# small to avoid instabilities due to vortices getting out from the 
# outlet. The script selects the case it will run randomply, from the 
# 16 predefined cases.
#--------------------------------------------------------------------------

#!/usr/bin/python

# Standard Python modules
from standard import *

# ScriNS modules
from Constants.all      import *
from Operators.all      import *
from Display.all        import *
from Discretization.all import *
from PhysicalModels.all import properties_for_air

#==========================================================================
#
# Define problem
#
#==========================================================================

planes = ('XY', 'XZ', 'YZ')
tests  = array([11,12,13,14, 21,22,23,24, 31,32,33,34, 41,42,43,44])
TEST   = tests[ floor(random()*16) ]
PLANE  = planes[ floor(random()*3) ]

# Node coordinates
xn = nodes(0, 1,    160)
yn = nodes(0, 1,    160)
zn = nodes(0, 0.025,  4)

# Cell coordinates
xc = avg(xn)
yc = avg(yn)
zc = avg(zn)

# Cell dimensions
nx,ny,nz, dx,dy,dz, rc,ru,rv,rw = cartesian_grid(xn,yn,zn)

# Set physical properties
rho, mu, cap, kappa = properties_for_air(rc)
     
# Time-stepping parameters
dt  =   0.15  # time step
ndt = 200     # number of time steps

# Create unknowns names, positions and sizes
uf = create_unknown('face-u-vel',  X, ru, DIRICHLET)
vf = create_unknown('face-v-vel',  Y, rv, DIRICHLET)
wf = create_unknown('face-w-vel',  Z, rw, DIRICHLET)
p  = create_unknown('pressure',    C, rc, NEUMANN)

print(TEST)

# Specify boundary conditions
if TEST == 11:
  for k in range(0,nz):
    uf.bnd[W].val[0,ny/4:3*ny/4,k] = +par(0.01, yn[ny/4:3*ny/4+1])    
    uf.bnd[E].typ[0,ny/4:3*ny/4,k] =  OUTLET  
    
elif TEST == 12:  # vertical mirror from 11
  for k in range(0,nz):
    uf.bnd[E].val[0,ny/4:3*ny/4,k] = -par(0.01, yn[ny/4:3*ny/4+1])    
    uf.bnd[W].typ[0,ny/4:3*ny/4,k] =  OUTLET   
    
elif TEST == 13:  # rotate 11
  for k in range(0,nz):
    vf.bnd[S].val[nx/4:3*nx/4,0,k] = +par(0.01, xn[nx/4:3*nx/4+1])    
    vf.bnd[N].typ[nx/4:3*nx/4,0,k] = OUTLET   
    
elif TEST == 14:  # horizontal mirror 13
  for k in range(0,nz):
    vf.bnd[N].val[nx/4:3*nx/4,0,k] = -par(0.01, xn[nx/4:3*nx/4+1])
    vf.bnd[S].typ[nx/4:3*nx/4,0,k] =  OUTLET   
    
elif TEST == 21:  # 2 exits
  for k in range(0,nz):
    uf.bnd[W].val[0,  ny/4:3*ny/4,k] = +par(0.01, yn[ny/4:3*ny/4+1])    
    uf.bnd[E].typ[0,       0:ny/4,k] = OUTLET   
    uf.bnd[E].typ[0,3*ny/4:ny,    k] = OUTLET   
    
elif TEST == 22:  # vertical mirror 21
  for k in range(0,nz):
    uf.bnd[E].val[0,  ny/4:3*ny/4,k] = -par(0.01, yn[ny/4:3*ny/4+1])    
    uf.bnd[W].typ[0,       0:ny/4,k] = OUTLET   
    uf.bnd[W].typ[0,3*ny/4:ny,    k] = OUTLET   
    
elif TEST == 23:  # rotated 21
  for k in range(0,nz):
    vf.bnd[S].val[  nx/4:3*nx/4,0,k] = +par(0.01, xn[nx/4:3*nx/4+1])
    vf.bnd[N].typ[        :nx/4,0,k] = OUTLET   
    vf.bnd[N].typ[3*nx/4:nx,    0,k] = OUTLET   
    
elif TEST == 24:  # horizontal mirror of 23
  for k in range(0,nz):
    vf.bnd[N].val[  nx/4:3*nx/4,0,k] = -par(0.01, xn[nx/4:3*nx/4+1])
    vf.bnd[S].typ[        :nx/4,0,k] =  OUTLET   
    vf.bnd[S].typ[3*nx/4:nx,    0,k] =  OUTLET   
    
elif TEST == 31:  # inlet and outlet at the same face
  for k in range(0,nz):
    uf.bnd[W].val[0,3*ny/4:ny,  k] = +par(0.01, yn[3*ny/4:ny+1])    
    uf.bnd[W].typ[0,      :ny/4,k] =  OUTLET   
    
elif TEST == 32:  # vertical mirror of 31
  for k in range(0,nz):
    uf.bnd[E].val[0,3*ny/4:ny,  k] = -par(0.01, yn[3*ny/4:ny+1])
    uf.bnd[E].typ[0,     0:ny/4,k] =  OUTLET   
    
elif TEST == 33:  # rotated 31
  for k in range(0,nz):
    vf.bnd[S].val[3*nx/4:nx,  0,k] = +par(0.01, xn[3*nx/4:nx+1])    
    vf.bnd[S].typ[      :nx/4,0,k] =  OUTLET   
    
elif TEST == 34:  # horizontal mirror of 33
  for k in range(0,nz):
    vf.bnd[N].val[3*nx/4:nx,  0,k] = -par(0.01, xn[3*nx/4:nx+1])    
    vf.bnd[N].typ[      :nx/4,0,k] =  OUTLET   
    
elif TEST == 41:  # inlet and outlet at the same face, one more outlet
  for k in range(0,nz):
    uf.bnd[W].val[0,3*ny/4:ny,  k] = +par(0.01, yn[3*ny/4:ny+1])
    uf.bnd[W].typ[0,      :ny/8,k] = OUTLET   
    uf.bnd[E].typ[0,      :ny/8,k] = OUTLET   
    
elif TEST == 42:  # vertical mirror of 41
  for k in range(0,nz):
    uf.bnd[E].val[0,3*ny/4:ny,  k] = -par(0.01, yn[3*ny/4:ny+1])
    uf.bnd[E].typ[0,      :ny/8,k] = OUTLET   
    uf.bnd[W].typ[0,      :ny/8,k] = OUTLET   
    
elif TEST == 43:  # rotated 41
  for k in range(0,nz):
    vf.bnd[S].val[3*nx/4:nx,  0,k] = +par(0.01, xn[3*nx/4:nx+1])
    vf.bnd[S].typ[     0:nx/8,0,k] =  OUTLET   
    vf.bnd[N].typ[     0:nx/8,0,k] =  OUTLET   
    
elif TEST == 44:  # horizontal mirror of 43
  for k in range(0,nz):
    vf.bnd[N].val[3*ny/4:nx,  0,k] = -par(0.01, xn[3*nx/4:nx+1])
    vf.bnd[N].typ[      :nx/8,0,k] =  OUTLET   
    vf.bnd[S].typ[      :nx/8,0,k] =  OUTLET 

for j in (B,T):
  uf.bnd[j].typ[:] = NEUMANN     
  vf.bnd[j].typ[:] = NEUMANN     
  wf.bnd[j].typ[:] = NEUMANN     

# Create a cylindrical obstacle in the middle just for kicks
obst = zeros(rc)
for k in range(0,nz):
  for j in range(0,ny):
    for i in range(0,nx):
      dist = sqrt( (j-ny/2+1)**2 + (i-nx/2+1)**2 )
      if dist < ny/4:
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

  if ts % 10 == 0:
    plot_isolines(p.val, (uf,vf,wf), (xn,yn,zn), Z)
