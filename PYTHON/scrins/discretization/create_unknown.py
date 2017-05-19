from collections import namedtuple

from numpy import ndarray, zeros

from scrins.constants.boundary_conditions import DIRICHLET, NEUMANN, OUTLET
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C

#==========================================================================
def create_unknown(name, pos, res, def_bc):
#--------------------------------------------------------------------------
# This function creates a new unkown; helping to shorten the main program.
#--------------------------------------------------------------------------
# Input parameters are:
#
# name   - string holding the name of the variable 
#          it should be used for post-processing
# pos    - character specifying if the variable is cell centered (C),
#          staggered in x direction (X) or in y direction (Y)
# res    - vector specifying resolutions in x and y direction
# def_bc - integer specifying if the default boundary condition is of 
#          dirichlet, neumann or outlet type
#--------------------------------------------------------------------------

  # Fetch resolutions 
  nx, ny, nz = res
  
  # Create boundary tuple
  key = namedtuple('key', 'typ val')  
  bnd = (key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),  \
         key(ndarray(shape=(1,ny,nz), dtype=int), zeros((1,ny,nz))),  \
         key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),  \
         key(ndarray(shape=(nx,1,nz), dtype=int), zeros((nx,1,nz))),  \
         key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))),  \
         key(ndarray(shape=(nx,ny,1), dtype=int), zeros((nx,ny,1))))
  
  # Create the unknown tuple
  unk = namedtuple('unk', 'name pos val old bnd')
  phi = unk(name, pos, zeros(res), zeros(res), bnd)
  
  #------------------
  # Set the position
  #------------------
  if pos != C and pos != X and pos != Y and pos != Z: 
    print('Variable must be defined at positions C, X, Y or Z')

  #---------------------------------------------------
  # Set default boundary conditions, types and values
  #---------------------------------------------------
  if def_bc != DIRICHLET and def_bc != NEUMANN and def_bc != OUTLET: 
    print('Variable must be defined with boundary conditions ',  \
          'DIRICHLET, NEUMANN or OUTLET!')

  phi.bnd[W].typ[0,:,:] = def_bc
  phi.bnd[E].typ[0,:,:] = def_bc
  phi.bnd[S].typ[:,0,:] = def_bc
  phi.bnd[N].typ[:,0,:] = def_bc
  phi.bnd[B].typ[:,:,0] = def_bc
  phi.bnd[T].typ[:,:,0] = def_bc

  phi.bnd[W].val[0,:,:] = 0
  phi.bnd[E].val[0,:,:] = 0 
  phi.bnd[S].val[:,0,:] = 0
  phi.bnd[N].val[:,0,:] = 0
  phi.bnd[B].val[:,:,0] = 0
  phi.bnd[T].val[:,:,0] = 0

  print("Created variable ", name)

  return phi  # end of function
