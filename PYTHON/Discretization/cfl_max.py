# ScriNS modules
from Constants.all      import *
from Operators.all      import *

#==========================================================================
def cfl_max(uvw, dt, dxyz):
#--------------------------------------------------------------------------

  # Unpack received tuples
  u,  v,  w  = uvw
  dx, dy, dz = dxyz

  # Take velocity's position
  d = u.pos
  
  # Mesh is cell-centered
  if d == C:
    cfl = dt * max( abs(u.val/dx).max(),   \
                    abs(v.val/dy).max(),   \
                    abs(w.val/dz).max() )
  
  # Mesh is staggered
  else:
    cfl = dt * max( abs(u.val/avg(X,dx)).max(),   \
                    abs(v.val/avg(Y,dy)).max(),   \
                    abs(w.val/avg(Z,dz)).max() )
  
  return cfl