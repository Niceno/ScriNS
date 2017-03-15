# Standard Python modules
from standard import *

# ScriNS modules
from Constants.all import *
from Operators.all import *

#==========================================================================
def obst_zero_val(d, val, obst):
#--------------------------------------------------------------------------
# Set value to zero inside obstacle.
# 
# d    - position of the variable, C, X, Y or Z
# val  - value to be zeroed in obstacle
# obst - matrix holding positions of obstacle
#--------------------------------------------------------------------------

  if d == C:  
    val = val * lnot(obst)
    
  elif d==X:
    obst_x = mx(obst[:-1,:,:], obst[1:,:,:])
    val = val * lnot(obst_x)
    
  elif d==Y: 
    obst_y = mx(obst[:,:-1,:], obst[:,1:,:])
    val = val * lnot(obst_y)
    
  elif d==Z: 
    obst_z = mx(obst[:,:,:-1], obst[:,:,1:])
    val = val * lnot(obst_z)

  return val  # end of function