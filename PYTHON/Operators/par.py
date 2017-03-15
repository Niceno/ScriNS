# Standard Python modules
from standard import *

# ScriNS modules
from Constants.all import *
from Display.all   import *

from Operators.all import avg

#==========================================================================
def par(mean_val, x_nodes):
#--------------------------------------------------------------------------
# A function to generate a parabolic profile over a set of cell centers,
# useful for specifying parabolic inlet velocity profiles.  It expects 
# nodal coordinates as input, but sends values at cell centers back.
#
# Input coordinates:    |-----|-----|-----|-----|-----|-----|-----|  
# Output values:           o-----o-----o-----o-----o-----o-----o
#
# Input parameters
#
# mean_val - mean value the parabola will have
# x_nodes  - nodal coordinatels
#--------------------------------------------------------------------------

  # It is known that maximum of a parabola is 3/2 of its mean value
  max_val = mean_val * 3/2

  # Normalized x coordinates (from -1 to +1)
  xn = copy(x_nodes)
  
  xn -= xn.min()                  
  xn /= (xn.max()-xn.min())
  xn *= 2                               
  xn -= 1  
  
  xc = avg(xn)

  yc = (1.0-xc*xc) * max_val

  return yc  # end of function