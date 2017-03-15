# Standard Python modules
from standard import *

# ScriNS modules
from constants.all      import *
from operators.all      import *

#==========================================================================
def nodes(*args):
#--------------------------------------------------------------------------
# Creates node coordinates in one dimension
#
# Input arguments:
#   s     - starting coordinate
#   e     - end coordinate
#   n     - number of cells (number of nodes is greater than one than this)
#   del_s - (optional) starting cell size
#   del_e - (optional) ending cell size 
#--------------------------------------------------------------------------
 
  #---------------------------------------------
  # Just a simple constant-spacing distribution
  #---------------------------------------------
  if len((args)) == 3:
    s = args[0]
    e = args[1]
    n = args[2]
    return linspace(s, e, n+1)

  #------------------------------------------------------------------------
  # Grid distribution is determined from the following function:
  #   x = a + b y + c y^2 + d y^3
  #
  # One should imagine this function as taking integer arguments.
  #
  # Boundary conditions:
  #
  #   x(0)   = s         => a                                   = s
  #   x(1)   = s + del_s => a + b       + c          + d        = s + del_s
  #   x(n-1) = e - del_e => a + b*(n-1) + c*(n-1)^2 + d*(n-1)^3 = e - del_e
  #   x(n)   = e         => a + b*n     + c*n^2     + d*n^3     = e
  #
  # It follows that:
  #
  #   |1   0      0        0     |  |a|  |s      |
  #   |1   1      1        1     |  |b|  |s+del_s|
  #   |1   n-1   (n-1)^2  (n-1)^3|  |c|  |e-del_e|
  #   |1   n      n^2      n^3   |  |d|  |e      |
  #------------------------------------------------------------------------
  elif len((args)) == 5:
    s     = args[0]
    e     = args[1]
    n     = args[2]
    del_s = args[3]
    del_e = args[4]

    # Form the system matrix    
    A = matrix( [ [1,  0,    0,           0         ],     \
                  [1,  1,    1,           1         ],     \
                  [1,  n-1,  pow(n-1,2),  pow(n-1,3)],     \
                  [1,  n,    pow(n,  2),  pow(n,  3)] ] )
    
    # Form the right hand side
    f = array( [ s        ,     \
                 s + del_s,     \
                 e - del_e,     \
                 e         ] )
    
    # Solve for coefficients
    abcd = solve(A, f)
    
    # Pre-allocate the array for coordinates ...
    x = zeros(n+1);
             
    # ... and fill it up.         
    for i in range(0, n+1):
      x[i] =   abcd[0]          \
             + abcd[1] * i      \
             + abcd[2] * i*i    \
             + abcd[3] * i*i*i  
    
    return x
  
  #--------------------------
  # Some error message might 
  # be printed if number of 
  # arguments is wrong
  #--------------------------
    
  return  # end of function