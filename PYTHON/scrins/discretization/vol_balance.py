
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.obst_zero_val  import obst_zero_val
from scrins.operators.cat import cat
from scrins.operators.dif import dif

#==========================================================================
def vol_balance(uvwf, dxyz, obst):
#--------------------------------------------------------------------------
# Computes the volume balance, which is essentially a right hand side in
# the Poisson's equation for pressure.
#
# Note that "obst" is an optional parameter.  If it is not sent, the source
# won't be zeroed inside the obstacle.  That is important for calculation
# of pressure, see function "calc_p".
#--------------------------------------------------------------------------

  # Unpack tuples
  dx, dy, dz = dxyz
  uf, vf, wf = uvwf

  # Compute it throughout the domain
  src = - dif(X, cat(X, (uf.bnd[W].val, uf.val, uf.bnd[E].val)))*dy*dz  \
        - dif(Y, cat(Y, (vf.bnd[S].val, vf.val, vf.bnd[N].val)))*dx*dz  \
        - dif(Z, cat(Z, (wf.bnd[B].val, wf.val, wf.bnd[T].val)))*dx*dy   
  
  # Zero it inside obstacles, if obstacle is sent as parameter
  if obst.any() != 0:
    src = obst_zero_val(C, src, obst)
  
  return src  # end of function