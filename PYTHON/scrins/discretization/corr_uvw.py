
from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.discretization.obst_zero_val import obst_zero_val
from scrins.operators.avg import avg
from scrins.operators.dif import dif

#==========================================================================
def corr_uvw(uvw, p, rho, dt, dxyz, obst):
#--------------------------------------------------------------------------

  # Unpack received tuples
  dx, dy, dz = dxyz

  # Compute pressure correction gradients
  p_x = dif(X, p.val) / avg(X, dx)  
  p_y = dif(Y, p.val) / avg(Y, dy)  
  p_z = dif(Z, p.val) / avg(Z, dz)  
  
  # Set to zero in obst
  if obst.any() != 0:
    p_x = obst_zero_val(X, p_x, obst)
    p_y = obst_zero_val(Y, p_y, obst)
    p_z = obst_zero_val(Z, p_z, obst)
  
  # Pad with boundary values by expanding from interior 
  # (This is done only for collocated formulation)
  if uvw[X].pos == C:
    p_x = avg(X, cat(X, (p_x[:1,:,:], p_x, p_x[-1:,:,:])))
    p_y = avg(Y, cat(Y, (p_y[:,:1,:], p_y, p_y[:,-1:,:])))
    p_z = avg(Z, cat(Z, (p_z[:,:,:1], p_z, p_z[:,:,-1:])))
  
  # Correct the velocities
  uvw[X].val[:] = uvw[X].val[:] - dt / avg(uvw[X].pos, rho) * p_x
  uvw[Y].val[:] = uvw[Y].val[:] - dt / avg(uvw[Y].pos, rho) * p_y
  uvw[Z].val[:] = uvw[Z].val[:] - dt / avg(uvw[Z].pos, rho) * p_z
  
  return  # end of function