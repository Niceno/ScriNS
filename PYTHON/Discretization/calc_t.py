# Standard Python modules
from standard import *

# ScriNS modules
from Constants.all      import *
from Operators.all      import *

from Discretization.adj_n_bnds     import adj_n_bnds
from Discretization.advection      import advection
from Discretization.create_matrix  import create_matrix

#==========================================================================
def calc_t(t, uvwf, rho_cap, kappa, dt, dxyz, obst):
#--------------------------------------------------------------------------

  # Unpack tuple(s)
  dx, dy, dz = dxyz

  # Fetch the resolution  
  rc = t.val.shape 

  # Discretize the diffusive part
  A_t, b_t = create_matrix(t, rho_cap/dt, kappa, dxyz, obst, 'n')

  # The advective fluxes 
  c_t = advection(rho_cap, t, uvwf, dxyz, dt, 'minmod')
  
  # Innertial term for enthalpy
  i_t = t.old * rho_cap * dx*dy*dz / dt

  # The entire source term
  f_t = b_t - c_t + i_t

  # Solve for temperature 
  res0 = bicgstab( A_t, reshape(f_t, prod(rc)), tol=TOL )
  t.val[:] = reshape(res0[0], rc) 

  adj_n_bnds(t)

  return  # end of function
