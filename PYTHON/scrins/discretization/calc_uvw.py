
from numpy import reshape
from scipy import prod
from scipy.sparse.linalg import bicgstab

from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.constants.solver import TOL
from scrins.discretization.adj_o_bnds     import adj_o_bnds
from scrins.discretization.advection      import advection
from scrins.discretization.create_matrix  import create_matrix
from scrins.discretization.obst_zero_val  import obst_zero_val
from scrins.operators.avg import avg
from scrins.operators.dif import dif



#==========================================================================
def calc_uvw(uvw, uvwf, rho, mu, p_tot, e_f, dt, dxyz, obst):
#--------------------------------------------------------------------------

  # Unpack tuples
  u,   v,   w   = uvw
  uf,  vf,  wf  = uvwf
  dx,  dy,  dz  = dxyz
  e_x, e_y, e_z = e_f

  # Fetch resolutions
  ru = u.val.shape
  rv = v.val.shape
  rw = w.val.shape
  
  # Pre-compute geometrical quantities
  dv = dx * dy * dz

  d = u.pos

  # Create linear systems
  A_u, b_u = create_matrix(u, rho/dt, mu, dxyz, obst, 'd')
  A_v, b_v = create_matrix(v, rho/dt, mu, dxyz, obst, 'd')
  A_w, b_w = create_matrix(w, rho/dt, mu, dxyz, obst, 'd')
 
  # Advection terms for momentum                            
  c_u = advection(rho, u, uvwf, dxyz, dt, 'superbee');
  c_v = advection(rho, v, uvwf, dxyz, dt, 'superbee');
  c_w = advection(rho, w, uvwf, dxyz, dt, 'superbee');

  # Innertial term for momentum (this works for collocated and staggered)
  i_u = u.old * avg(u.pos, rho) * avg(u.pos, dv) / dt
  i_v = v.old * avg(v.pos, rho) * avg(v.pos, dv) / dt
  i_w = w.old * avg(w.pos, rho) * avg(w.pos, dv) / dt
  
  # Compute staggered pressure gradients
  p_tot_x = dif(X, p_tot) / avg(X, dx)
  p_tot_y = dif(Y, p_tot) / avg(Y, dy)
  p_tot_z = dif(Z, p_tot) / avg(Z, dz)

  # Make pressure gradients cell-centered
  if d == C:
    p_tot_x = avg(X, cat(X, (p_tot_x[:1,:,:], p_tot_x, p_tot_x[-1:,:,:])))
    p_tot_y = avg(Y, cat(Y, (p_tot_y[:,:1,:], p_tot_y, p_tot_y[:,-1:,:])))
    p_tot_z = avg(Z, cat(Z, (p_tot_z[:,:,:1], p_tot_z, p_tot_z[:,:,-1:])))

  # Total pressure gradients (this works for collocated and staggered)
  p_st_u = p_tot_x * avg(u.pos, dv)
  p_st_v = p_tot_y * avg(v.pos, dv)
  p_st_w = p_tot_z * avg(w.pos, dv)
  
  # Full force terms for momentum equations (collocated and staggered)
  f_u = b_u - c_u + i_u - p_st_u + e_x * avg(u.pos, dv)
  f_v = b_v - c_v + i_v - p_st_v + e_y * avg(v.pos, dv)
  f_w = b_w - c_w + i_w - p_st_w + e_z * avg(w.pos, dv)

  # Take care of obsts in the domian
  if obst.any() != 0:
    f_u = obst_zero_val(u.pos, f_u, obst)
    f_v = obst_zero_val(v.pos, f_v, obst)
    f_w = obst_zero_val(w.pos, f_w, obst)

  # Solve for velocities
  res0 = bicgstab( A_u, reshape(f_u, prod(ru)), tol=TOL )
  res1 = bicgstab( A_v, reshape(f_v, prod(rv)), tol=TOL )
  res2 = bicgstab( A_w, reshape(f_w, prod(rw)), tol=TOL )
  u.val[:] = reshape(res0[0], ru)
  v.val[:] = reshape(res1[0], rv)
  w.val[:] = reshape(res2[0], rw)

  # Update velocities in boundary cells
  adj_o_bnds((u,v,w), (dx,dy,dz), dt)

  # Update face velocities (also substract cell-centered pressure gradients 
  #                         and add staggered pressure gradients)
  if d == C:
    uf.val[:] = avg(X,u.val + dt /       rho  * (      p_tot_x     ))       \
                            - dt / avg(X,rho) * (dif(X,p_tot) / avg(X,dx))  
    vf.val[:] = avg(Y,v.val + dt /       rho  * (      p_tot_y     ))       \
                            - dt / avg(Y,rho) * (dif(Y,p_tot) / avg(Y,dy))  
    wf.val[:] = avg(Z,w.val + dt /       rho  * (      p_tot_z     ))       \
                            - dt / avg(Z,rho) * (dif(Z,p_tot) / avg(Z,dz))  

    for j in (W,E):
      uf.bnd[j].val[:] = u.bnd[j].val[:]  
      vf.bnd[j].val[:] = avg(Y, v.bnd[j].val[:])
      wf.bnd[j].val[:] = avg(Z, w.bnd[j].val[:])  
    for j in (S,N):
      uf.bnd[j].val[:] = avg(X, u.bnd[j].val[:])
      vf.bnd[j].val[:] = v.bnd[j].val[:]
      wf.bnd[j].val[:] = avg(Z, w.bnd[j].val[:])  
    for j in (B,T):  
      uf.bnd[j].val[:] = avg(X, u.bnd[j].val[:])
      vf.bnd[j].val[:] = avg(Y, v.bnd[j].val[:])
      wf.bnd[j].val[:] = w.bnd[j].val[:]  

  else:
    uf.val[:] = u.val[:]
    vf.val[:] = v.val[:]
    wf.val[:] = w.val[:]
    for j in (W,E,S,N,B,T):
      uf.bnd[j].val[:] = u.bnd[j].val[:]
      vf.bnd[j].val[:] = v.bnd[j].val[:]
      wf.bnd[j].val[:] = w.bnd[j].val[:]

  if obst.any() != 0:
    uf.val[:] = obst_zero_val(X, uf.val, obst)
    vf.val[:] = obst_zero_val(Y, vf.val, obst)
    wf.val[:] = obst_zero_val(Z, wf.val, obst)

  return  # end of function
