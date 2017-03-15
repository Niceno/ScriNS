# Standard Python modules
from standard import *

# ScriNS modules
from constants.all import *
from operators.all import *

#==========================================================================
def adj_o_bnds(uvw, dxyz, dt):
#--------------------------------------------------------------------------

  # Unpack tuples  
  u,  v,  w  = uvw
  dx, dy, dz = dxyz

  # Local variables used in this function 
  area_in   = 0.0  # area of the inlet
  area_out  = 0.0  # area of the outlet 
  vol_in    = 0.0  # inlet volume flux; positive for inflow
  vol_out_1 = 0.0  # outlet volume flux; positive for outflow
  vol_out_2 = 0.0  # outlet volume flux; positive for outflow
  
  verbatim = False
  
  sx = dy * dz
  sy = dx * dz
  sz = dx * dy
  
  sx = sx[:1,:,:]
  sy = sy[:,:1,:]
  sz = sz[:,:,:1]
  
  #------------------------------------------------------------------
  # Compute the volume flowing in (v_in), volume flowing out (v_out) 
  # as well as inlet and outlet areas (a_in, a_out)
  #------------------------------------------------------------------

  # Inlets: these arrays will hold values true (1) in cells \
  # with inlet boundary conditions, and false (0) otherwise
  if_w_in = (u.bnd[W].typ[:1,:,:]==DIRICHLET) & (u.bnd[W].val[:1,:,:]>+TINY)
  if_e_in = (u.bnd[E].typ[:1,:,:]==DIRICHLET) & (u.bnd[E].val[:1,:,:]<-TINY)
  if_s_in = (v.bnd[S].typ[:,:1,:]==DIRICHLET) & (v.bnd[S].val[:,:1,:]>+TINY)
  if_n_in = (v.bnd[N].typ[:,:1,:]==DIRICHLET) & (v.bnd[N].val[:,:1,:]<-TINY)
  if_b_in = (w.bnd[B].typ[:,:,:1]==DIRICHLET) & (w.bnd[B].val[:,:,:1]>+TINY)
  if_t_in = (w.bnd[T].typ[:,:,:1]==DIRICHLET) & (w.bnd[T].val[:,:,:1]<-TINY)
  
  # Using the arrays defined above, compute inlet surface area
  area_in += (if_w_in * sx).sum()
  area_in += (if_e_in * sx).sum()
  area_in += (if_s_in * sy).sum()
  area_in += (if_n_in * sy).sum()
  area_in += (if_b_in * sz).sum()
  area_in += (if_t_in * sz).sum()
  
  # If there is no inlet, nothing to do here any longer
  if area_in < TINY:
    return u, v, w  # one end of function
  
  # Using the arrays defined above, compute inlet volume flux
  vol_in += (if_w_in * u.bnd[W].val[:1,:,:] * sx).sum()
  vol_in -= (if_e_in * u.bnd[E].val[:1,:,:] * sx).sum()
  vol_in += (if_s_in * v.bnd[S].val[:,:1,:] * sy).sum()
  vol_in -= (if_n_in * v.bnd[N].val[:,:1,:] * sy).sum()
  vol_in += (if_b_in * w.bnd[B].val[:,:,:1] * sz).sum()
  vol_in -= (if_t_in * w.bnd[T].val[:,:,:1] * sz).sum()
  
  # Outlets: these arrays will hold values true (1) in cells ...
  # with outlet boundary conditions, and false (0) otherwise
  if_w_out_u = ( u.bnd[W].typ[:1,:,:] == OUTLET )
  if_w_out_v = avg(v.pos, if_w_out_u*1) > 0.5
  if_w_out_w = avg(w.pos, if_w_out_u*1) > 0.5

  if_e_out_u = ( u.bnd[E].typ[:1,:,:] == OUTLET )
  if_e_out_v = avg(v.pos, if_e_out_u*1) > 0.5
  if_e_out_w = avg(w.pos, if_e_out_u*1) > 0.5

  if_s_out_v = ( v.bnd[S].typ[:,:1,:] == OUTLET )
  if_s_out_u = avg(u.pos, if_s_out_v*1) > 0.5
  if_s_out_w = avg(w.pos, if_s_out_v*1) > 0.5

  if_n_out_v = ( v.bnd[N].typ[:,:1,:] == OUTLET )
  if_n_out_u = avg(u.pos, if_n_out_v*1) > 0.5
  if_n_out_w = avg(w.pos, if_n_out_v*1) > 0.5

  if_b_out_w = ( w.bnd[B].typ[:,:,:1] == OUTLET )
  if_b_out_u = avg(u.pos, if_b_out_w*1) > 0.5
  if_b_out_v = avg(v.pos, if_b_out_w*1) > 0.5

  if_t_out_w = ( w.bnd[T].typ[:,:,:1] == OUTLET )
  if_t_out_u = avg(u.pos, if_t_out_w*1) > 0.5
  if_t_out_v = avg(v.pos, if_t_out_w*1) > 0.5
                    
  # Using the arrays defined above, compute outlet surface area
  area_out += (if_w_out_u * sx).sum()
  area_out += (if_e_out_u * sx).sum()
  area_out += (if_s_out_v * sy).sum()
  area_out += (if_n_out_v * sy).sum()
  area_out += (if_b_out_w * sz).sum()
  area_out += (if_t_out_w * sz).sum()
  
  # Using the arrays defined above, compute outlet volume flux
  vol_out_1 -= (if_w_out_u * u.bnd[W].val[:1,:,:] * sx).sum()
  vol_out_1 += (if_e_out_u * u.bnd[E].val[:1,:,:] * sx).sum()
  vol_out_1 -= (if_s_out_v * v.bnd[S].val[:,:1,:] * sy).sum()
  vol_out_1 += (if_n_out_v * v.bnd[N].val[:,:1,:] * sy).sum()
  vol_out_1 -= (if_b_out_w * w.bnd[B].val[:,:,:1] * sz).sum()
  vol_out_1 += (if_t_out_w * w.bnd[T].val[:,:,:1] * sz).sum()
  
  #---------------------------------
  # Check and calculate corrections
  #---------------------------------
  
  if (area_in == 0):
    ub_in = 0
  else:
    ub_in  = vol_in / area_in
  
  if (area_out == 0):
    ub_out = 0
  else:
    ub_out = vol_out_1 / area_out
  
  #----------------------------------------------
  # If nothing comes out, make a bulk correction
  #----------------------------------------------
  if(ub_out < TINY):  
    u_bulk_corr = ub_in * area_in / area_out
    
    u.bnd[W].val[:1,:,:] = u.bnd[W].val[:1,:,:] * lnot(if_w_out_u) \
                         - u_bulk_corr          *      if_w_out_u
    u.bnd[E].val[:1,:,:] = u.bnd[E].val[:1,:,:] * lnot(if_e_out_u) \
                         + u_bulk_corr          *      if_e_out_u
    v.bnd[S].val[:,:1,:] = v.bnd[S].val[:,:1,:] * lnot(if_s_out_v) \
                         - u_bulk_corr          *      if_s_out_v
    v.bnd[N].val[:,:1,:] = v.bnd[N].val[:,:1,:] * lnot(if_n_out_v) \
                         + u_bulk_corr          *      if_n_out_v
    w.bnd[B].val[:,:,:1] = w.bnd[B].val[:,:,:1] * lnot(if_b_out_w) \
                         - u_bulk_corr          *      if_b_out_w
    w.bnd[T].val[:,:,:1] = w.bnd[T].val[:,:,:1] * lnot(if_t_out_w) \
                         + u_bulk_corr          *      if_t_out_w
                         
  #--------------------------------------------------------------
  # Correction outflow by applying convective boundary condition
  #--------------------------------------------------------------
  else:
    du_dx_w = (u.val[ :1,:,:]-u.bnd[W].val[:1,:,:]) / dx[ :1,:,:]
    dv_dx_w = (v.val[ :1,:,:]-v.bnd[W].val[:1,:,:]) / avg(v.pos, dx[ :1,:,:])
    dw_dx_w = (w.val[ :1,:,:]-w.bnd[W].val[:1,:,:]) / avg(w.pos, dx[ :1,:,:])
  
    du_dx_e = (u.val[-1:,:,:]-u.bnd[E].val[:1,:,:]) / dx[-1:,:,:]
    dv_dx_e = (v.val[-1:,:,:]-v.bnd[E].val[:1,:,:]) / avg(v.pos, dx[-1:,:,:])
    dw_dx_e = (w.val[-1:,:,:]-w.bnd[E].val[:1,:,:]) / avg(w.pos, dx[-1:,:,:])
  
    du_dy_s = (u.val[:, :1,:]-u.bnd[S].val[:,:1,:]) / avg(u.pos, dy[:, :1,:])
    dv_dy_s = (v.val[:, :1,:]-v.bnd[S].val[:,:1,:]) / dy[:, :1,:]
    dw_dy_s = (w.val[:, :1,:]-w.bnd[S].val[:,:1,:]) / avg(w.pos, dy[:, :1,:])
  
    du_dy_n = (u.val[:,-1:,:]-u.bnd[N].val[:,:1,:]) / avg(u.pos, dy[:,-1:,:])
    dv_dy_n = (v.val[:,-1:,:]-v.bnd[N].val[:,:1,:]) / dy[:,-1:,:]
    dw_dy_n = (w.val[:,-1:,:]-w.bnd[N].val[:,:1,:]) / avg(w.pos, dy[:,-1:,:])
  
    du_dz_b = (u.val[:,:, :1]-u.bnd[B].val[:,:,:1]) / avg(u.pos, dz[:,:, :1])
    dv_dz_b = (v.val[:,:, :1]-v.bnd[B].val[:,:,:1]) / avg(v.pos, dz[:,:, :1])
    dw_dz_b = (w.val[:,:, :1]-w.bnd[B].val[:,:,:1]) / dz[:,:, :1]
  
    du_dz_t = (u.val[:,:,-1:]-u.bnd[T].val[:,:,:1]) / avg(u.pos, dz[:,:,-1:])
    dv_dz_t = (v.val[:,:,-1:]-v.bnd[T].val[:,:,:1]) / avg(v.pos, dz[:,:,-1:])
    dw_dz_t = (w.val[:,:,-1:]-w.bnd[T].val[:,:,:1]) / dz[:,:,-1:]
 
    u_bnd_w_corr = (u.bnd[W].val[:1,:,:] + ub_out * dt * du_dx_w)
    v_bnd_w_corr = (v.bnd[W].val[:1,:,:] + ub_out * dt * dv_dx_w)
    w_bnd_w_corr = (w.bnd[W].val[:1,:,:] + ub_out * dt * dw_dx_w)
  
    u.bnd[W].val[:1,:,:] = u.bnd[W].val[:1,:,:] * lnot(if_w_out_u) \
                         + u_bnd_w_corr         *      if_w_out_u
    v.bnd[W].val[:1,:,:] = v.bnd[W].val[:1,:,:] * lnot(if_w_out_v) \
                         + v_bnd_w_corr         *      if_w_out_v
    w.bnd[W].val[:1,:,:] = w.bnd[W].val[:1,:,:] * lnot(if_w_out_w) \
                         + w_bnd_w_corr         *      if_w_out_w
  
    u_bnd_e_corr = (u.bnd[E].val[:1,:,:] + ub_out * dt * du_dx_e)
    v_bnd_e_corr = (v.bnd[E].val[:1,:,:] + ub_out * dt * dv_dx_e)
    w_bnd_e_corr = (w.bnd[E].val[:1,:,:] + ub_out * dt * dw_dx_e)
  
    u.bnd[E].val[:1,:,:] = u.bnd[E].val[:1,:,:] * lnot(if_e_out_u) \
                         + u_bnd_e_corr         *      if_e_out_u
    v.bnd[E].val[:1,:,:] = v.bnd[E].val[:1,:,:] * lnot(if_e_out_v) \
                         + v_bnd_e_corr         *      if_e_out_v
    w.bnd[E].val[:1,:,:] = w.bnd[E].val[:1,:,:] * lnot(if_e_out_w) \
                         + w_bnd_e_corr         *      if_e_out_w
  
    u_bnd_s_corr = (u.bnd[S].val[:,:1,:] + ub_out * dt * du_dy_s)
    v_bnd_s_corr = (v.bnd[S].val[:,:1,:] + ub_out * dt * dv_dy_s)
    w_bnd_s_corr = (w.bnd[S].val[:,:1,:] + ub_out * dt * dw_dy_s)
  
    u.bnd[S].val[:,:1,:] = u.bnd[S].val[:,:1,:] * lnot(if_s_out_u) \
                         + u_bnd_s_corr         *      if_s_out_u
    v.bnd[S].val[:,:1,:] = v.bnd[S].val[:,:1,:] * lnot(if_s_out_v) \
                         + v_bnd_s_corr         *      if_s_out_v
    w.bnd[S].val[:,:1,:] = w.bnd[S].val[:,:1,:] * lnot(if_s_out_w) \
                         + w_bnd_s_corr         *      if_s_out_w
  
    u_bnd_n_corr = (u.bnd[N].val[:,:1,:] + ub_out * dt * du_dy_n)
    v_bnd_n_corr = (v.bnd[N].val[:,:1,:] + ub_out * dt * dv_dy_n)
    w_bnd_n_corr = (w.bnd[N].val[:,:1,:] + ub_out * dt * dw_dy_n)
  
    u.bnd[N].val[:,:1,:] = u.bnd[N].val[:,:1,:] * lnot(if_n_out_u) \
                         + u_bnd_n_corr         *     if_n_out_u
    v.bnd[N].val[:,:1,:] = v.bnd[N].val[:,:1,:] * lnot(if_n_out_v) \
                         + v_bnd_n_corr         *      if_n_out_v
    w.bnd[N].val[:,:1,:] = w.bnd[N].val[:,:1,:] * lnot(if_n_out_w) \
                         + w_bnd_n_corr         *      if_n_out_w
  
    u_bnd_b_corr = (u.bnd[B].val[:,:,:1] + ub_out * dt * du_dz_b)
    v_bnd_b_corr = (v.bnd[B].val[:,:,:1] + ub_out * dt * dv_dz_b)
    w_bnd_b_corr = (w.bnd[B].val[:,:,:1] + ub_out * dt * dw_dz_b)
  
    u.bnd[B].val[:,:,:1] = u.bnd[B].val[:,:,:1] * lnot(if_b_out_u) \
                         + u_bnd_b_corr         *      if_b_out_u
    v.bnd[B].val[:,:,:1] = v.bnd[B].val[:,:,:1] * lnot(if_b_out_v) \
                         + v_bnd_b_corr         *      if_b_out_v
    w.bnd[B].val[:,:,:1] = w.bnd[B].val[:,:,:1] * lnot(if_b_out_w) \
                         + w_bnd_b_corr         *      if_b_out_w
  
    u_bnd_t_corr = (u.bnd[T].val[:,:,:1] + ub_out *dt * du_dz_t)
    v_bnd_t_corr = (v.bnd[T].val[:,:,:1] + ub_out *dt * dv_dz_t)
    w_bnd_t_corr = (w.bnd[T].val[:,:,:1] + ub_out *dt * dw_dz_t)
  
    u.bnd[T].val[:,:,:1] = u.bnd[T].val[:,:,:1] * lnot(if_t_out_u) \
                         + u_bnd_t_corr         *      if_t_out_u
    v.bnd[T].val[:,:,:1] = v.bnd[T].val[:,:,:1] * lnot(if_t_out_v) \
                         + v_bnd_t_corr         *      if_t_out_v
    w.bnd[T].val[:,:,:1] = w.bnd[T].val[:,:,:1] * lnot(if_t_out_w) \
                         + w_bnd_t_corr         *      if_t_out_w
  
  if verbatim == True:
    print('+----------------------------+'     )
    print('|  ub_in     = %12.5e  |' %ub_in    )
    print('|  a_in      = %12.5e  |' %area_in  )
    print('|  v_in      = %12.5e  |' %vol_in   )
    print('|  ub_out    = %12.5e  |' %ub_out   )
    print('|  a_out     = %12.5e  |' %area_out )
    print('|  v_out_1   = %12.5e  |' %vol_out_1)
    
  #----------------------------------------------
  # Scaling correction to whatever you did above 
  # (bulk correction or convective outflow)
  #----------------------------------------------
  vol_out_2 = 0.0
  vol_out_2 -= (if_w_out_u * u.bnd[W].val[:1,:,:] * sx).sum()
  vol_out_2 += (if_e_out_u * u.bnd[E].val[:1,:,:] * sx).sum()
  vol_out_2 -= (if_s_out_v * v.bnd[S].val[:,:1,:] * sy).sum()
  vol_out_2 += (if_n_out_v * v.bnd[N].val[:,:1,:] * sy).sum()
  vol_out_2 -= (if_b_out_w * w.bnd[B].val[:,:,:1] * sz).sum()
  vol_out_2 += (if_t_out_w * w.bnd[T].val[:,:,:1] * sz).sum()
  
  if vol_out_2 > TINY:
    factor = vol_in / vol_out_2
  else:
    factor = 1.0  
  
  if verbatim == True:
    print('+----------------------------+') 
    print('|  v_out_2   = %12.5e  |' %vol_out_2)
    print('|  factor    = %12.5e  |' %factor   )
    print('+----------------------------+') 
  
  #--------------------------------------
  # Correction to satisfy volume balance
  #--------------------------------------
  u.bnd[W].val[:1,:,:] = u.bnd[W].val[:1,:,:] * lnot(if_w_out_u)          \
                       + u.bnd[W].val[:1,:,:] *      if_w_out_u * factor
  u.bnd[E].val[:1,:,:] = u.bnd[E].val[:1,:,:] * lnot(if_e_out_u)          \
                       + u.bnd[E].val[:1,:,:] *      if_e_out_u * factor
  v.bnd[S].val[:,:1,:] = v.bnd[S].val[:,:1,:] * lnot(if_s_out_v)          \
                       + v.bnd[S].val[:,:1,:] *      if_s_out_v * factor
  v.bnd[N].val[:,:1,:] = v.bnd[N].val[:,:1,:] * lnot(if_n_out_v)          \
                       + v.bnd[N].val[:,:1,:] *      if_n_out_v * factor
  w.bnd[B].val[:,:,:1] = w.bnd[B].val[:,:,:1] * lnot(if_b_out_w)          \
                       + w.bnd[B].val[:,:,:1] *      if_b_out_w * factor
  w.bnd[T].val[:,:,:1] = w.bnd[T].val[:,:,:1] * lnot(if_t_out_w)          \
                       + w.bnd[T].val[:,:,:1] *      if_t_out_w * factor

  return  # another end of function
