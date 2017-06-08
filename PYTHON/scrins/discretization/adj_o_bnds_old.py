# Standard Python modules
from standard import *

# ScriNS modules
from constants.all import *

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
  
  verbatim = True
  
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

  if_w_in = (u.bnd[W].typ[:1,:,:] == DIRICHLET)  \
          & (u.bnd[W].val[:1,:,:] > +TINY)
  # Inlets: these arrays will hold values true (1) in cells \
  # with inlet boundary conditions, and false (0) otherwise
  if_w_in = (u.bnd[W].typ[:1,:,:]==DIRICHLET) & (u.bnd[W].val[:1,:,:]>+TINY)
  if_e_in = (u.bnd[E].typ[:1,:,:]==DIRICHLET) & (u.bnd[E].val[:1,:,:]<-TINY)
  if_s_in = (v.bnd[S].typ[:,:1,:]==DIRICHLET) & (v.bnd[S].val[:,:1,:]>+TINY)
  if_n_in = (v.bnd[N].typ[:,:1,:]==DIRICHLET) & (v.bnd[N].val[:,:1,:]<-TINY)
  if_b_in = (w.bnd[B].typ[:,:,:1]==DIRICHLET) & (w.bnd[B].val[:,:,:1]>+TINY)
  if_t_in = (w.bnd[T].typ[:,:,:1]==DIRICHLET) & (w.bnd[T].val[:,:,:1]<-TINY)
  
  # Using the arrays defined above, compute inlet surface area
  area_in = area_in + (if_w_in * sx).sum()
  area_in = area_in + (if_e_in * sx).sum()
  area_in = area_in + (if_s_in * sy).sum()
  area_in = area_in + (if_n_in * sy).sum()
  area_in = area_in + (if_b_in * sz).sum()
  area_in = area_in + (if_t_in * sz).sum()
  
  # If there is no inlet, nothing to do here any longer
  if area_in < TINY:
    return u, v, w  # one end of function
  
  # Using the arrays defined above, compute inlet volume flux
  vol_in = vol_in + (if_w_in * u.bnd[W].val[:1,:,:] * sx).sum()
  vol_in = vol_in - (if_e_in * u.bnd[E].val[:1,:,:] * sx).sum()
  vol_in = vol_in + (if_s_in * v.bnd[S].val[:,:1,:] * sy).sum()
  vol_in = vol_in - (if_n_in * v.bnd[N].val[:,:1,:] * sy).sum()
  vol_in = vol_in + (if_b_in * w.bnd[B].val[:,:,:1] * sz).sum()
  vol_in = vol_in - (if_t_in * w.bnd[T].val[:,:,:1] * sz).sum()
  
  # Outlets: these arrays will hold values true (1) in cells ...
  # with outlet boundary conditions, and false (0) otherwise
  if_w_out = ( u.bnd[W].typ[:1,:,:] == OUTLET )
  if_e_out = ( u.bnd[E].typ[:1,:,:] == OUTLET )
  if_s_out = ( v.bnd[S].typ[:,:1,:] == OUTLET )
  if_n_out = ( v.bnd[N].typ[:,:1,:] == OUTLET )
  if_b_out = ( w.bnd[B].typ[:,:,:1] == OUTLET )
  if_t_out = ( w.bnd[T].typ[:,:,:1] == OUTLET )
  
  # Using the arrays defined above, compute outlet surface area
  area_out = area_out + (if_w_out * sx).sum()
  area_out = area_out + (if_e_out * sx).sum()
  area_out = area_out + (if_s_out * sy).sum()
  area_out = area_out + (if_n_out * sy).sum()
  area_out = area_out + (if_b_out * sz).sum()
  area_out = area_out + (if_t_out * sz).sum()
  
  # Using the arrays defined above, compute outlet volume flux
  vol_out_1 = vol_out_1 - (if_w_out * u.bnd[W].val[:1,:,:] * sx).sum()
  vol_out_1 = vol_out_1 + (if_e_out * u.bnd[E].val[:1,:,:] * sx).sum()
  vol_out_1 = vol_out_1 - (if_s_out * v.bnd[S].val[:,:1,:] * sy).sum()
  vol_out_1 = vol_out_1 + (if_n_out * v.bnd[N].val[:,:1,:] * sy).sum()
  vol_out_1 = vol_out_1 - (if_b_out * w.bnd[B].val[:,:,:1] * sz).sum()
  vol_out_1 = vol_out_1 + (if_t_out * w.bnd[T].val[:,:,:1] * sz).sum()
  
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
  
  bulk_corr = 0.0   # bulk correction to velocity
  corr      = 0.0   # switch between bulk correction and scaling
  
  if(ub_out < TINY):  # bulk correction makes sense if nothing comes out
    bulk_corr = ub_in * area_in / area_out
    corr      = 0.0   
  else:               # scaling factor makes sense if something comes out
    bulk_corr = 0.0
    corr      = 1.0
  
  if verbatim == True:
    print('+----------------------------+'     )
    print('|  ub_in     = %12.5e  |' %ub_in    )
    print('|  a_in      = %12.5e  |' %area_in  )
    print('|  v_in      = %12.5e  |' %vol_in   )
    print('|  ub_out    = %12.5e  |' %ub_out   )
    print('|  a_out     = %12.5e  |' %area_out )
    print('|  v_out_1   = %12.5e  |' %vol_out_1)
    print('|  bulk_corr = %12.5e  |' %bulk_corr)
    print('|  corr      = %12.5e  |' %corr     )

  #--------------------------------------------------------------
  # Correction outflow by applying convective boundary condition
  #--------------------------------------------------------------
  u_bnd_w_corr = - bulk_corr + corr * \
   (u.bnd[W].val[:1,:,:] + ub_out*dt*(u.val[:1,:,:]-u.bnd[W].val[:1,:,:])/dx[:1,:,:])
  u.bnd[W].val[:1,:,:] = u.bnd[W].val[:1,:,:] * lnot(if_w_out) \
                      + u_bnd_w_corr        *     if_w_out
  
  u_bnd_e_corr = + bulk_corr + corr * \
   (u.bnd[E].val[:1,:,:] + ub_out*dt*(u.val[-1:,:,:]-u.bnd[E].val[:1,:,:])/dx[-1:,:,:])
  u.bnd[E].val[:1,:,:] = u.bnd[E].val[:1,:,:] * lnot(if_e_out) \
                      + u_bnd_e_corr        *     if_e_out
  
  v_bnd_s_corr = - bulk_corr + corr * \
   (v.bnd[S].val[:,:1,:] + ub_out*dt*(v.val[:,:1,:]-v.bnd[S].val[:,:1,:])/dy[:,:1,:])
  v.bnd[S].val[:,:1,:] = v.bnd[S].val[:,:1,:] * lnot(if_s_out) \
                      + v_bnd_s_corr        *     if_s_out
  
  v_bnd_n_corr = + bulk_corr + corr * \
   (v.bnd[N].val[:,:1,:] + ub_out*dt*(v.val[:,-1:,:]-v.bnd[N].val[:,:1,:])/dy[:,-1:,:])
  v.bnd[N].val[:,:1,:] = v.bnd[N].val[:,:1,:] * lnot(if_n_out) \
                      + v_bnd_n_corr        *     if_n_out
  
  w_bnd_b_corr = - bulk_corr + corr * \
   (w.bnd[B].val[:,:,:1] + ub_out*dt*(w.val[:,:,:1]-w.bnd[B].val[:,:,:1])/dz[:,:,:1])
  w.bnd[B].val[:,:,:1] = w.bnd[B].val[:,:,:1] * lnot(if_b_out) \
                      + w_bnd_b_corr        *     if_b_out
  
  w_bnd_t_corr = + bulk_corr + corr * \
   (w.bnd[T].val[:,:,:1] + ub_out*dt*(w.val[:,:,-1:]-w.bnd[T].val[:,:,:1])/dz[:,:,-1:])
  w.bnd[T].val[:,:,:1] = w.bnd[T].val[:,:,:1] * lnot(if_t_out) \
                      + w_bnd_t_corr        *     if_t_out
  
  #----------------------------------------------
  # Scaling correction to whatever you did above 
  # (bulk correction or convective outflow)
  #----------------------------------------------
  vol_out_2 = 0.0
  vol_out_2 = vol_out_2 - (if_w_out * u.bnd[W].val[:1,:,:] * sx).sum()
  vol_out_2 = vol_out_2 + (if_e_out * u.bnd[E].val[:1,:,:] * sx).sum()
  vol_out_2 = vol_out_2 - (if_s_out * v.bnd[S].val[:,:1,:] * sy).sum()
  vol_out_2 = vol_out_2 + (if_n_out * v.bnd[N].val[:,:1,:] * sy).sum()
  vol_out_2 = vol_out_2 - (if_b_out * w.bnd[B].val[:,:,:1] * sz).sum()
  vol_out_2 = vol_out_2 + (if_t_out * w.bnd[T].val[:,:,:1] * sz).sum()
  
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
  u.bnd[W].val[:1,:,:] = u.bnd[W].val[:1,:,:] * lnot(if_w_out)          \
                       + u.bnd[W].val[:1,:,:] *      if_w_out * factor
  u.bnd[E].val[:1,:,:] = u.bnd[E].val[:1,:,:] * lnot(if_e_out)          \
                       + u.bnd[E].val[:1,:,:] *      if_e_out * factor
  v.bnd[S].val[:,:1,:] = v.bnd[S].val[:,:1,:] * lnot(if_s_out)          \
                       + v.bnd[S].val[:,:1,:] *      if_s_out * factor
  v.bnd[N].val[:,:1,:] = v.bnd[N].val[:,:1,:] * lnot(if_n_out)          \
                       + v.bnd[N].val[:,:1,:] *      if_n_out * factor
  w.bnd[B].val[:,:,:1] = w.bnd[B].val[:,:,:1] * lnot(if_b_out)          \
                       + w.bnd[B].val[:,:,:1] *      if_b_out * factor
  w.bnd[T].val[:,:,:1] = w.bnd[T].val[:,:,:1] * lnot(if_t_out)          \
                       + w.bnd[T].val[:,:,:1] *      if_t_out * factor

  return u, v, w  # another end of function
