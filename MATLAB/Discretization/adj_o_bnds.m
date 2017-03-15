%==========================================================================
function [u, v, w] = adj_o_bnds(u, v, w, dx, dy, dz, dt)
%--------------------------------------------------------------------------
% Computes values at the outlet boundary using the convective outlet.
%--------------------------------------------------------------------------

% Local variables used in this function 
area_in   = 0.0;  % area of the inlet
area_out  = 0.0;  % area of the outlet 
vol_in    = 0.0;  % inlet volume flux; positive for inflow
vol_out_1 = 0.0;  % outlet volume flux; positive for outflow
vol_out_2 = 0.0;  % outlet volume flux; positive for outflow

verbatim = false;

sx = dy .* dz;
sy = dx .* dz;
sz = dx .* dy;

sx = sx(1,:,:);
sy = sy(:,1,:);
sz = sz(:,:,1);

%------------------------------------------------------------------
% Compute the volume flowing in (v_in), volume flowing out (v_out) 
% as well as inlet and outlet areas (a_in, a_out)
%------------------------------------------------------------------

% Inlets: these arrays will hold values true (1) in cells ...
% with inlet boundary conditions, and false (0) otherwise
if_w_in = ( u.bnd(W).type(1,:,:) == 'd' & u.bnd(W).val(1,:,:) > +TINY );
if_e_in = ( u.bnd(E).type(1,:,:) == 'd' & u.bnd(E).val(1,:,:) < -TINY );
if_s_in = ( v.bnd(S).type(:,1,:) == 'd' & v.bnd(S).val(:,1,:) > +TINY );
if_n_in = ( v.bnd(N).type(:,1,:) == 'd' & v.bnd(N).val(:,1,:) < -TINY );
if_b_in = ( w.bnd(B).type(:,:,1) == 'd' & w.bnd(B).val(:,:,1) > +TINY );
if_t_in = ( w.bnd(T).type(:,:,1) == 'd' & w.bnd(T).val(:,:,1) < -TINY );

% Using the arrays defined above, compute inlet surface area
area_in = area_in + sum( sum( if_w_in .* sx ) );
area_in = area_in + sum( sum( if_e_in .* sx ) );
area_in = area_in + sum( sum( if_s_in .* sy ) );
area_in = area_in + sum( sum( if_n_in .* sy ) );
area_in = area_in + sum( sum( if_b_in .* sz ) );
area_in = area_in + sum( sum( if_t_in .* sz ) );

% If there is no inlet, nothing to do here any longer
if area_in < TINY  
  return;
end

% Using the arrays defined above, compute inlet volume flux
vol_in = vol_in + sum(sum(if_w_in .* u.bnd(W).val(1,:,:) .* sx));
vol_in = vol_in - sum(sum(if_e_in .* u.bnd(E).val(1,:,:) .* sx));
vol_in = vol_in + sum(sum(if_s_in .* v.bnd(S).val(:,1,:) .* sy));
vol_in = vol_in - sum(sum(if_n_in .* v.bnd(N).val(:,1,:) .* sy));
vol_in = vol_in + sum(sum(if_b_in .* w.bnd(B).val(:,:,1) .* sz));
vol_in = vol_in - sum(sum(if_t_in .* w.bnd(T).val(:,:,1) .* sz));

% Outlets: these arrays will hold values true (1) in cells ...
% with outlet boundary conditions, and false (0) otherwise
if_w_out_x = ( u.bnd(W).type(1,:,:) == 'o' );
if_e_out_x = ( u.bnd(E).type(1,:,:) == 'o' );
if_s_out_y = ( v.bnd(S).type(:,1,:) == 'o' );
if_n_out_y = ( v.bnd(N).type(:,1,:) == 'o' );
if_b_out_z = ( w.bnd(B).type(:,:,1) == 'o' );
if_t_out_z = ( w.bnd(T).type(:,:,1) == 'o' );
if_w_out_y = avg(v.pos, if_w_out_x) > 0.5;
if_w_out_z = avg(w.pos, if_w_out_x) > 0.5;
if_e_out_y = avg(v.pos, if_e_out_x) > 0.5;
if_e_out_z = avg(w.pos, if_e_out_x) > 0.5;
if_s_out_x = avg(u.pos, if_s_out_y) > 0.5;
if_s_out_z = avg(w.pos, if_s_out_y) > 0.5;
if_n_out_x = avg(u.pos, if_n_out_y) > 0.5;
if_n_out_z = avg(w.pos, if_n_out_y) > 0.5;
if_b_out_x = avg(u.pos, if_b_out_z) > 0.5;
if_b_out_y = avg(v.pos, if_b_out_z) > 0.5;
if_t_out_x = avg(u.pos, if_t_out_z) > 0.5;
if_t_out_y = avg(v.pos, if_t_out_z) > 0.5;

% Using the arrays defined above, compute outlet surface area
area_out = area_out + sum( sum( if_w_out_x .* sx ) );
area_out = area_out + sum( sum( if_e_out_x .* sx ) );
area_out = area_out + sum( sum( if_s_out_y .* sy ) );
area_out = area_out + sum( sum( if_n_out_y .* sy ) );
area_out = area_out + sum( sum( if_b_out_z .* sz ) );
area_out = area_out + sum( sum( if_t_out_z .* sz ) );

% Using the arrays defined above, compute outlet volume flux
vol_out_1 = vol_out_1 - sum(sum(if_w_out_x .* u.bnd(W).val(1,:,:) .* sx));
vol_out_1 = vol_out_1 + sum(sum(if_e_out_x .* u.bnd(E).val(1,:,:) .* sx));
vol_out_1 = vol_out_1 - sum(sum(if_s_out_y .* v.bnd(S).val(:,1,:) .* sy));
vol_out_1 = vol_out_1 + sum(sum(if_n_out_y .* v.bnd(N).val(:,1,:) .* sy));
vol_out_1 = vol_out_1 - sum(sum(if_b_out_z .* w.bnd(B).val(:,:,1) .* sz));
vol_out_1 = vol_out_1 + sum(sum(if_t_out_z .* w.bnd(T).val(:,:,1) .* sz));

%---------------------------------
% Check and calculate corrections
%---------------------------------

if (area_in == 0)
  ub_in = 0;
else
  ub_in  = vol_in / area_in;
end

if (area_out == 0)
  ub_out = 0;
else
  ub_out = vol_out_1 / area_out;
end

%--------------------------------------------------
% Bulk correction makes sense if nothing comes out
%--------------------------------------------------
if(ub_out < TINY)  
  u_bulk_corr = ub_in * area_in / area_out;

  u.bnd(W).val(1,:,:) = u.bnd(W).val(1,:,:) .* not(if_w_out_x) ...
                      - u_bulk_corr         .*     if_w_out_x;
  u.bnd(E).val(1,:,:) = u.bnd(E).val(1,:,:) .* not(if_e_out_x) ...
                      + u_bulk_corr         .*     if_e_out_x;
  v.bnd(S).val(:,1,:) = v.bnd(S).val(:,1,:) .* not(if_s_out_y) ... 
                      - u_bulk_corr         .*     if_s_out_y;
  v.bnd(N).val(:,1,:) = v.bnd(N).val(:,1,:) .* not(if_n_out_y) ...
                      + u_bulk_corr         .*     if_n_out_y;
  w.bnd(B).val(:,:,1) = w.bnd(B).val(:,:,1) .* not(if_b_out_z) ...
                      - u_bulk_corr         .*     if_b_out_z;
  w.bnd(T).val(:,:,1) = w.bnd(T).val(:,:,1) .* not(if_t_out_z) ... 
                      + u_bulk_corr         .*     if_t_out_z; 
                  
%--------------------------------------------------------------
% Correction outflow by applying convective boundary condition
%--------------------------------------------------------------
else               
  du_dx_w = (u.val(  1,:,:) - u.bnd(W).val(1,:,:)) ./ dx(  1,:,:);
  dv_dx_w = (v.val(  1,:,:) - v.bnd(W).val(1,:,:)) ./ avg(v.pos, dx(  1,:,:));
  dw_dx_w = (w.val(  1,:,:) - w.bnd(W).val(1,:,:)) ./ avg(w.pos, dx(  1,:,:));

  du_dx_e = (u.val(end,:,:) - u.bnd(E).val(1,:,:)) ./ dx(end,:,:);
  dv_dx_e = (v.val(end,:,:) - v.bnd(E).val(1,:,:)) ./ avg(v.pos, dx(end,:,:));
  dw_dx_e = (w.val(end,:,:) - w.bnd(E).val(1,:,:)) ./ avg(w.pos, dx(end,:,:));

  du_dy_s = (u.val(:,  1,:) - u.bnd(S).val(:,1,:)) ./ avg(u.pos, dy(:,  1,:));
  dv_dy_s = (v.val(:,  1,:) - v.bnd(S).val(:,1,:)) ./ dy(:,  1,:);
  dw_dy_s = (w.val(:,  1,:) - w.bnd(S).val(:,1,:)) ./ avg(w.pos, dy(:,  1,:));

  du_dy_n = (u.val(:,end,:) - u.bnd(N).val(:,1,:)) ./ avg(u.pos, dy(:,end,:));
  dv_dy_n = (v.val(:,end,:) - v.bnd(N).val(:,1,:)) ./ dy(:,end,:);
  dw_dy_n = (w.val(:,end,:) - w.bnd(N).val(:,1,:)) ./ avg(w.pos, dy(:,end,:));

  du_dz_b = (u.val(:,:,  1) - u.bnd(B).val(:,:,1)) ./ avg(u.pos, dz(:,:,  1));
  dv_dz_b = (v.val(:,:,  1) - v.bnd(B).val(:,:,1)) ./ avg(v.pos, dz(:,:,  1));
  dw_dz_b = (w.val(:,:,  1) - w.bnd(B).val(:,:,1)) ./ dz(:,:,  1);

  du_dz_t = (u.val(:,:,end) - u.bnd(T).val(:,:,1)) ./ avg(u.pos, dz(:,:,end));
  dv_dz_t = (v.val(:,:,end) - v.bnd(T).val(:,:,1)) ./ avg(v.pos, dz(:,:,end));
  dw_dz_t = (w.val(:,:,end) - w.bnd(T).val(:,:,1)) ./ dz(:,:,end);

  u_bnd_w_corr = (u.bnd(W).val(1,:,:) + ub_out * dt * du_dx_w);
  v_bnd_w_corr = (v.bnd(W).val(1,:,:) + ub_out * dt * dv_dx_w);
  w_bnd_w_corr = (w.bnd(W).val(1,:,:) + ub_out * dt * dw_dx_w);

  u.bnd(W).val(1,:,:) = u.bnd(W).val(1,:,:) .* not(if_w_out_x) ...
                      + u_bnd_w_corr        .*     if_w_out_x;
  v.bnd(W).val(1,:,:) = v.bnd(W).val(1,:,:) .* not(if_w_out_y) ...
                      + v_bnd_w_corr        .*     if_w_out_y;
  w.bnd(W).val(1,:,:) = w.bnd(W).val(1,:,:) .* not(if_w_out_z) ...
                      + w_bnd_w_corr        .*     if_w_out_z;

  u_bnd_e_corr = (u.bnd(E).val(1,:,:) + ub_out * dt * du_dx_e);
  v_bnd_e_corr = (v.bnd(E).val(1,:,:) + ub_out * dt * dv_dx_e);
  w_bnd_e_corr = (w.bnd(E).val(1,:,:) + ub_out * dt * dw_dx_e);

  u.bnd(E).val(1,:,:) = u.bnd(E).val(1,:,:) .* not(if_e_out_x) ...
                      + u_bnd_e_corr        .*     if_e_out_x;
  v.bnd(E).val(1,:,:) = v.bnd(E).val(1,:,:) .* not(if_e_out_y) ...
                      + v_bnd_e_corr        .*     if_e_out_y;
  w.bnd(E).val(1,:,:) = w.bnd(E).val(1,:,:) .* not(if_e_out_z) ...
                      + w_bnd_e_corr        .*     if_e_out_z;

  u_bnd_s_corr = (u.bnd(S).val(:,1,:) + ub_out * dt * du_dy_s);
  v_bnd_s_corr = (v.bnd(S).val(:,1,:) + ub_out * dt * dv_dy_s);
  w_bnd_s_corr = (w.bnd(S).val(:,1,:) + ub_out * dt * dw_dy_s);

  u.bnd(S).val(:,1,:) = u.bnd(S).val(:,1,:) .* not(if_s_out_x) ... 
                      + u_bnd_s_corr        .*     if_s_out_x;
  v.bnd(S).val(:,1,:) = v.bnd(S).val(:,1,:) .* not(if_s_out_y) ... 
                      + v_bnd_s_corr        .*     if_s_out_y;
  w.bnd(S).val(:,1,:) = w.bnd(S).val(:,1,:) .* not(if_s_out_z) ... 
                      + w_bnd_s_corr        .*     if_s_out_z;

  u_bnd_n_corr = (u.bnd(N).val(:,1,:) + ub_out * dt * du_dy_n);
  v_bnd_n_corr = (v.bnd(N).val(:,1,:) + ub_out * dt * dv_dy_n);
  w_bnd_n_corr = (w.bnd(N).val(:,1,:) + ub_out * dt * dw_dy_n);

  u.bnd(N).val(:,1,:) = u.bnd(N).val(:,1,:) .* not(if_n_out_x) ...
                      + u_bnd_n_corr        .*     if_n_out_x;
  v.bnd(N).val(:,1,:) = v.bnd(N).val(:,1,:) .* not(if_n_out_y) ...
                      + v_bnd_n_corr        .*     if_n_out_y;
  w.bnd(N).val(:,1,:) = w.bnd(N).val(:,1,:) .* not(if_n_out_z) ...
                      + w_bnd_n_corr        .*     if_n_out_z;

  u_bnd_b_corr = (u.bnd(B).val(:,:,1) + ub_out * dt * du_dz_b);
  v_bnd_b_corr = (v.bnd(B).val(:,:,1) + ub_out * dt * dv_dz_b);
  w_bnd_b_corr = (w.bnd(B).val(:,:,1) + ub_out * dt * dw_dz_b);

  u.bnd(B).val(:,:,1) = u.bnd(B).val(:,:,1) .* not(if_b_out_x) ...
                      + u_bnd_b_corr        .*     if_b_out_x;
  v.bnd(B).val(:,:,1) = v.bnd(B).val(:,:,1) .* not(if_b_out_y) ...
                      + v_bnd_b_corr        .*     if_b_out_y;
  w.bnd(B).val(:,:,1) = w.bnd(B).val(:,:,1) .* not(if_b_out_z) ...
                      + w_bnd_b_corr        .*     if_b_out_z;

  u_bnd_t_corr = (u.bnd(T).val(:,:,1) + ub_out *dt * du_dz_t);
  v_bnd_t_corr = (v.bnd(T).val(:,:,1) + ub_out *dt * dv_dz_t);
  w_bnd_t_corr = (w.bnd(T).val(:,:,1) + ub_out *dt * dw_dz_t);

  u.bnd(T).val(:,:,1) = u.bnd(T).val(:,:,1) .* not(if_t_out_x) ... 
                      + u_bnd_t_corr        .*     if_t_out_x;
  v.bnd(T).val(:,:,1) = v.bnd(T).val(:,:,1) .* not(if_t_out_y) ... 
                      + v_bnd_t_corr        .*     if_t_out_y;
  w.bnd(T).val(:,:,1) = w.bnd(T).val(:,:,1) .* not(if_t_out_z) ... 
                      + w_bnd_t_corr        .*     if_t_out_z;
end

if verbatim == true
  disp( sprintf('+----------------------------+') );
  disp( sprintf('|  ub_in     = %12.5e  |', ub_in    ) );
  disp( sprintf('|  a_in      = %12.5e  |', area_in     ) );
  disp( sprintf('|  v_in      = %12.5e  |', vol_in     ) );
  disp( sprintf('|  ub_out    = %12.5e  |', ub_out   ) );
  disp( sprintf('|  a_out     = %12.5e  |', area_out    ) );
  disp( sprintf('|  v_out_1   = %12.5e  |', vol_out_1  ) );
end

%----------------------------------------------
% Scaling correction to whatever you did above 
% (bulk correction or convective outflow)
%----------------------------------------------
vol_out_2 = 0.0;
vol_out_2 = vol_out_2 - sum(sum(if_w_out_x .* u.bnd(W).val(1,:,:) .* sx));
vol_out_2 = vol_out_2 + sum(sum(if_e_out_x .* u.bnd(E).val(1,:,:) .* sx));
vol_out_2 = vol_out_2 - sum(sum(if_s_out_y .* v.bnd(S).val(:,1,:) .* sy));
vol_out_2 = vol_out_2 + sum(sum(if_n_out_y .* v.bnd(N).val(:,1,:) .* sy));
vol_out_2 = vol_out_2 - sum(sum(if_b_out_z .* w.bnd(B).val(:,:,1) .* sz));
vol_out_2 = vol_out_2 + sum(sum(if_t_out_z .* w.bnd(T).val(:,:,1) .* sz));

if(vol_out_2 > TINY)
  factor = vol_in / vol_out_2;
else
  factor = 1.0;  
end    

if verbatim == true
  disp( sprintf('+----------------------------+') );
  disp( sprintf('|  v_out_2   = %12.5e  |', vol_out_2  ) );
  disp( sprintf('|  factor    = %12.5e  |', factor   ) );
  disp( sprintf('+----------------------------+') );
end

%--------------------------------------
% Correction to satisfy volume balance
%--------------------------------------
u.bnd(W).val(1,:,:) = u.bnd(W).val(1,:,:) .* not(if_w_out_x) ...
                    + u.bnd(W).val(1,:,:) .*     if_w_out_x * factor;
u.bnd(E).val(1,:,:) = u.bnd(E).val(1,:,:) .* not(if_e_out_x) ...
                    + u.bnd(E).val(1,:,:) .*     if_e_out_x * factor;
v.bnd(S).val(:,1,:) = v.bnd(S).val(:,1,:) .* not(if_s_out_y) ...
                    + v.bnd(S).val(:,1,:) .*     if_s_out_y * factor;
v.bnd(N).val(:,1,:) = v.bnd(N).val(:,1,:) .* not(if_n_out_y) ...
                    + v.bnd(N).val(:,1,:) .*     if_n_out_y * factor;
w.bnd(B).val(:,:,1) = w.bnd(B).val(:,:,1) .* not(if_b_out_z) ...
                    + w.bnd(B).val(:,:,1) .*     if_b_out_z * factor;
w.bnd(T).val(:,:,1) = w.bnd(T).val(:,:,1) .* not(if_t_out_z) ...
                    + w.bnd(T).val(:,:,1) .*     if_t_out_z * factor;

end