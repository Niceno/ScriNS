%==========================================================================
function [u,v,w, uf,vf,wf] = ...
    calc_uvw(u,  v,  w,      ... % transported velocity (coll. or stag.)  
             uf, vf, wf,     ... % advection velocity (staggered)
             rho, mu,        ... % physical properties 
             p_tot,          ... % total pressure 
             e_x, e_y, e_z,  ... % external forces in x, y and z direction
             dt, ts, dx, dy, dz, obst )
%--------------------------------------------------------------------------

% Variable to keep the old time step
persistent ts_old;  if isempty(ts_old)  ts_old = 0;  end

% If you entered here first time in a time step, set the old value
if ts ~= ts_old
  u.old  = u.val;
  v.old  = v.val;
  w.old  = w.val;
  ts_old = ts;
end

r_u = size( u.val ); n_u = prod( r_u );  % resolutions
r_v = size( v.val ); n_v = prod( r_v );  % resolutions
r_w = size( w.val ); n_w = prod( r_w );  % resolutions

dv = dx .* dy .* dz;

pos = u.pos;

% Create linear systems
[A_u, b_u] = create_matrix(u,           ... % variable
                           rho./dt,     ... % inertial term
                           mu,          ... % diffusive coefficient
                           dx, dy, dz,  ... % dimensions
                           obst, 'd');      % obstacles
 
[A_v, b_v] = create_matrix(v,           ... % variable
                           rho./dt,     ... % inertial term
                           mu,          ... % diffusive coefficient
                           dx, dy, dz,  ... % dimensions
                           obst, 'd');      % obstacles
                           
[A_w, b_w] = create_matrix(w,           ... % variable
                           rho./dt,     ... % inertial term
                           mu,          ... % diffusive coefficient
                           dx, dy, dz,  ... % dimensions
                           obst, 'd');      % obstacles

% Advection terms for momentum                            
c_u = advection(rho, u, uf, vf, wf, dx, dy, dz, dt, 'superbee');
c_v = advection(rho, v, uf, vf, wf, dx, dy, dz, dt, 'superbee');
c_w = advection(rho, w, uf, vf, wf, dx, dy, dz, dt, 'superbee');

% Innertial term for momentum
i_u = u.old .* avg(u.pos, rho) .* avg(u.pos, dv) ./ dt;
i_v = v.old .* avg(v.pos, rho) .* avg(v.pos, dv) ./ dt;
i_w = w.old .* avg(w.pos, rho) .* avg(w.pos, dv) ./ dt;
  
% Compute staggered pressure gradients
p_tot_x = dif(X, p_tot) ./ avg(X, dx);  
p_tot_y = dif(Y, p_tot) ./ avg(Y, dy);  
p_tot_z = dif(Z, p_tot) ./ avg(Z, dz);

% Make pressure gradients cell-centered
if pos == C 
  p_tot_x = avg(X, cat(X, p_tot_x(1,:,:), p_tot_x, p_tot_x(end,:,:)));
  p_tot_y = avg(Y, cat(Y, p_tot_y(:,1,:), p_tot_y, p_tot_y(:,end,:)));
  p_tot_z = avg(Z, cat(Z, p_tot_z(:,:,1), p_tot_z, p_tot_z(:,:,end)));
end

% Total pressure gradients
p_st_u = p_tot_x .* avg(u.pos, dv);
p_st_v = p_tot_y .* avg(v.pos, dv);
p_st_w = p_tot_z .* avg(w.pos, dv);
  
% Full force terms for momentum equations
f_u = b_u - c_u + i_u - p_st_u + e_x .* avg(u.pos, dv);
f_v = b_v - c_v + i_v - p_st_v + e_y .* avg(v.pos, dv);
f_w = b_w - c_w + i_w - p_st_w + e_z .* avg(w.pos, dv); 

% Take care of obsts in the domian
if size(obst,1) ~= 0
  f_u = obst_zero_val(u.pos, f_u, obst);
  f_v = obst_zero_val(v.pos, f_v, obst);
  f_w = obst_zero_val(w.pos, f_w, obst);
end

% Solve velocities
if SOL == 'd'
  u.val = reshape( A_u\(reshape(f_u, n_u, 1)), r_u );
  v.val = reshape( A_v\(reshape(f_v, n_v, 1)), r_v );
  w.val = reshape( A_w\(reshape(f_w, n_w, 1)), r_w );
else
  u.val = reshape( bicgstab(A_u, reshape(f_u, n_u, 1), TOL), r_u );
  v.val = reshape( bicgstab(A_v, reshape(f_v, n_v, 1), TOL), r_v );
  w.val = reshape( bicgstab(A_w, reshape(f_w, n_w, 1), TOL), r_w );
end
  
% Update velocities in boundary cells
[u, v, w] = adj_o_bnds(u, v, w, dx, dy, dz, dt);

% Update face velocities (also substract cell-centered pressure gradients 
%                         and add staggered pressure gradients)
if pos == C
  uf.val = avg(X,u.val + dt ./       rho  .* (      p_tot_x     )) ...
                       - dt ./ avg(X,rho) .* (dif(X,p_tot) ./ avg(X,dx));  
  vf.val = avg(Y,v.val + dt ./       rho  .* (      p_tot_y     )) ...
                       - dt ./ avg(Y,rho) .* (dif(Y,p_tot) ./ avg(Y,dy));  
  wf.val = avg(Z,w.val + dt ./       rho  .* (      p_tot_z     )) ...
                       - dt ./ avg(Z,rho) .* (dif(Z,p_tot) ./ avg(Z,dz));  
  for j=1:6
    uf.bnd(j).val  = u.bnd(j).val;  
    vf.bnd(j).val  = v.bnd(j).val;  
    wf.bnd(j).val  = w.bnd(j).val;  
  end
else
  uf.val = u.val;  uf.bnd = u.bnd;
  vf.val = v.val;  vf.bnd = v.bnd;
  wf.val = w.val;  wf.bnd = w.bnd; 
end

if size(obst,1) ~= 0
  uf.val = obst_zero_val(X, uf.val, obst);
  vf.val = obst_zero_val(Y, vf.val, obst);
  wf.val = obst_zero_val(Z, wf.val, obst);
end
  
end