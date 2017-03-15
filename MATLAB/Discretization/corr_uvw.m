%==========================================================================
function [u,v,w] = corr_uvw(u,v,w, p, rho, dt, dx,dy,dz, obst)
%--------------------------------------------------------------------------

% Compute pressure correction gradients
p_x = dif(X, p.val) ./ avg(X, dx);  
p_y = dif(Y, p.val) ./ avg(Y, dy);  
p_z = dif(Z, p.val) ./ avg(Z, dz);  

% Set to zero in obst
if size(obst,1) ~= 0
  p_x = obst_zero_val(X, p_x, obst);
  p_y = obst_zero_val(Y, p_y, obst);
  p_z = obst_zero_val(Z, p_z, obst);
end

% Pad with boundary values by expanding from interior 
% (This is done only for collocated formulation)
if u.pos == C
  p_x = avg(X, cat(X, p_x(1,:,:), p_x, p_x(end,:,:)));
  p_y = avg(Y, cat(Y, p_y(:,1,:), p_y, p_y(:,end,:)));
  p_z = avg(Z, cat(Z, p_z(:,:,1), p_z, p_z(:,:,end)));
end

% Correct the velocities
u.val = u.val - dt ./ avg(u.pos, rho) .* p_x;
v.val = v.val - dt ./ avg(v.pos, rho) .* p_y;
w.val = w.val - dt ./ avg(w.pos, rho) .* p_z;
  
end