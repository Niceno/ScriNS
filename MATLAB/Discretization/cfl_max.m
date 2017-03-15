%==========================================================================
function [cfl] = cfl_max( u, v, w, dt, dx, dy, dz )
%--------------------------------------------------------------------------
% Computes the maximum CFL number in a computational domain.
%
% u, v, w    - velocity components, either collocated or staggered
% dt         - time step
% dx, dy, dz - cell dimensions
%--------------------------------------------------------------------------

% Take velocity's position
pos = u.pos;

% Mesh is cell-centered
if pos == C
  cfl = dt * max([ max(max(max(abs(u.val)./dx))), 
                   max(max(max(abs(v.val)./dy))),
                   max(max(max(abs(w.val)./dz))) ]);

% Mesh is staggered
else
  cfl = dt * max([ max(max(max(abs(u.val)./avg(X,dx)))), 
                   max(max(max(abs(v.val)./avg(Y,dy)))),
                   max(max(max(abs(w.val)./avg(Z,dz)))) ]);
end

end