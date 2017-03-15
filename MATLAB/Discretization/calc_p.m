%==========================================================================
function [p] = calc_p( p,           ... 
                       uf, vf, wf,  ...
                       rho,         ... 
                       dt, dx, dy, dz, obst )
%--------------------------------------------------------------------------

r_c = size( p.val ); n_c = prod( r_c );  % resolutions

% Create linear system
[A_p, b_p] = create_matrix(p,             ... % variable
                           zeros( r_c ),  ... % inertial term
                           dt./rho,       ... % diffusive coefficient
                           dx, dy, dz,    ... % cell dimensions
                           obst, 'n');    ... % obstacles
            
% Compute the source for the pressure
% Important: don't send "obst" as a parameter here, because you don't want
% to take it into account at this stage.  After velocity corrections, yes!
b_p = vol_balance(uf, vf, wf, dx, dy, dz);
    
disp( sprintf('maximum volume error before correction = %12.5e', ...
      max(max(max(abs(b_p))))));
disp( sprintf('volume imbalance before correction     = %12.5e', ...
      sum(sum(sum(b_p)))));

% Solve pressure
if SOL == 'd'
  p.val = reshape( A_p\(reshape(b_p, n_c, 1)), r_c );
else
  L = ichol(A_p);
  p.val = reshape( bicgstab(A_p, reshape(b_p, n_c, 1), ... % A, b
                            TOL,                       ... % tolerance 
                            n_c,                       ... % max iter
                            L, L',                     ... % preconditioner
                            reshape(p.val, n_c, 1)),   ... % initial value
                            r_c);
end               

% Anchor it to values around zero (the absolute value of pressure
% correction can get really volatile.  Although it is in prinicple not
% important for incompressible flows, it is ugly for post-processing.
p.val = p.val - mean(mean(mean(p.val)));

% Set to zero in obstacle (it can get strange 
% values during the iterative solution procedure)
if size(obst,1) ~= 0
  p.val = obst_zero_val(C, p.val, obst);
end

% Finally adjust the boundary values
p = adj_n_bnds(p);

end

