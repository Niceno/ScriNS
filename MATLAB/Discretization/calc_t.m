%==========================================================================
function [t] = calc_t( t,                   ... 
                       fu, fv, fw,          ...
                       rho_cap, lambda,     ... 
                       dt, ts, dx, dy, dz, obst )
%--------------------------------------------------------------------------

% Variable to keep the old time step
persistent ts_old;  if isempty(ts_old)  ts_old = 0;  end

% If you entered here first time in a time step, set the old value
if ts ~= ts_old
  t.old  = t.val;
  ts_old = ts;
end

r_c = size( t.val ); n_c = prod( r_c );  % resolutions
                                
[A_t, b_t] = create_matrix(t,            ... % variable
                           rho_cap./dt,  ... % innertial term
                           lambda,       ... % diffusive coefficent
                           dx, dy, dz,   ... % dimensions
                           obst);            % obstacles
                    
% Advection terms for temperatures
c_t = advection(rho_cap, t, fu, fv, fw, dx, dy, dz, dt, 'smart');

% Innertial term for enthalpy
i_t = t.old .* rho_cap .* dx .* dy .* dz / dt;

% The entire source term
f_t = b_t - c_t + i_t;
  
% Solve temperature
if SOL == 'd'
  t.val = reshape( A_t\(reshape(f_t, n_c, 1)), r_c );
else  
  t.val = reshape( bicgstab(A_t, reshape(f_t, n_c, 1), TOL), r_c );
end  

t = adj_n_bnds(t);

end

