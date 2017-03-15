%==========================================================================
function [A, b] = create_matrix(phi, in, mu, dx, dy, dz, obst, obc)
%--------------------------------------------------------------------------
% pos   - position of variable (C - central, 
%                               X - staggered in x direction,
%                               Y - staggered in y direction,
%                               Z - staggered in z direction)
% in         - innertial term
% mu         - viscous coefficient
% dx, dy, dz - cell size in x, y and z directions
% obc        - obstacles's boundary condition, ('n' - Neumann,
%                                               'd' - Dirichlet)
%--------------------------------------------------------------------------

%-------------------------------
% Create right hand side vector
%-------------------------------
[nx ny nz] = size(phi.val);  

b = zeros(nx, ny, nz);

%------------------------------------
% Create default matrix coefficients 
%------------------------------------
D = phi.pos;

% Handle central coefficient due to innertia
c.P = avg(D, in) .* avg(D, dx.*dy.*dz);

% Pre-compute geometrical quantities
sx = dy .* dz;
sy = dx .* dz;
sz = dx .* dy;

% Compute default matrix coeffiecient values  
c.W = cat(X,                                                            ...
    avg(D,mu(  1,:,:)) .* avg(D,sx(  1,:,:)) ./ avg(D,(dx(  1,:,:))/2), ...
    avg(D,avg(X, mu))  .* avg(D,avg(X, sx))  ./ avg(D,avg(X, dx))); 
c.E = cat(X,                                                            ...
    avg(D,avg(X, mu))  .* avg(D,avg(X, sx))  ./ avg(D,avg(X, dx)),      ...
    avg(D,mu(end,:,:)) .* avg(D,sx(end,:,:)) ./ avg(D,(dx(end,:,:)/2)));
c.S = cat(Y,                                                            ...
    avg(D,mu(:,  1,:)) .* avg(D,sy(:,  1,:)) ./ avg(D,(dy(:,  1,:)/2)), ...
    avg(D,avg(Y, mu))  .* avg(D,avg(Y, sy))  ./ avg(D,avg(Y, dy))); 
c.N = cat(Y,                                                            ...
    avg(D,avg(Y, mu))  .* avg(D,avg(Y, sy))  ./ avg(D,avg(Y, dy)),      ...
    avg(D,mu(:,end,:)) .* avg(D,sy(:,end,:)) ./ avg(D,(dy(:,end,:)/2)));
c.B = cat(Z,                                                            ...
    avg(D,mu(:,:,  1)) .* avg(D,sz(:,:,  1)) ./ avg(D,(dz(:,:,  1)/2)), ...
    avg(D,avg(Z, mu))  .* avg(D,avg(Z, sz))  ./ avg(D,avg(Z, dz))); 
c.T = cat(Z,                                                            ...
    avg(D,avg(Z, mu))  .* avg(D,avg(Z, sz))  ./ avg(D,avg(Z, dz)),      ...
    avg(D,mu(:,:,end)) .* avg(D,sz(:,:,end)) ./ avg(D,(dz(:,:,end)/2))); 

% Correct for staggered variables
if     D == X
  c.W = mu(1:end-1,:,:) .* sx(1:end-1,:,:) ./ dx(1:end-1,:,:);
  c.E = mu(2:end,  :,:) .* sx(2:end,  :,:) ./ dx(2:end,  :,:);
elseif D == Y
  c.S = mu(:,1:end-1,:) .* sy(:,1:end-1,:) ./ dy(:,1:end-1,:);
  c.N = mu(:,2:end,  :) .* sy(:,2:end,  :) ./ dy(:,2:end,  :);
elseif D == Z
  c.B = mu(:,:,1:end-1) .* sz(:,:,1:end-1) ./ dz(:,:,1:end-1);
  c.T = mu(:,:,2:end  ) .* sz(:,:,2:end  ) ./ dz(:,:,2:end  );
end    

%-----------------------------------------------------------------------
% Zero them (correct them) for vanishing derivative boundary condition.
%-----------------------------------------------------------------------

% The values defined here will be false (0) wherever there is no inlet.
if_w_d = ( phi.bnd(W).type(1,:,:) == 'd' );
if_e_d = ( phi.bnd(E).type(1,:,:) == 'd' );
if_s_d = ( phi.bnd(S).type(:,1,:) == 'd' );
if_n_d = ( phi.bnd(N).type(:,1,:) == 'd' );
if_b_d = ( phi.bnd(B).type(:,:,1) == 'd' );
if_t_d = ( phi.bnd(T).type(:,:,1) == 'd' );

c.W(  1,  :,  :) = c.W(  1,  :,  :) .* if_w_d;
c.E(end,  :,  :) = c.E(end,  :,  :) .* if_e_d;
c.S(  :,  1,  :) = c.S(  :,  1,  :) .* if_s_d;
c.N(  :,end,  :) = c.N(  :,end,  :) .* if_n_d;
c.B(  :,  :,  1) = c.B(  :,  :,  1) .* if_b_d;
c.T(  :,  :,end) = c.T(  :,  :,end) .* if_t_d;

%--------------------------------------------
% Fill the source terms with boundary values
%--------------------------------------------
b(  1,  :,  :) = b(  1,  :,  :) + c.W(  1,  :,  :).*phi.bnd(W).val(1,:,:);
b(end,  :,  :) = b(end,  :,  :) + c.E(end,  :,  :).*phi.bnd(E).val(1,:,:);
b(  :,  1,  :) = b(  :,  1,  :) + c.S(  :,  1,  :).*phi.bnd(S).val(:,1,:);
b(  :,end,  :) = b(  :,end,  :) + c.N(  :,end,  :).*phi.bnd(N).val(:,1,:);
b(  :,  :,  1) = b(  :,  :,  1) + c.B(  :,  :,  1).*phi.bnd(B).val(:,:,1);
b(  :,  :,end) = b(  :,  :,end) + c.T(  :,  :,end).*phi.bnd(T).val(:,:,1);

%---------------------------------------
% Correct system matrices for obstacles
%---------------------------------------
if size(obst,1) ~= 0
  c = obst_mod_matrix(phi, c, obst, obc);
end

%-----------------------------------------------
% Add all neighbours to the central matrix,
% and zero the coefficients towards boundaries
%-----------------------------------------------
c.P = c.P + c.W + c.E + c.S + c.N + c.B + c.T;

c.W(  1,  :,  :) = 0;  c.E(end,  :,  :) = 0;
c.S(  :,  1,  :) = 0;  c.N(  :,end,  :) = 0;
c.B(  :,  :,  1) = 0;  c.T(  :,  :,end) = 0;

%----------------------
% Create sparse matrix 
%----------------------
n = nx * ny * nz;

c.P = reshape(c.P, n, 1);
c.W = reshape(c.W, n, 1);  c.E = reshape(c.E, n, 1);
c.S = reshape(c.S, n, 1);  c.N = reshape(c.N, n, 1);
c.B = reshape(c.B, n, 1);  c.T = reshape(c.T, n, 1);

A = spdiags( [-c.T    -c.N    -c.E    c.P    -c.W    -c.S    -c.B  ], ...
             [-nx*ny  -nx     -1      0      +1      +nx     +nx*ny], ...
             n, n);
         
end