%==========================================================================
function [c] = advection(rho, phi, uf, vf, wf, dx, dy, dz, dt, lim_name)
%--------------------------------------------------------------------------
% rho        - advection coefficient
% phi        - transported variable
% uf, vf, wf - staggered advection velocity components
% dx, dy, dz - cell size in x, y and z direction
% lim_name   - limiter name
%--------------------------------------------------------------------------

[nx ny nz] = size(phi.val);

pos = phi.pos;

sx = dy .* dz;
sy = dx .* dz;
sz = dx .* dy;

%-------------------------------------------------
% Specific for cell-centered transported variable
%-------------------------------------------------
if pos==C  

  % Facial values of physical properties including boundary cells
  rho_x_fac = cat(X, rho(1,:,:), avg(X, rho), rho(nx,:,:));  % nxp,ny, nz
  rho_y_fac = cat(Y, rho(:,1,:), avg(Y, rho), rho(:,ny,:));  % nx, nyp,nz 
  rho_z_fac = cat(Z, rho(:,:,1), avg(Z, rho), rho(:,:,nz));  % nx, ny, nzp 

  % Facial values of areas including boundary cells
  a_x_fac = cat(X, sx(1,:,:), avg(X, sx), sx(nx,:,:));
  a_y_fac = cat(Y, sy(:,1,:), avg(Y, sy), sy(:,ny,:));
  a_z_fac = cat(Z, sz(:,:,1), avg(Z, sz), sz(:,:,nz));
  
  del_x = avg(X, dx);
  del_y = avg(Y, dy);
  del_z = avg(Z, dz);
  
  % Facial values of velocities without boundary values
  u_fac = uf.val;  % nxm,ny, nz
  v_fac = vf.val;  % nx, nym,nz
  w_fac = wf.val;  % nx, ny, nzm

  % Boundary velocity values
  u_bnd_W = uf.bnd(W).val;  u_bnd_E = uf.bnd(E).val;
  v_bnd_S = vf.bnd(S).val;  v_bnd_N = vf.bnd(N).val;
  w_bnd_B = wf.bnd(B).val;  w_bnd_T = wf.bnd(T).val;
end

%------------------------------------------------------------
% Specific for transported variable staggered in x direction
%------------------------------------------------------------
if pos==X  
    
  % Facial values of physical properties including boundary cells
  rho_x_fac = rho;                             % nx, ny, nz
  rho_nod_y = avg(X, avg(Y, rho) );            % nxm,nym,nz
  rho_y_fac = cat(Y, rho_nod_y(:,  1,:), ...   % nxm,nyp,nz
                     rho_nod_y(:,  :,:), ...
                     rho_nod_y(:,end,:));                      
  rho_nod_z = avg(X, avg(Z, rho) );            % nxm,ny,nzm
  rho_z_fac = cat(Z, rho_nod_z(:,:,  1), ...   % nxm,ny,nzp
                     rho_nod_z(:,:,  :), ...
                     rho_nod_z(:,:,end));                      

  % Facial values of areas including boundary cells
  a_x_fac = sx;
  a_y_fac = cat(Y, avg(X,sy(:,1,:)),avg(X, avg(Y, sy)),avg(X, sy(:,ny,:)));
  a_z_fac = cat(Z, avg(X,sz(:,:,1)),avg(X, avg(Z, sz)),avg(X, sz(:,:,nz)));

  del_x = dx(2:end-1,:,:);
  del_y = avg(X, avg(Y, dy));
  del_z = avg(X, avg(Z, dz));
  
  % Facial values of velocities without boundary values
  u_fac = avg(X, uf.val);  % nxmm,ny, nz
  v_fac = avg(X, vf.val);  % nxm, nym,nz
  w_fac = avg(X, wf.val);  % nxm, ny, nzm

  % Boundary velocity values
  u_bnd_W = uf.bnd(W).val;          u_bnd_E = uf.bnd(E).val;
  v_bnd_S = avg(X, vf.bnd(S).val);  v_bnd_N = avg(X, vf.bnd(N).val);
  w_bnd_B = avg(X, wf.bnd(B).val);  w_bnd_T = avg(X, wf.bnd(T).val);
end

%------------------------------------------------------------
% Specific for transported variable staggered in y direction
%------------------------------------------------------------
if pos==Y  
    
  % Facial values of physical properties including boundary cells
  rho_nod_x = avg(Y, avg(X, rho) );           % nxm,nym,nz
  rho_x_fac = cat(X, rho_nod_x(  1,:,:), ...  % nxp,nym,nz
                     rho_nod_x(  :,:,:), ...
                     rho_nod_x(end,:,:));   
  rho_y_fac = rho;                            % nx, ny, nz
  rho_nod_z = avg(Y, avg(Z, rho) );           % nx, nym,nzm
  rho_z_fac = cat(Z, rho_nod_z(:,:,  1), ...  % nx, nym,nzp
                     rho_nod_z(:,:,  :), ...
                     rho_nod_z(:,:,end));   

  % Facial values of areas including boundary cells
  a_x_fac = cat(X, avg(Y, sx(1,:,:)),avg(Y, avg(X,sx)),avg(Y, sx(nx,:,:)));
  a_y_fac = sy;
  a_z_fac = cat(Z, avg(Y, sz(:,:,1)),avg(Y, avg(Z,sz)),avg(Y, sz(:,:,nz)));

  del_x = avg(Y, avg(X, dx));
  del_y = dy(:,2:end-1,:);
  del_z = avg(Y, avg(Z, dz));
    
  % Facial values of velocities without boundary values
  u_fac = avg(Y, uf.val);  % nxm,nym, nz
  v_fac = avg(Y, vf.val);  % nx, nymm,nz
  w_fac = avg(Y, wf.val);  % nx, nym, nzm

  % Facial values of velocities with boundary values
  u_bnd_W = avg(Y, uf.bnd(W).val);  u_bnd_E = avg(Y, uf.bnd(E).val);
  v_bnd_S = vf.bnd(S).val;          v_bnd_N = vf.bnd(N).val;
  w_bnd_B = avg(Y, wf.bnd(B).val);  w_bnd_T = avg(Y, wf.bnd(T).val);
end

%------------------------------------------------------------
% Specific for transported variable staggered in z direction
%------------------------------------------------------------
if pos==Z  
    
  % Facial values of physical properties including boundary cells
  rho_nod_x = avg(Z, avg(X, rho) );           % nxm,ny, nzm
  rho_x_fac = cat(X, rho_nod_x(  1,:,:), ...  % nxp,ny, nzm
                     rho_nod_x(  :,:,:), ...
                     rho_nod_x(end,:,:));   
  rho_nod_y = avg(Z, avg(Y, rho) );           % nx, nym,nzm
  rho_y_fac = cat(Y, rho_nod_y(:,  1,:), ...  % nx, nyp,nzm
                     rho_nod_y(:,  :,:), ...
                     rho_nod_y(:,end,:));   
  rho_z_fac = rho;                            % nx, ny, nz

  % Facial values of areas including boundary cells
  a_x_fac = cat(X, avg(Z, sx(1,:,:)),avg(Z, avg(X,sx)),avg(Z, sx(nx,:,:)));
  a_y_fac = cat(Y, avg(Z, sy(:,1,:)),avg(Z, avg(Y,sy)),avg(Z, sy(:,ny,:)));
  a_z_fac = sz;

  del_x = avg(Z, avg(X,dx));
  del_y = avg(Z, avg(Y,dy));
  del_z = dz(:,:,2:end-1);
    
  % Facial values of velocities without boundary values
  u_fac = avg(Z, uf.val);  % nxm,ny,  nzm
  v_fac = avg(Z, vf.val);  % nx, nym, nzm
  w_fac = avg(Z, wf.val);  % nx, ny,  nzmm

  % Facial values of velocities with boundary values
  u_bnd_W = avg(Z, uf.bnd(W).val);  u_bnd_E = avg(Z, uf.bnd(E).val);
  v_bnd_S = avg(Z, vf.bnd(S).val);  v_bnd_N = avg(Z, vf.bnd(N).val);
  w_bnd_B = wf.bnd(B).val;          w_bnd_T = wf.bnd(T).val;
end

%------------------------------
% Common part of the algorithm
%------------------------------

%------------------------------------------------------------
%
%    |-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|-o-|
%      1   2   3   4   5   6   7   8   9   10     phi
%        x---x---x---x---x---x---x---x---x      
%        1   2   3   4   5   6   7   8   9        d_x initial
%    0---x---x---x---x---x---x---x---x---x---0      
%    1   2   3   4   5   6   7   8   9  10  11    d_x padded
%
%------------------------------------------------------------

% Compute consecutive differences (and avoid division by zero)
d_x = dif(X, phi.val);  % nxm, ny, nz  
d_x((d_x >  -TINY) & (d_x <=   0.0)) = -TINY; 
d_x((d_x >=   0.0) & (d_x <  +TINY)) = +TINY; 
d_x = cat(X, d_x(1,:,:), d_x, d_x(end,:,:));

d_y = dif(Y, phi.val);  % nx, nym, nz  
d_y((d_y >  -TINY) & (d_y <=   0.0)) = -TINY; 
d_y((d_y >=   0.0) & (d_y <  +TINY)) = +TINY; 
d_y = cat(Y, d_y(:,1,:), d_y, d_y(:,end,:));
  
d_z = dif(Z, phi.val);  % nx, ny, nzm  
d_z((d_z >  -TINY) & (d_z <   0.0)) = -TINY; 
d_z((d_z >=   0.0) & (d_z < +TINY)) = +TINY; 
d_z = cat(Z, d_z(:,:,1), d_z, d_z(:,:,end));
  
% Ratio of consecutive gradients for positive and negative flow
r_x_we = d_x(2:end-1,:,:) ./ d_x(1:end-2,:,:);  % nxm,ny, nz
r_x_ew = d_x(3:end,  :,:) ./ d_x(2:end-1,:,:);  % nxm,ny, nz

r_y_sn = d_y(:,2:end-1,:) ./ d_y(:,1:end-2,:);  % nx, nym,nz
r_y_ns = d_y(:,3:end,  :) ./ d_y(:,2:end-1,:);  % nx, nym,nz

r_z_bt = d_z(:,:,2:end-1) ./ d_z(:,:,1:end-2);  % nx, ny, nzm
r_z_tb = d_z(:,:,3:end  ) ./ d_z(:,:,2:end-1);  % nx, ny, nzm

flow_we = u_fac >= 0;    flow_ew = not(flow_we);
flow_sn = v_fac >= 0;    flow_ns = not(flow_sn);
flow_bt = w_fac >= 0;    flow_tb = not(flow_bt);

r_x = r_x_we .* flow_we + r_x_ew .* flow_ew;
r_y = r_y_sn .* flow_sn + r_y_ns .* flow_ns;
r_z = r_z_bt .* flow_bt + r_z_tb .* flow_tb;

% Apply a limiter
if( strcmp(lim_name, 'upwind') )
  psi_x = r_x * 0.0;                                         % upwind
  psi_y = r_y * 0.0;                                         % upwind
  psi_z = r_z * 0.0;                                         % upwind
elseif( strcmp(lim_name, 'minmod') )
  psi_x = max     (0.0, min(r_x, 1.0));                      % minmod
  psi_y = max     (0.0, min(r_y, 1.0));                      % minmod
  psi_z = max     (0.0, min(r_z, 1.0));                      % minmod
elseif( strcmp(lim_name, 'superbee') )
  psi_x = max_of_3(0.0, min(2.0*r_x, 1.0), min(r_x, 2.0));   % superbee
  psi_y = max_of_3(0.0, min(2.0*r_y, 1.0), min(r_y, 2.0));   % superbee
  psi_z = max_of_3(0.0, min(2.0*r_z, 1.0), min(r_z, 2.0));   % superbee
elseif( strcmp(lim_name, 'smart') )
  psi_x = max(0.0, min_of_3(2.0*r_x, (0.25+0.75*r_x), 4.0)); % smart
  psi_y = max(0.0, min_of_3(2.0*r_y, (0.25+0.75*r_y), 4.0)); % smart
  psi_z = max(0.0, min_of_3(2.0*r_z, (0.25+0.75*r_z), 4.0)); % smart
elseif( strcmp(lim_name, 'koren') )
  psi_x = max(0.0, min_of_3(2.0*r_x, (2.0+r_x)/3.0, 2.0));   % koren
  psi_y = max(0.0, min_of_3(2.0*r_y, (2.0+r_y)/3.0, 2.0));   % koren
  psi_z = max(0.0, min_of_3(2.0*r_z, (2.0+r_z)/3.0, 2.0));   % koren
end

flux_fac_lim_x =   phi.val(1:end-1,:,:) .* u_fac .* flow_we           ...
               +   phi.val(2:end,  :,:) .* u_fac .* flow_ew           ...
               +   0.5 * abs(u_fac) .* (1 - abs(u_fac) * dt ./ del_x) ...
               .*  (   psi_x(:,:,:) .* d_x(1:nx-1,:,:) .* flow_we     ...
                     + psi_x(:,:,:) .* d_x(2:nx,  :,:) .* flow_ew );  

flux_fac_lim_y =   phi.val(:,1:end-1,:) .* v_fac .* flow_sn           ...
               +   phi.val(:,2:end,  :) .* v_fac .* flow_ns           ...
               +   0.5 * abs(v_fac) .* (1 - abs(v_fac) * dt ./ del_y) ...
               .*  (   psi_y(:,:,:) .* d_y(:,1:ny-1,:) .* flow_sn     ...
                     + psi_y(:,:,:) .* d_y(:,2:ny,  :) .* flow_ns );

flux_fac_lim_z =   phi.val(:,:,1:end-1) .* w_fac .* flow_bt           ...
               +   phi.val(:,:,2:end)   .* w_fac .* flow_tb           ...
               +   0.5 * abs(w_fac) .* (1 - abs(w_fac) * dt ./ del_z) ...
               .*  (   psi_z(:,:,:) .* d_z(:,:,1:nz-1) .* flow_bt     ...
                     + psi_z(:,:,:) .* d_z(:,:,2:nz)   .* flow_tb );

% Pad with boundary values
flux_fac_lim_x = cat(X, phi.bnd(W).val .* u_bnd_W, ...
                        flux_fac_lim_x,            ...
                        phi.bnd(E).val .* u_bnd_E);
flux_fac_lim_y = cat(Y, phi.bnd(S).val .* v_bnd_S, ...
                        flux_fac_lim_y,            ...
                        phi.bnd(N).val .* v_bnd_N);
flux_fac_lim_z = cat(Z, phi.bnd(B).val .* w_bnd_B, ...
                        flux_fac_lim_z,            ...
                        phi.bnd(T).val .* w_bnd_T);

% Multiply with face areas                     
flux_fac_lim_x = rho_x_fac .* flux_fac_lim_x .* a_x_fac;
flux_fac_lim_y = rho_y_fac .* flux_fac_lim_y .* a_y_fac;
flux_fac_lim_z = rho_z_fac .* flux_fac_lim_z .* a_z_fac;
  
% Sum contributions from all directions up
c = dif(X, flux_fac_lim_x) + ...
    dif(Y, flux_fac_lim_y) + ...
    dif(Z, flux_fac_lim_z);

end
 
