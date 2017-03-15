%==========================================================================
%                                                       o ... scalars
%                          (n)                          - ... u velocities      
%                                                       | ... v velocities
%       +-------+-------+-------+-------+-------+          
%       |       |       |       |       |       |            
%       |   o   -   o   -   o   -   o   -   o   | j=ny     
%       |       |       |       |       |       | 
%       +---|---+---|---+---|---+---|---+---|---+     j=nym     
%       |       |       |       |       |       |
%       |   o   -   o   -   o   -   o   -   o   | ...
%       |       |       |       |       |       |
%  (w)  +---|---+---|---+---|---+---|---+---|---+    j=2        (e)
%       |       |       |       |       |       |
%       |   o   -   o   -   o   -   o   -   o   | j=2
%       |       |       |       |       |       |
%       +---|---+---|---+---|---+---|---+---|---+    j=1 (v-velocity)
%       |       |       |       |       |       |
%       |   o   -   o   -   o   -   o   -   o   | j=1   (scalar cell)
%       |       |       |       |       |       |
%       +-------+-------+-------+-------+-------+
%  y       i=1     i=2     ...     ...     i=nx      (scalar cells)
% ^            i=1      i=2    ...    i=nxm      (u-velocity cells)
% |
% +---> x                  (s)
%  
%--------------------------------------------------------------------------
clear;

%------------------------------
% set path to common functions
%------------------------------
path(path, '../Constants');
path(path, '../Discretization');
path(path, '../InputOutput');
path(path, '../Operators');
path(path, '../PhysicalModels');

%----------------
% define problem
%----------------

% Set domain dimensions and grid resolution
xn = nodes(0, 1.25,  256);
yn = nodes(0, 0.125,  32);
zn = nodes(0, 0.125,  32);

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% set physical properties
[rho, mu, cap, lambda] = air_properties(1:nx, 1:ny, 1:nz);

% time-stepping parameters
dt     =    0.005;    % time steps
ndt    = 2;       % number of time steps

% create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('u-velocity',  X, r.u, 'd');
vf = create_unknown('v-velocity',  Y, r.v, 'd');
wf = create_unknown('w-velocity',  Z, r.w, 'd');
p  = create_unknown('pressure',    C, r.c, 'n');
p_tot = zeros(r.c);

% specify boundary conditions
for k=1:nz
  uf.bnd(W).val(1,:,k)=avg(par(0.1, ny+1));    
end
uf.bnd(E).type(1,:,:) = 'o'; 

uf.bnd(B).type(:,:,1) = 'n'; 
uf.bnd(T).type(:,:,1) = 'n'; 
vf.bnd(B).type(:,:,1) = 'n'; 
vf.bnd(T).type(:,:,1) = 'n'; 
wf.bnd(B).type(:,:,1) = 'n'; 
wf.bnd(T).type(:,:,1) = 'n'; 

p = adj_n_bnds(p);

% create obstacles
th = 5;
block(1).im  =  3*nx/16;           % i minus
block(1).ip  =  block(1).im + th;  % i plus
block(1).jm  =  1;                 % j minus
block(1).jp  =  3*ny/4;            % j plus
block(1).km  =  1;                 % k minus
block(1).kp  =  3*ny/4;            % k plus

block(2).im  =  5*nx/16;           % i minus
block(2).ip  =  block(2).im + th;  % i plus
block(2).jm  =  ny/4 + 1;          % j minus
block(2).jp  =  ny;                % j plus
block(2).km  =  ny/4 + 1;          % k minus
block(2).kp  =  ny;                % k plus
  
block(3).im  =  7*nx/16;           % i minus
block(3).ip  =  block(3).im + th;  % i plus
block(3).jm  =  1;                 % j minus
block(3).jp  =  3*ny/4;            % j plus
block(3).km  =  1;                 % k minus
block(3).kp  =  3*ny/4;            % k plus

block(4).im  =  9*nx/16;           % i minus
block(4).ip  =  block(4).im + th;  % i plus
block(4).jm  =  ny/4 + 1;          % j minus
block(4).jp  =  ny;                % j plus
block(4).km  =  ny/4 + 1;          % k minus
block(4).kp  =  ny;                % k plus

obst = zeros(nx,ny,nz);
for o=1:4
  for i=block(o).im:block(o).ip  
    for j=block(o).jm:block(o).jp
      for k=block(o).km:block(o).kp
        obst(i,j,k) = 1;
      end  
    end
  end
end

%-----------
%
% Time loop 
%
%-----------
for ts = 1:ndt

  disp_time_step(ts);
  
  %----------
  % Momentum
  %----------
  [uf,vf,wf, uf,vf,wf] = calc_uvw(uf,vf,wf, uf,vf,wf,                  ...
                                  rho, mu,                             ...
                                  p_tot,                               ...
                                  zeros(r.u), zeros(r.v), zeros(r.w),  ...
                                  dt, ts, dx,dy,dz, obst);

  %----------
  % Pressure
  %----------
  p = calc_p(p, uf,vf,wf,          ... % variable and advection velocity
             rho,                  ... % density
             dt, dx,dy,dz, obst);      % numerical domain
         
  p_tot = p_tot + p.val;
  
  %---------------------
  % Velocity correction
  %---------------------
  [uf,vf,wf] = corr_uvw(uf,vf,wf, p, rho, dt, dx,dy,dz, obst);

  % Compute the source for pressure again to check mass error
  err = vol_balance(uf,vf,wf, dx,dy,dz, obst);
  disp( sprintf('maximum mass error after correction  = %12.5e', ...
        max(max(max(abs(err))))));
  
  cfl = cfl_max( uf,vf,wf, dt, dx,dy,dz );
  disp( sprintf('cfl max = %5.3f', cfl ) );

  %------
  % Plot
  %------
  if mod(ts,10) == 0
    k_plot = 16;
    xp = avg(xn);  
    yp = avg(yn); 
    up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
    vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));

    subplot(3,1,1)
      [cnt,hnd] = contourf(xp,yp,p.val(:,:,k_plot)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,k_plot)', vp(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    subplot(3,1,2)
      [cnt,hnd] = contourf(xp,yp,p_tot(:,:,k_plot)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,k_plot)', vp(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Total pressure at time step %d ', ts),   ...
             sprintf('out of %d\n',                     ndt),  ...
             sprintf('CFL = %5.3f\n',                   cfl)]);
    subplot(3,1,3)
      [cnt,hnd] = contourf(xp,yp,sqrt([up(:,:,k_plot)'.*up(:,:,k_plot)' ...
                                     + vp(:,:,k_plot)'.*vp(:,:,k_plot)'])); 
      hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,k_plot)', vp(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Velocity magnitude at time step %d ', ts),   ...
             sprintf('out of %d\n',                         ndt),  ...
             sprintf('CFL = %5.3f\n',                       cfl)]);
    drawnow;
  end
  
  if mod(ts,100) == 0
    %--------------------------
    % Export to Tecplot format
    %--------------------------
    container.var(1:4) = [uf vf wf p];
    file_name = sprintf('results-obstacles-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);
  end
  
end  

