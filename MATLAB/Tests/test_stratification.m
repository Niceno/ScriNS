%==========================================================================
% Solves flow in a channel with stable or unstable stratification.  
% The type of stratification is set by parameter "STRATIFICATION"
%
% Uses non-Boussinesq model for buoyancy term; i.e. density depends on
% temperature.
%
% Gravity term is engaged gradually to avoid vortex at the outlet.
%--------------------------------------------------------------------------
clear;

%-------------------------------------------------
% Chose between stable or unstable stratification
%-------------------------------------------------
STRATIFICATION = 's'; % 's' or 'u'

%------------------------------
% Set path to common functions
%------------------------------
path(path, '../Constants');
path(path, '../Discretization');
path(path, '../InputOutput');
path(path, '../Operators');
path(path, '../PhysicalModels');

%==========================================================================
%
% Define problem
%
%==========================================================================

% Set domain dimensions and grid resolution
xn = nodes(0, 10,   80);  xc = avg(xn);
yn = nodes(0,  1,   20);  yc = avg(yn);
zn = nodes(0,  0.2,  5);  zc = avg(zn);   

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);


r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% Set physical properties
rho   (1:nx,1:ny,1:nz) = 1.25;    % density 
mu    (1:nx,1:ny,1:nz) = 0.1;     % viscosity
lambda(1:nx,1:ny,1:nz) = 0.001;   % thermal conductivity
cap   (1:nx,1:ny,1:nz) = 1.0;     % thermal capacity

% Time-stepping parameters
dt     =   0.20;    % time step
ndt    = 500;       % number of time steps

% Create unknowns; names, positions and sizes
u = create_unknown('u-velocity',  X, r.u, 'd');
v = create_unknown('v-velocity',  Y, r.v, 'd');
w = create_unknown('w-velocity',  Z, r.w, 'd');
p = create_unknown('pressure',    C, r.c, 'n');
t = create_unknown('temperature', C, r.c, 'n');
p_tot = zeros( r.c );

% Specify boundary conditions
u.bnd(W).type(1,:,:)='d'; 
for k=1:nz
  u.bnd(W).val(1,:,k)=avg(par(0.1, ny+1));    
end
u.bnd(E).type(1,:,:)='o'; 
u.bnd(E).val (1,:,:)=0.1;

u.bnd(B).type(:,:,1) = 'n'; 
u.bnd(T).type(:,:,1) = 'n'; 
v.bnd(B).type(:,:,1) = 'n'; 
v.bnd(T).type(:,:,1) = 'n'; 
w.bnd(B).type(:,:,1) = 'n'; 
w.bnd(T).type(:,:,1) = 'n'; 
t.bnd(B).type(:,:,1) = 'n'; 
t.bnd(T).type(:,:,1) = 'n'; 

t.val(1:nx, 1:ny, 1:nz) = 70;
if STRATIFICATION == 'u'
  dtemp = (60-80)/ny;
  t.bnd(W).type(1,1:ny,1:nz) = 'd'; 
  for k=1:nz
    t.bnd(W).val (1,1:ny,k) = linspace(60-dtemp/2,80+dtemp/2,ny);
  end
  t.bnd(S).type(1:nx,1,1:nz) = 'd'; 
  t.bnd(S).val (1:nx,1,1:nz) = 60.0;
  t.bnd(N).type(1:nx,1,1:nz) = 'd'; 
  t.bnd(N).val (1:nx,1,1:nz) = 80.0;
elseif STRATIFICATION == 's'
  dtemp = (80-60)/ny;
  t.bnd(W).type(1,1:ny) = 'd'; 
  for k=1:nz
    t.bnd(W).val(1,1:ny,k) = linspace(80-dtemp/2,60+dtemp/2,ny);
  end  
  t.bnd(S).type(1:nx,1,1:nz) = 'd'; 
  t.bnd(S).val (1:nx,1,1:nz) = 80.0;
  t.bnd(N).type(1:nx,1,1:nz) = 'd'; 
  t.bnd(N).val (1:nx,1,1:nz) = 60.0;
end
t.bnd(E).type(1,1:ny,1:nz) = 'n'; 
t.bnd(E).val (1,1:ny,1:nz) = 70.0;

t = adj_n_bnds(t); 
p = adj_n_bnds(p);

obst = [];

%==========================================================================
%
% Solution algorithm
%
%==========================================================================

%-----------
%
% Time loop 
%
%-----------
for ts=1:ndt

  disp_time_step(ts);
  
  rho = 1.06 + (t.val - 60.0) * (0.99 - 1.06) / 20.0;
  
  %------------------------
  % Temperature (enthalpy)
  %------------------------
  t = calc_t(t, u, v, w,                 ... % variable and advection velocity
             (rho.*cap), lambda,         ... % physical properties
             dt, ts, dx, dy, dz, obst);      % numerical domain
  
  %----------
  % Momentum
  %----------

  % Gravity term for momentum 
  g_v = - G * avg(Y, rho);
  
  [u, v, w] = calc_uvw(u, v, w, u, v, w,             ...
                       rho, mu,                      ...
                       p_tot,                        ...
                       zeros(r.u), g_v, zeros(r.w),  ...
                       dt, ts, dx, dy, dz, obst);

  %----------
  % Pressure
  %----------
  p = calc_p(p, u, v, w,         ... % variable and advection velocity
             rho,                ... % density
             dt, dx, dy, dz, obst);  % numerical domain

  p_tot = p_tot + p.val;
  
  %---------------------
  % Velocity correction
  %---------------------
  [u, v, w] = corr_uvw(u, v, w, p, rho, dt, dx, dy, dz, obst);

  % Compute the source for pressure again to check mass error
  err = vol_balance(u,v,w, dx,dy,dz, obst);
  disp( sprintf('maximum mass error after correction  = %12.5e', ...
        max(max(max(abs(err))))));
  
  cfl = cfl_max(u,v,w, dt, dx,dy,dz);
  disp( sprintf('cfl max = %5.3f', cfl ) );
 
%==========================================================================
%
% Visualisation
%
%==========================================================================
  
  %------
  % Plot
  %------
  if mod(ts,50) == 0
    k_plot = 3;
    uc = avg(X, cat(X, u.bnd(W).val,u.val,u.bnd(E).val));
    vc = avg(Y, cat(Y, v.bnd(S).val,v.val,v.bnd(N).val));
    subplot(3,1,1)
      [cnt,hnd] = contourf(xc, yc, t.val(:,:,k_plot)', ...
                           linspace(59.99,80.01,10)); 
      hold on;
      quiver(xc, yc, uc(:,:,k_plot)', vc(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title(sprintf('Temperature at time step %d out of %d\n', ts, ndt))
    subplot(3,1,2)
      [cnt,hnd] = contourf(xc,yc,rho(:,:,k_plot)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xc, yc, uc(:,:,k_plot)', vc(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title(sprintf('Density at time step %d out of %d\n', ts, ndt))
    subplot(3,1,3)
      [cnt,hnd] = contourf(xc,yc,p_tot(:,:,k_plot)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xc, yc, uc(:,:,k_plot)', vc(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title(sprintf('Total pressure at time step %d out of %d\n', ts, ndt))
    drawnow;

    %--------------------------
    % Export to Tecplot format
    %--------------------------
    container.var(1:5) = [u v w t p];
    file_name = sprintf('results-stratification-staggered-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);
          
  end
end

