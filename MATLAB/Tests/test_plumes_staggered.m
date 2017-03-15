%==========================================================================
% This script is to reproduce two-dimensional mixed convection case, with
% the aim of testing the outflow boundary, particularly the "convective"
% boundary condition which allows eddies to leave the domain.
%
% There is also a collocated version of this script, called
% "demo_plums_collocated.m".  It would be good to keep both version as
% similar as possible to each other, to test the differences between
% staggered and collocated arrangements always possible.
%--------------------------------------------------------------------------
clear
clc;

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
xn = nodes(0, 10, 300);                xc = avg(xn);
yn = nodes(0,  1,  40, 1/500, 1/500);  yc = avg(yn);
zn = nodes(0,  3,   3);                zc = avg(zn);

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% Set physical properties
prop.rho   (1:nx,1:ny,1:nz) = 1;     % density 
prop.mu    (1:nx,1:ny,1:nz) = 0.1;   % viscosity
prop.lambda(1:nx,1:ny,1:nz) = 0.15;  % thermal conductivity
prop.cap   (1:nx,1:ny,1:nz) = 1.0;   % thermal capacity

% Time-stepping parameters
dt     =    0.003;  % time step
ndt    = 150;      % number of time steps

% Create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('u-velocity',  X, r.u, 'd');
vf = create_unknown('v-velocity',  Y, r.v, 'd');
wf = create_unknown('w-velocity',  Z, r.w, 'd');
p  = create_unknown('pressure',    C, r.c, 'n');
t  = create_unknown('temperature', C, r.c, 'n');

% Specify boundary conditions
uf.bnd(W).type(1,:,:) = 'd'; 
for k=1:nz
  uf.bnd(W).val(1,:,k)  = avg(par(1.0,ny+1));
end  
uf.bnd(E).type(1,:,:) = 'o'; uf.bnd(E).val(1,:,:)  = 1.0;

uf.bnd(B).type(:,:,1) = 'n'; 
uf.bnd(T).type(:,:,1) = 'n'; 
vf.bnd(B).type(:,:,1) = 'n'; 
vf.bnd(T).type(:,:,1) = 'n'; 
wf.bnd(B).type(:,:,1) = 'n'; 
wf.bnd(T).type(:,:,1) = 'n'; 

t.bnd(W).type(1,:,:) = 'd';
for k=1:nz
  t.bnd(W).val(1,:,k) = 1.0-yc;
end  
t.bnd(S).type(:,1,:) = 'd'; t.bnd(S).val(:,1,:) = +1.0;
t.bnd(N).type(:,1,:) = 'd'; t.bnd(N).val(:,1,:) =  0.0;

t = adj_n_bnds(t);
p = adj_n_bnds(p);

% Specify initial conditions
uf.val(:,:,:) = 1.0;
t.val(:,:,:) = 0;

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

  disp_time_step( ts );

  %------------------------
  % Temperature (enthalpy)
  %------------------------
  t = calc_t(t, uf,vf,wf,                        ... % advection velocity
             (prop.rho.*prop.cap), prop.lambda,  ... % physical properties
             dt, ts, dx,dy,dz, obst);                % numerical domain

  %----------
  % Momentum
  %----------

  % Gravity term for momentum
  g_v = 150 * avg(Y, t.val);
   
  [uf,vf,wf, uf,vf,wf] = calc_uvw(uf,vf,wf, uf,vf,wf,           ...
                                  prop.rho, prop.mu,            ...
                                  zeros(r.c),                   ...
                                  zeros(r.u), g_v, zeros(r.w),  ...
                                  dt, ts, dx,dy,dz, obst);
  
  %----------
  % Pressure
  %----------
  p = calc_p(p, uf,vf,wf,          ... % variable and advection velocity
             prop.rho,             ... % density
             dt, dx,dy,dz, obst);      % numerical domain

  %---------------------
  % Velocity correction
  %---------------------
  [uf,vf,wf] = corr_uvw(uf,vf,wf, p, prop.rho, dt, dx,dy,dz, obst);

  % Compute the source for pressure again to check mass error
  err = vol_balance(uf,vf,wf, dx,dy,dz, obst);
  disp( sprintf('maximum mass error after correction  = %12.5e', ...
        max(max(max(abs(err))))));
  
  cfl = cfl_max(uf,vf,wf, dt, dx,dy,dz);
  disp( sprintf('cfl max = %5.3f', cfl ) );
  
%==========================================================================
%
% Visualisation
%
%==========================================================================
  
  %-------------------------
  % Plot in Matlab's figure
  %-------------------------
  if mod(ts,10) == 0

    kp = 2;  
    xp = avg(xn);  
    yp = avg(yn); 
    up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
    vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));
    
    subplot(2,1,1)
      [cnt,hnd] = contourf(xp, yp, t.val(:,:,kp)', linspace(-0.0, 1.0, 11)); 
      hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Temperature at time step %d ', ts),   ...
             sprintf('out of %d\n',                  ndt),  ...
             sprintf('CFL = %5.3f\n',                cfl)]);
      drawnow;
    subplot(2,1,2)
      [cnt,hnd] = contourf(xp, yp, p.val(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
      drawnow;
  end
      
  %--------------------------
  % Export to Tecplot format
  %--------------------------
  if mod(ts,100) == 0
      
    container.var(1:5) = [uf vf wf t p];
    file_name = sprintf('results-plumes-staggered-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);

  end  
  
end
 
