%==========================================================================
% This program solves thermally driven cavity at Ra = 1.0e6, 
% in dimensional and non-dimensional forms.
%
% Equations in dimensional form:
%
% D(rho u)/Dt = nabla(mu (nabla u)^T) - nabla p + g 
% D(rho cp T)/Dt = nabla(lambda (nabla T)^T) 
%
% Equations in non-dimensional form, for natural convection problems
%
% DU/Dt = nabla(1/sqrt(Gr) (nabla U)^T) - nabla P + theta 
% D theta/Dt = nabla(1/(Pr*sqrt(Gr)) (nabla theta)^T) 
%
%--------------------------------------------------------------------------
% For thermally driven cavity, with properties of air at 60 deg:
%
% nu   =  1.89035E-05;
% beta =  0.003;
% dT   = 17.126;
% L    =  0.1;
% Pr   = 0.709;
%
% characteristic non-dimensional numbers are:
% Gr = 1.4105E+06
% Ra = 1.0000E+06
%--------------------------------------------------------------------------
clear;

%----------------------------------------------------
% Chose between dimensional or non-dimensional forms
%----------------------------------------------------
FORM = 'n'; % 'd' or 'n' 

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
if FORM == 'd'
  L  = 0.1;
else  
  L  = 1.0;
end
xn = nodes(0, L,    64, L/256, L/256);  xc = avg(xn);
yn = nodes(0, L,    64, L/256, L/256);  yc = avg(yn);
zn = nodes(0, L/10,  5);                zc = avg(zn);   

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% Set physical properties
if FORM == 'd'
  [rho, mu, cap, lambda] = air_properties(1:nx, 1:ny, 1:nz);
else
  grashof = 1.4105E+06;
  prandtl = 0.7058;
  rho   (1:nx, 1:ny, 1:nz) = 1.0;
  mu    (1:nx, 1:ny, 1:nz) = 1.0 / sqrt(grashof);
  cap   (1:nx, 1:ny, 1:nz) = 1.0;
  lambda(1:nx, 1:ny, 1:nz) = 1.0 / (prandtl * sqrt(grashof));
end

% Time-stepping parameters
if FORM == 'd'
  dt     =    0.01;  % time steps
  ndt    = 333;     % number of time steps
else
  dt     =    0.02;  % time steps
  ndt    = 333;     % number of time steps
end

% Create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('u-velocity',  X, r.u, 'd');
vf = create_unknown('v-velocity',  Y, r.v, 'd');
wf = create_unknown('w-velocity',  Z, r.w, 'd');
p  = create_unknown('pressure',    C, r.c, 'n');
t  = create_unknown('temperature', C, r.c, 'n');
p_tot = zeros(nx, ny, nz);

if FORM == 'd'
  t.val(:,:,:) = 60.0;
  t.bnd(W).type(1,:,:) = 'd';  t.bnd(W).val(1,:,:) = 68.563;
  t.bnd(E).type(1,:,:) = 'd';  t.bnd(E).val(1,:,:) = 51.437;
  t.bnd(S).val (:,1,:) = 60.0;
  t.bnd(N).val (:,1,:) = 60.0;
else  
  t.val(:,:,:) = 0.0;
  t.bnd(W).type(1,:,:) = 'd';  t.bnd(W).val(1,:,:) =  0.5;
  t.bnd(E).type(1,:,:) = 'd';  t.bnd(E).val(1,:,:) = -0.5;
end  

uf.bnd(B).type(:,:,1) = 'n'; 
uf.bnd(T).type(:,:,1) = 'n'; 
vf.bnd(B).type(:,:,1) = 'n'; 
vf.bnd(T).type(:,:,1) = 'n'; 
wf.bnd(B).type(:,:,1) = 'n'; 
wf.bnd(T).type(:,:,1) = 'n'; 

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
tic
for ts=1:ndt
  
  disp_time_step(ts);

  if FORM == 'd'
    rho = 1.127 + (t.val - 40.0) * (1.00 - 1.127) / 40.0;
  end
    
  %------------------------
  % Temperature (enthalpy)
  %------------------------
  t = calc_t(t, uf,vf,wf,         ... % variable and advection velocity
             (rho.*cap), lambda,  ... % physical properties
             dt, ts,              ... % time stepping,
             dx,dy,dz, obst);         % numerical domain
                    
  %----------
  % Momentum
  %----------
 
  % Gravity term for momentum
  if FORM == 'd'
    g_v = -G * avg(Y,rho);
  else
    g_v = avg(Y,t.val);
  end    
    
  [uf,vf,wf, uf,vf,wf] = calc_uvw(uf,vf,wf, uf,vf,wf,           ...
                                  rho, mu,                      ...
                                  p_tot,                        ...
                                  zeros(r.u), g_v, zeros(r.w),  ...
                                  dt, ts,                       ...
                                  dx,dy,dz, obst);
                                       
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
  
  cfl = cfl_max(uf,vf,wf, dt, dx,dy,dz);
  disp( sprintf('cfl max = %5.3f', cfl ) );
  
%==========================================================================
%
% Visualization
%
%==========================================================================

  %-------------------------
  % Plot in Matlab's figure
  %-------------------------
  if mod(ts,1) == 0
    kp = 3;
    xp = avg(xn);  
    yp = avg(yn); 
    up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
    vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));
    
    subplot(1,3,1)
      [cnt,hnd] = contourf(xp, yp, t.val(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Temperature at time step %d ', ts),   ...
             sprintf('out of %d\n',                  ndt),  ...
             sprintf('CFL = %5.3f\n',                cfl)]);
    subplot(1,3,2)
      [cnt,hnd] = contourf(xp, yp, p.val(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    subplot(1,3,3)
      [C,hnd] = contourf(xp, yp, p_tot(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure at time step %d ', ts),   ...
             sprintf('out of %d\n',               ndt),  ...
             sprintf('CFL = %5.3f\n',             cfl)]);
    drawnow;
  end
  
  %--------------------------
  % Export to Tecplot format
  %--------------------------
  if mod(ts,100) == 0
      
    container.var(1:5) = [uf vf wf t p];
    file_name = sprintf('results-thermal-cavity-staggered-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);
    
  end
    
end
toc

