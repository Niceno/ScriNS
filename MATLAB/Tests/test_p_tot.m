%==========================================================================
% Demonstrates the variation of projection algorythm which computes total
% pressure, as a sum of all pressure corrections.  
%
% The total pressure being built up in this way counter-balances the 
% gravity term in momentum equations.
%
% It seems that such an approach is important for bouyancy dominated flows.
%
% Gravity term is under-relaxed here, but it works even without it.
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
xn = nodes(0, 10,   80);
yn = nodes(0,  1,   20);
zn = nodes(0,  0.5,  5);

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% Set physical properties
rho   (1:nx,1:ny,1:nz) = 1.25;    % density 
mu    (1:nx,1:ny,1:nz) = 0.1;     % viscosity

% Time-stepping parameters
dt     =    0.1;    % time step
ndt    = 200;       % number of time steps

% Create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('u-velocity', X, r.u, 'd');
vf = create_unknown('v-velocity', Y, r.v, 'd');
wf = create_unknown('w-velocity', Z, r.w, 'd');
p  = create_unknown('pressure',   C, r.c, 'n');
ptot = zeros( r.c );

% Specific boundary conditions
uf.bnd(W).type(1,:,:)='d'; 
for k=1:nz
  uf.bnd(W).val(1,:,k)=avg(par(0.1, ny+1));    
end
uf.bnd(E).type(1,:,:)='o'; 
uf.bnd(E).val (1,:,:)=0.1;

uf.bnd(B).type(:,:,1) = 'n'; 
uf.bnd(T).type(:,:,1) = 'n'; 
vf.bnd(B).type(:,:,1) = 'n'; 
vf.bnd(T).type(:,:,1) = 'n'; 
wf.bnd(B).type(:,:,1) = 'n'; 
wf.bnd(T).type(:,:,1) = 'n'; 

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
for ts = 1:ndt
  
  disp_time_step(ts);

  %----------
  % Momentum
  %----------

  % Gravity term for momentum
  g_v = - G * avg(Y, rho) * min(ts/100,1);
  
  [uf,vf,wf, uf,vf,wf] = calc_uvw(uf,vf,wf, uf,vf,wf,           ...
                                  rho, mu,                      ...
                                  ptot,                        ...
                                  zeros(r.u), g_v, zeros(r.w),  ...
                                  dt, ts, dx,dy,dz, obst);
  
  %----------
  % Pressure
  %----------
  p = calc_p(p, uf,vf,wf,          ... % variable and advection velocity
             rho,                    ... % density
             dt, dx,dy,dz, obst);      % numerical domain
            
  ptot = ptot + p.val;
  
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
% Visualisation
%
%==========================================================================
  
  %-------------------------
  % Plot to Matlab's figure
  %-------------------------
  if mod(ts,20) == 0
    kp = 3;
    xp = avg(xn);
    yp = avg(yn); 
    up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
    vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));
    subplot(2,1,1)
      [cnt,hnd] = contourf(xp,yp,ptot(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Total pressureat time step %d ', ts),   ...
             sprintf('out of %d\n',                    ndt),  ...
             sprintf('CFL = %5.3f\n',                  cfl)]);
    subplot(2,1,2)
      [cnt,hnd] = contourf(xp,yp,p.val(:,:,kp)'); hold on;
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
    
    container.var(1:4) = [uf,vf,wf, p];
    file_name = sprintf('results-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);
  end
  
end


