%==========================================================================
%
%
%
%--------------------------------------------------------------------------
clear;

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
xn = nodes(0, 0.1,    120);
yn = nodes(0, 0.1,    120);
zn = nodes(0, 0.003,    4);

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];  
r.u = [nx-1 ny   nz  ];  
r.v = [nx   ny-1 nz  ];  
r.w = [nx   ny   nz-1];  

% Set physical properties
[prop.rho,  ...
 prop.mu,   ...
 prop.cap,  ...
 prop.lambda] = air_properties(1:nx, 1:ny, 1:nz);

% time-stepping parameters
dt     =     0.0005;  % time steps
ndt    =  1000;     % number of time steps

% create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('face-u-vel',  X, r.u, 'd');
vf = create_unknown('face-v-vel',  Y, r.v, 'd');
wf = create_unknown('face-w-vel',  Z, r.w, 'd');
p  = create_unknown('pressure',    C, r.c, 'n');

uf.bnd(B).type(:,:,1) = 'n'; 
uf.bnd(T).type(:,:,1) = 'n'; 
vf.bnd(B).type(:,:,1) = 'n'; 
vf.bnd(T).type(:,:,1) = 'n'; 
wf.bnd(B).type(:,:,1) = 'n'; 
wf.bnd(T).type(:,:,1) = 'n'; 

p = adj_n_bnds(p);

% Create a disc
obst = ones(r.c);
for k=1:nz
  for j=1:ny
    for i=1:nx
      dist = sqrt( (j-ny/2)^2 + (i-nx/2)^2 );
      if dist < 0.52 * ny
        obst(i,j,k) = 0;       
      end
    end
  end
end

% Specify boundary conditions
for k=1:nz
  uf.bnd(W).val (1, 3*ny/8+1:5*ny/8, k) = +avg(par(0.5,ny/4+1));    
  uf.bnd(E).type(1, 3*ny/8+1:5*ny/8, k) = 'o';  
  vf.bnd(S).type(3*nx/8+1:5*nx/8, 1, k) = 'o';   
  vf.bnd(N).type(3*nx/8+1:5*nx/8, 1, k) = 'o';   
end  
    
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
  
  %----------
  % Momentum
  %----------
  [uf, vf, wf, uf, vf, wf] =                                      ...
                     calc_uvw(uf, vf, wf, uf, vf, wf,             ...
                              prop.rho, prop.mu,                  ...
                              zeros(r.c),                         ...
                              zeros(r.u), zeros(r.v), zeros(r.w), ...
                              dt, ts, dx, dy, dz, obst);

  %----------
  % Pressure
  %----------
  p = calc_p(p, uf, vf, wf,          ... % variable and advection velocity
             prop.rho,                    ... % density
             dt, dx, dy, dz, obst);      % numerical domain

  %---------------------
  % Velocity correction
  %---------------------
  [uf,vf,wf] = corr_uvw(uf,vf,wf, p, prop.rho, dt, dx,dy,dz, obst);

  % Compute the source for pressure again to check mass error
  err = vol_balance(uf,vf,wf, dx,dy,dz, obst);
  disp( sprintf('maximum mass error after correction  = %12.5e', ...
        max(max(max(abs(err))))));
  
  cfl = cfl_max( uf, vf, wf, dt, dx, dy, dz );
  disp( sprintf('cfl max = %5.3f', cfl ) );
  
%==========================================================================
%
% Visualisation
%
%==========================================================================

  %-----------------------
  % Plot in Matlab format
  %-----------------------
  if mod(ts,10) == 0
    kp = 2;
    xp = avg(xn);  
    yp = avg(yn); 
    up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
    vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));
    subplot(1,2,1)
      [cnt,hnd] = contourf(xp,yp,p.val(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 4, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    subplot(1,2,2)
      [cnt,hnd] = contourf(xp, yp, sqrt([   up(:,:,kp)'.*up(:,:,kp)' ...
                                          + vp(:,:,kp)'.*vp(:,:,kp)'])); 
      hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 4, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Velocity magnitude at time step %d ', ts),   ...
             sprintf('out of %d\n',                         ndt),  ...
             sprintf('CFL = %5.3f\n',                       cfl)]);
    drawnow;
  end
end
