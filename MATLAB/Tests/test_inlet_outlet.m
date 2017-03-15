%==========================================================================
% This scripts tests inlet and outlet conditons at various places in
% computational domain.
%
% Computational domain is a simple box, and Reynolds number is rather 
% small to avoid instabilities due to vortices getting out from the 
% outlet. The script selects the case it will run randomply, from the 
% 16 predefined cases.
%--------------------------------------------------------------------------
clear;

tests = [11:14,21:24,31:34,41:44];
TEST = tests( ceil(rand*16) );

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
xn = nodes(0, 1,    160);
yn = nodes(0, 1,    160);
zn = nodes(0, 0.025,  4);

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];  
r.u = [nx-1 ny   nz  ];  
r.v = [nx   ny-1 nz  ];  
r.w = [nx   ny   nz-1];  

% set physical properties
[prop.rho,  ...
 prop.mu,   ...
 prop.cap,  ...
 prop.lambda] = air_properties(1:nx, 1:ny, 1:nz);


% time-stepping parameters
dt     =     0.15;  % time steps
ndt    =   200;     % number of time steps

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

% Create a cylindrical obstacle in the middle just for kicks
obst = zeros(r.c);
for k=1:nz
  for j=1:ny
    for i=1:nx
      dist = sqrt( (j-ny/2)^2 + (i-nx/2)^2 );
      if dist < ny/4
        obst(i,j,k) = 1;       
      end
    end
  end
end

% specify boundary conditions
if TEST == 11
  for k=1:nz
    uf.bnd(W).val (1,ny/4+1:3*ny/4,k) = +avg(par(0.01,ny/2+1));    
    uf.bnd(E).type(1,ny/4+1:3*ny/4,k) = 'o';  
  end  
elseif TEST == 12 % vertical mirror from 11
  for k=1:nz
    uf.bnd(E).val (1,ny/4+1:3*ny/4,k) = -avg(par(0.01,ny/2+1));    
    uf.bnd(W).type(1,ny/4+1:3*ny/4,k) = 'o';   
  end  
elseif TEST == 13 % rotate 11
  for k=1:nz
    vf.bnd(S).val (nx/4+1:3*nx/4,1,k) = +avg(par(0.01,nx/2+1));    
    vf.bnd(N).type(nx/4+1:3*nx/4,1,k) = 'o';   
  end  
elseif TEST == 14 % horizontal mirror 13
  for k=1:nz
    vf.bnd(N).val (nx/4+1:3*nx/4,1,k) = -avg(par(0.01,nx/2+1));    
    vf.bnd(S).type(nx/4+1:3*nx/4,1,k) = 'o';   
  end  
elseif TEST == 21 % 2 exits
  for k=1:nz
    uf.bnd(W).val (1,  ny/4+1:3*ny/4,k) = +avg(par(0.01,ny/2+1));    
    uf.bnd(E).type(1,       1:ny/4+1,k) = 'o';   
    uf.bnd(E).type(1,3*ny/4+1:ny,    k) = 'o';   
  end  
elseif TEST == 22 % vertical mirror 21
  for k=1:nz
    uf.bnd(E).val (1,  ny/4+1:3*ny/4,k) = -avg(par(0.01,ny/2+1));    
    uf.bnd(W).type(1,       1:ny/4+1,k) = 'o';   
    uf.bnd(W).type(1,3*ny/4+1:ny,    k) = 'o';   
  end  
elseif TEST == 23 % rotated 21
  for k=1:nz
    vf.bnd(S).val (  nx/4+1:3*nx/4,1,k) = +avg(par(0.01,nx/2+1));    
    vf.bnd(N).type(       1:nx/4+1,1,k) = 'o';   
    vf.bnd(N).type(3*nx/4+1:nx,    1,k) = 'o';   
  end  
elseif TEST == 24 % horizontal mirror of 23
  for k=1:nz
    vf.bnd(N).val (  nx/4+1:3*nx/4,1,k) = -avg(par(0.01,nx/2+1));    
    vf.bnd(S).type(       1:nx/4+1,1,k) = 'o';   
    vf.bnd(S).type(3*nx/4+1:nx,    1,k) = 'o';   
  end  
elseif TEST == 31 % inlet and outlet at the same face
  for k=1:nz
    uf.bnd(W).val (1,3*ny/4+1:ny,  k) = +avg(par(0.01,ny/4+1));    
    uf.bnd(W).type(1,       1:ny/4,k) = 'o';   
  end  
elseif TEST == 32 % vertical mirror of 31
  for k=1:4
    uf.bnd(E).val (1,3*ny/4+1:ny,  k) = -avg(par(0.01,ny/4+1));    
    uf.bnd(E).type(1,       1:ny/4,k) = 'o';   
  end  
elseif TEST == 33 % rotated 31
  for k=1:4
    vf.bnd(S).val (3*nx/4+1:nx,  1,k) = +avg(par(0.01,nx/4+1));    
    vf.bnd(S).type(       1:nx/4,1,k) = 'o';   
  end  
elseif TEST == 34 % horizontal mirror of 33
  for k=1:4
    vf.bnd(N).val (3*nx/4+1:nx,  1,k) = -avg(par(0.01,nx/4+1));    
    vf.bnd(N).type(       1:nx/4,1,k) = 'o';   
  end  
elseif TEST == 41 % inlet and outlet at the same face, one more outlet
  for k=1:4
    uf.bnd(W).val (1,3*ny/4+1:ny,  k) = +avg(par(0.01,ny/4+1));    
    uf.bnd(W).type(1,       1:ny/8,k) = 'o';   
    uf.bnd(E).type(1,       1:ny/8,k) = 'o';   
  end  
elseif TEST == 42 % vertical mirror of 41
  for k=1:4
    uf.bnd(E).val (1,3*ny/4+1:ny,  k) = -avg(par(0.01,ny/4+1));    
    uf.bnd(E).type(1,       1:ny/8,k) = 'o';   
    uf.bnd(W).type(1,       1:ny/8,k) = 'o';   
  end  
elseif TEST == 43 % rotated 41
  for k=1:4
    vf.bnd(S).val (3*nx/4+1:nx,  1,k) = +avg(par(0.01,nx/4+1));    
    vf.bnd(S).type(       1:nx/8,1,k) = 'o';   
    vf.bnd(N).type(       1:nx/8,1,k) = 'o';   
  end  
elseif TEST == 44 % horizontal mirror of 43
  for k=1:4
    vf.bnd(N).val (3*ny/4+1:nx,  1,k) = -avg(par(0.01,nx/4+1));    
    vf.bnd(N).type(       1:nx/8,1,k) = 'o';   
    vf.bnd(S).type(       1:nx/8,1,k) = 'o';   
  end  
end    
    
%==========================================================================
%
% Solution algorithm
%
%==========================================================================

%-----------
%
% time loop 
%
%-----------
for ts=1:ndt

  disp_time_step(ts);
  
  %----------
  % Momentum
  %----------
  [uf, vf, wf, uf, vf, wf] =                                            ...
                    calc_uvw(uf, vf, wf, uf, vf, wf,                    ...
                             prop.rho, prop.mu,                         ...
                             zeros( r.c ),                              ...
                             zeros( r.u ), zeros( r.v ), zeros( r.w ),  ...
                             dt, ts, dx, dy, dz, obst);

  %----------
  % Pressure
  %----------
  p = calc_p(p, uf, vf, wf,          ... % variable and advection velocity
             prop.rho,               ... % density
             dt, dx, dy, dz, obst);      % numerical domain

  %---------------------
  % Velocity correction
  %---------------------
  [uf, vf, wf] = corr_uvw(uf,vf,wf, p, prop.rho, dt, dx,dy,dz, obst);

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

  %----------------
  % Plot in Matlab
  %----------------
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
      title([sprintf('Test %2d: ',                           TEST),  ...
             sprintf('Pressure correction at time step %d ', ts),    ...
             sprintf('out of %d\n',                          ndt),   ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    subplot(1,2,2)
      [cnt,hnd] = contourf(xp, yp, sqrt(  [up(:,:,kp)'.*up(:,:,kp)' ...
                                         + vp(:,:,kp)'.*vp(:,:,kp)'])); 
      hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 4, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Test %2d: ',                          TEST),  ...
             sprintf('Velocity magnitude at time step %d ', ts),    ...
             sprintf('out of %d\n',                         ndt),   ...
             sprintf('CFL = %5.3f\n',                       cfl)]);
    drawnow;
  end
end
