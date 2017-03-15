%==========================================================================
% Tests the rotation of a shape (defined in file shape.dat or in the script 
% itself) using advection schemes implemented in the code.
%
% Name of the scheme to be used is set with the variable "NAME"
%--------------------------------------------------------------------------
clear;

NAME = 'superbee'

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
xn = nodes(-1, 1,   120);
yn = nodes(-1, 1,   120);
zn = nodes( 0, 0.1,   4);

% Cell dimensions
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];  
r.u = [nx-1 ny   nz  ];  
r.v = [nx   ny-1 nz  ];  
r.w = [nx   ny   nz-1];  

% Set physical properties
prop.rho   (1:nx, 1:ny, 1:nz) = 1.0e3;   % density 
prop.mu    (1:nx, 1:ny, 1:nz) = 0.01;    % viscosity
prop.lambda(1:nx, 1:ny, 1:nz) = 0.00000; % thermal conductivity
prop.cap   (1:nx, 1:ny, 1:nz) = 1.0;     % thermal capacity

% Time-stepping parameters
dt     = pi/600;    % time steps
ndt    = 300;       % number of time steps

% Create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('face-u-vel',  X, r.u, 'd');
vf = create_unknown('face-v-vel',  Y, r.v, 'd');
wf = create_unknown('face-w-vel',  Z, r.w, 'd');
tc = create_unknown('temperature', C, r.c, 'd');

% Create rotational velocity field
for k=1:nz  
  for j=1:ny  
    for i=1:nx-1  
      uf.val(i,j,k) = -yn(j+1);
    end  
  end  
  for j=1:ny-1  
    for i=1:nx  
      vf.val(i,j,k) =  xn(i+1);
    end  
  end  
end  

% Initial and boundary conditions
SHAPE = ...
[ [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]; ... 
  [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ... 
  [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0]; ...
  [0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0]; ... 
  [0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0]; ...
  [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0]; ...
  [0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0]; ...
  [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0]; ...
  [0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0] ];
% load('shape.dat', 'shape', '-ascii')
[ni nj] = size(SHAPE);

for k=1:nz
  tx.val(nx/2-ni/2+1:nx/2+ni/2, 3*ny/4-ni/2+1 : 3*ny/4+ni/2, k) = SHAPE;
  ty.val(nx/2-ni/2+1:nx/2+ni/2, 3*ny/4-ni/2+1 : 3*ny/4+ni/2, k) = SHAPE;
  tc.val(nx/2-ni/2+1:nx/2+ni/2, 3*ny/4-ni/2+1 : 3*ny/4+ni/2, k) = SHAPE;
end

obst = [];

%-----------
%
% Time loop 
%
%-----------
for ts=1:ndt

  disp_time_step( ts );

  %---------------------
  % Collocated variable
  %---------------------
  tc = calc_t(tc, uf, vf, wf,             ... % variable and advection velocity
              prop.rho, prop.mu,          ... % physical properties
              dt, ts, dx, dy, dz, obst);      % numerical domain      
  
  cfl = cfl_max( uf, vf, wf, dt, dx, dy, dz );
  disp( sprintf('cfl max = %5.3f', cfl ) );

  %------
  % Plot
  %------
  kp = 2;
  up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
  vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));
  xp = avg(xn);
  yp = avg(yn);
  tp = tc.val(:,:,kp);
  [cnt,hnd] = contourf(xp, yp, tp', linspace(-0.01,1.01,3));
  hold on;
  set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
  quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
  hold off;
  axis equal
  title(sprintf('Time step %d out of %d\nCFL = %5.3f\n', ts, ndt, cfl))
  drawnow;
end

