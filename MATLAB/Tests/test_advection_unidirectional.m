%==========================================================================
% Program to test implementation of advection schemes in the code, using
% one-dimensoinal transport of a step function in X, Y or Z 
% direction, either in positive or negative sense.
%
% The coordinate direction is specified with the local variable "TEST",
% which can be either X, Y or Z.
%
% Sense is specified with the variable "FLOW", which can assume values 'p'
% for positive, and 'n' for negative sense.
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

%---------------------
% Set local constants
%---------------------
test = [X, Y, Z];  TEST = test( ceil(rand*3) );
flow = ['p','n'];  FLOW = flow( ceil(rand*2) );

NAME = 'koren'
FLOW = 'n'
MIN = 10;
MAX = 20;

%==========================================================================
%
% Define problem
%
%==========================================================================

% Set domain dimensions and grid resolution
if TEST == X
  xn = nodes(0,  1,   200);
  yn = nodes(0,  0.1,   4);
  zn = nodes(0,  0.1,   4);
  if FLOW == 'p'
    u_bulk = +1.0;  v_bulk =  0.0;  w_bulk =  0.0;
  elseif FLOW == 'n'
    u_bulk = -1.0;  v_bulk =  0.0;  w_bulk =  0.0;
  end  
elseif TEST == Y
  xn = nodes(0,  0.1,   4);
  yn = nodes(0,  1,   200);
  zn = nodes(0,  0.1,   4);
  if FLOW == 'p'
    u_bulk =  0.0;  v_bulk = +1.0;  w_bulk =  0.0;
  elseif FLOW == 'n'
    u_bulk =  0.0;  v_bulk = -1.0;  w_bulk =  0.0;
  end  
elseif TEST == Z
  xn = nodes(0,  0.1,   4);
  yn = nodes(0,  0.1,   4);
  zn = nodes(0,  1,   200);
  if FLOW == 'p'
    u_bulk =  0.0;  v_bulk =  0.0;  w_bulk = +1.0;
  elseif FLOW == 'n'
    u_bulk =  0.0;  v_bulk =  0.0;  w_bulk = -1.0;
  end  
end

% Cell dimensions
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% Set physical properties
rho   (1:nx,1:ny,1:nz) = 1.0; 
mu    (1:nx,1:ny,1:nz) = 0.0;
lambda(1:nx,1:ny,1:nz) = 0.0; 
cap   (1:nx,1:ny,1:nz) = 1.0; 

% time-stepping parameters
dt     =   0.002;    % time step
ndt    = 350;        % number of time steps

% create unknowns; names, positions and sizes
uf = create_unknown('face-u-vel',  X, r.u, 'd');
vf = create_unknown('face-v-vel',  Y, r.v, 'd');
wf = create_unknown('face-w-vel',  Z, r.w, 'd');
t  = create_unknown('temperature', C, r.c, 'n');

% specify inital and boundary conditions
uf.val(:,:,:) = u_bulk;
vf.val(:,:,:) = v_bulk;
wf.val(:,:,:) = w_bulk;
for j=1:6
  uf.bnd(j).val(:,:,:) = u_bulk;    
  vf.bnd(j).val(:,:,:) = v_bulk;    
  wf.bnd(j).val(:,:,:) = w_bulk;
end

t.val(:,:,:) = MIN;
for j=1:6
  t.bnd(j).val(:,:,:) = MIN;
end  
if TEST == X  t.val(nx/2-25:nx/2+25, :, :) = MAX;  end
if TEST == Y  t.val(:, ny/2-25:ny/2+25, :) = MAX;  end
if TEST == Z  t.val(:, :, nz/2-25:nz/2+25) = MAX;  end
  
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
  t = calc_t(t, uf, vf, wf,              ... % variable and adv. velocity
             (rho.*cap), lambda,         ... % physical properties
             dt, ts, dx, dy, dz, obst);      % numerical domain
  
%==========================================================================
%
% Visualization
%
%==========================================================================

  %------
  % Plot
  %------
  clf
  xp = avg(xn);  
  yp = avg(yn);
  zp = avg(zn);
  if TEST == X
    plot(xp, reshape( t.val(:,ny/2,nz/2), nx,1 ), 'r.-')
    axis([-0.1 1.1 MIN*0.8 MAX*1.2])
  end  
  if TEST == Y
    plot(yp, reshape( t.val(nx/2,:,nz/2), ny,1 ), 'g.-')
    axis([-0.1 1.1 MIN*0.8 MAX*1.2])
  end  
  if TEST == Z
    plot(zp, reshape( t.val(nx/2,ny/2,:), nz,1 ), 'b.-')
    axis([-0.1 +1.1 MIN*0.8 MAX*1.2])
  end  
  title([sprintf('Transport with %s ', NAME), ...
         sprintf('in %s direction, ',  TEST), ...
         sprintf('%s sense',           FLOW)]);
  drawnow

  cfl = cfl_max( uf, vf, wf, dt, dx, dy, dz );
  disp( sprintf('cfl max = %5.3f', cfl ) );
end



