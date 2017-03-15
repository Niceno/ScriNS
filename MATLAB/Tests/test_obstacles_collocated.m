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
xn = nodes(0, 1,     256);
yn = nodes(0, 0.125,  32);
zn = nodes(0, 0.125,   4);

% Expand dx, dy and dz over whole domain
[dx dy dz nx ny nz] = cartesian_grid(xn, yn, zn);

% Resolution ranges
r.c = [nx   ny   nz  ];
r.u = [nx-1 ny   nz  ];
r.v = [nx   ny-1 nz  ];
r.w = [nx   ny   nz-1];

% Set physical properties
[rho, mu, cap, lambda] = air_properties(1:nx, 1:ny, 1:nz);

% Time-stepping parameters
dt     =    0.003;   % time steps
ndt    = 5000;       % number of time steps

% Create unknowns; names, positions, sizes and default boundary conditions
uc = create_unknown('u-velocity',  C, r.c, 'd');
vc = create_unknown('v-velocity',  C, r.c, 'd');
wc = create_unknown('w-velocity',  C, r.c, 'd');
uf = create_unknown('face-u-vel',  X, r.u, 'd');
vf = create_unknown('face-v-vel',  Y, r.v, 'd');
wf = create_unknown('face-w-vel',  Z, r.w, 'd');
p  = create_unknown('pressure',    C, r.c, 'n');

% Specify boundary conditions
for k=1:nz
  uc.bnd(W).val(1,:,k)=avg(par(0.1, ny+1));    
end
uc.bnd(E).type(1,:,:) = 'o'; 

uc.bnd(B).type(:,:,1) = 'n'; 
uc.bnd(T).type(:,:,1) = 'n'; 
vc.bnd(B).type(:,:,1) = 'n'; 
vc.bnd(T).type(:,:,1) = 'n'; 
wc.bnd(B).type(:,:,1) = 'n'; 
wc.bnd(T).type(:,:,1) = 'n'; 

p = adj_n_bnds(p);

% Copy the values to face velocities 
uf.val = avg(X, uc.val);  uf.bnd = uc.bnd;
vf.val = avg(Y, vc.val);  vf.bnd = vc.bnd;
wf.val = avg(Z, wc.val);  wf.bnd = wc.bnd;
 
% Create a slanted obstacle in the flow domain
obst = zeros(r.c);
for o=1:4
  for j=1:3*ny/4;
    for i=nx/4+j:nx/4+3*ny/4; 
      for k=1:nz  
        obst(i,j,k) = 1;
      end  
    end
  end
end

%-----------
%
% time loop 
%
%-----------
for ts = 1:ndt

  disp_time_step(ts);
  
  %----------
  % Momentum
  %----------
  [uc,vc,wc, uf,vf,wf] = calc_uvw(uc,vc,wc, uf,vf,wf,                  ...
                                  rho, mu,                             ...
                                  zeros(r.c),                          ...
                                  zeros(r.c), zeros(r.c), zeros(r.c),  ...
                                  dt, ts, dx,dy,dz, obst);

  %----------
  % Pressure
  %----------
  p = calc_p(p, uf,vf,wf,            ... % variable and advection velocity
             rho,                    ... % density
             dt, dx, dy, dz, obst);      % numerical domain
         
  %---------------------
  % Velocity correction
  %---------------------
  
  % Correct the cell centered velocites
  [uc,vc,wc] = corr_uvw(uc,vc,wc, p, rho, dt, dx,dy,dz, obst);
    
  % Correct face velocities in staggered fashion 
  [uf,vf,wf] = corr_uvw(uf,vf,wf, p, rho, dt, dx,dy,dz, obst);

  % Compute the source for pressure again to check mass error
  err = vol_balance(uf,vf,wf, dx,dy,dz, obst);
  disp( sprintf('maximum mass error after correction  = %12.5e', ...
        max(max(max(abs(err))))));
  
  cfl = cfl_max(uc,vc,wc, dt, dx,dy,dz);
  disp( sprintf('cfl max = %5.3f', cfl ) );
  
  %----------------
  % Plot in Matlab
  %----------------
  if mod(ts,20) == 0
    kp = 2;
    xp = avg(xn);  
    yp = avg(yn); 
    up = uc.val;
    vp = vc.val;

    subplot(2,1,1)
%      [cnt,hnd] = contourf(xp, yp, p.val(:,:,kp)'); hold on;
      [cnt,hnd] = contourf(xp, yp, err(:,:,kp)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    subplot(2,1,2)
      [cnt,hnd] = contourf(xp, yp, sqrt(  [up(:,:,kp)'.*up(:,:,kp)'  ...
                                         + vp(:,:,kp)'.*vp(:,:,kp)'])); 
      hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,kp)', vp(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Velocity magnitude at time step %d ', ts),   ...
             sprintf('out of %d\n',                         ndt),  ...
             sprintf('CFL = %5.3f\n',                       cfl)]);
    drawnow;
  end
  
  %--------------------------
  % Export to Tecplot format
  %--------------------------
  if mod(ts,200) == 0
    container.var(1:4) = [uc vc wc p];
    file_name = sprintf('results-obstacles-collocated-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);
  end
    
end


