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

% set physical properties
[rho, mu, cap, lambda] = air_properties(1:nx, 1:ny, 1:nz);

% time-stepping parameters
dt     =    0.003;   % time steps
ndt    = 5000;       % number of time steps

% create unknowns; names, positions, sizes and default boundary conditions
uf = create_unknown('u-velocity', X, r.u, 'd');
vf = create_unknown('v-velocity', Y, r.v, 'd');
wf = create_unknown('w-velocity', Z, r.w, 'd');
p  = create_unknown('pressure',   C, r.c, 'n');

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

obst = zeros( r.c );
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
                                  zeros(r.c),                          ...
                                  zeros(r.u), zeros(r.v), zeros(r.w),  ...
                                  dt, ts, dx, dy, dz, obst);

  %----------
  % Pressure
  %----------
  p = calc_p(p, uf, vf, wf,          ... % variable and advection velocity
             rho,                    ... % density
             dt, dx, dy, dz, obst);      % numerical domain
         
  %---------------------
  % Velocity correction
  %---------------------
  [uf, vf, wf] = corr_uvw(uf, vf, wf, p, rho, dt, dx, dy, dz, obst);

  % Compute the source for pressure again to check mass error
  err = vol_balance(uf,vf,wf, dx,dy,dz, obst);
  disp( sprintf('maximum mass error after correction  = %12.5e', ...
        max(max(max(abs(err))))));
  
  cfl = cfl_max( uf, vf, wf, dt, dx, dy, dz );
  disp( sprintf('cfl max = %5.3f', cfl ) );
  
  %----------------
  % Plot in Matlab
  %----------------
  if mod(ts,20) == 0
    k_plot = 2;
    xp = avg(xn);  
    yp = avg(yn); 
    up = avg(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val));
    vp = avg(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val));

    subplot(2,1,1)
      [cnt,hnd] = contourf(xp,yp,p.val(:,:,k_plot)'); hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up(:,:,k_plot)', vp(:,:,k_plot)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    subplot(2,1,2)
      [cnt,hnd] = contourf(xp, yp, sqrt([up(:,:,k_plot)'.*up(:,:,k_plot)' ...
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
    
  %--------------------------
  % Export to Tecplot format
  %--------------------------
  if mod(ts,200) == 0
    container.var(1:4) = [uf vf wf p];
    file_name = sprintf('results-obstacles-staggered-%6.6d.dat', ts);
    export_tecplot(file_name, container, r.c, xn, yn, zn);
  end  

end

