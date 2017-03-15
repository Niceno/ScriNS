%==========================================================================
% Demonstrates the membrane for Kerstin.  Two domains are computed
% independently, but linked through boundary conditions in an inner loop
% within each time step.  This seems the most practical approach of all
% because the membrane model implementations are very obvious, at one
% place, in the main function.  The convergence of the conditions at the
% membrane seems to be rather fast too.
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

%------------------
% Phase indicators
%------------------
AIR = 1
H2O = 2

%==========================================================================
%
% Define problem
%
%==========================================================================

% Set domain dimensions and grid resolution for separate domains
x(AIR).nod = nodes(0, 0.1,   80);
y(AIR).nod = nodes(0, 0.01,  20);
z(AIR).nod = nodes(0, 0.1,   80);

x(H2O).nod = x(AIR).nod;
y(H2O).nod = nodes(0.01, 0.03,  30);
z(H2O).nod = z(AIR).nod;

% Expand dx, dy and dz over whole domain
for c=AIR:H2O
  [del(c).x del(c).y del(c).z   ...
   ts(c).x  ts(c).y  ts(c).z] = cartesian_grid(x(c).nod,y(c).nod,z(c).nod);
end

% Resolution ranges
for c=AIR:H2O
  r(c).c = [ts(c).x   ts(c).y   ts(c).z  ];  
  r(c).u = [ts(c).x-1 ts(c).y   ts(c).z  ];
  r(c).v = [ts(c).x   ts(c).y-1 ts(c).z  ];
  r(c).w = [ts(c).x   ts(c).y   ts(c).z-1];
end

% Set physical properties for separate domains
[prop(AIR).rho,  ...
 prop(AIR).mu,   ...
 prop(AIR).cp,   ...
 prop(AIR).lambda] = air_properties( 1:ts(AIR).x, 1:ts(AIR).y, 1:ts(AIR).z );

[prop(H2O).rho,  ...
 prop(H2O).mu,   ...
 prop(H2O).cp,   ...
 prop(H2O).lambda] = water_properties( 1:ts(H2O).x, 1:ts(H2O).y, 1:ts(H2O).z );

% Time-stepping parameters
dt     =    0.001;    % time step
ndt    =  200;        % number of time steps

% Create unknowns; names, positions, sizes and default boundary conditions
for c=AIR:H2O
  u(c)     = create_unknown('u-velocity',     X, r(c).u, 'd');
  v(c)     = create_unknown('v-velocity',     Y, r(c).v, 'd');
  w(c)     = create_unknown('w-velocity',     Z, r(c).w, 'd');
  t(c)     = create_unknown('temperature',    C, r(c).c, 'n');
  p(c)     = create_unknown('pressure',       C, r(c).c, 'n');
  p_tot(c) = create_unknown('total-pressure', C, r(c).c, 'n');
end
  
% Specific boundary conditions
for k=1:ts(AIR).z
  v(AIR).bnd(N).val(:,1,k)=-0.005;    
end
u(AIR).bnd(E).type(1,:,:)='o'; 
u(AIR).bnd(E).val (1,:,:)=0.1;

for k=1:ts(H2O).z
  u(H2O).bnd(W).val(1,:,k)=avg(par(0.1, ts(H2O).y+1));    
end
u(H2O).bnd(E).type(1,:,:)='o'; 
u(H2O).bnd(E).val (1,:,:)=0.1;

for c=AIR:H2O
  u(c).bnd(B).type(:,:,1) = 'n'; 
  u(c).bnd(T).type(:,:,1) = 'n'; 
  v(c).bnd(B).type(:,:,1) = 'n'; 
  v(c).bnd(T).type(:,:,1) = 'n'; 
  w(c).bnd(B).type(:,:,1) = 'n'; 
  w(c).bnd(T).type(:,:,1) = 'n'; 
end

for c=1:H2O
  p(c) = adj_n_bnds(p(c));
end

t(H2O).bnd(W).type(1,:,:) = 'd';  t(H2O).bnd(W).val (1,:,:) =  80;
t(AIR).bnd(S).type(:,1,:) = 'd';  t(AIR).bnd(S).val (:,1,:) =  60;
t(H2O).bnd(S).type(:,1,:) = 'd';  t(H2O).bnd(S).val (:,1,:) =  70;  
t(AIR).bnd(N).type(:,1,:) = 'd';  t(AIR).bnd(N).val (:,1,:) =  70;

t(H2O).val(:,:,:) = 70;
t(AIR).val(:,:,:) = 50;

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

  %------------------------
  % Temperature (enthalpy)
  %------------------------
  for c=AIR:H2O
    t(c).old = t(c).val;
  end
  
  % Initialize temperatures of the membrane
  const =    prop(H2O).lambda(:,  1,:) .* del(AIR).y(:,end,:)  ...
          ./ prop(AIR).lambda(:,end,:) ./ del(H2O).y(:,  1,:);
  t_mem_old = (const.*t(H2O).val(:,1,:)+t(AIR).val(:,end,:)) ./ (1+const);     
  
  while 1

    for c=AIR:H2O
      t(c) = calc_t(t(c), u(c),v(c),w(c),     ...
                    prop(c).rho .* prop(c).cp,  ...
                    prop(c).lambda,             ... 
                    dt, ts, del(c).x,del(c).y,del(c).z, obst);          
    end

    % Compute new temperature of the membrane
    const =    prop(H2O).lambda(:,  1,:) .* del(AIR).y(:,end,:)  ...
            ./ prop(AIR).lambda(:,end,:) ./ del(H2O).y(:,  1,:);
    t_mem = (const.*t(H2O).val(:,1,:)+t(AIR).val(:,end,:)) ./ (1+const);     
      
    % Update boundary conditions with membrane's temperature
    t(H2O).bnd(S).val (:,1,:) =  t_mem;  
    t(AIR).bnd(N).val (:,1,:) =  t_mem;
      
    % Check if convergence has been reached
    eps = max(max(abs(t_mem-t_mem_old)));
    disp( sprintf('eps max = %5.3e', eps ) );
    if eps < TOL
      break
    end  
      
    t_mem_old = t_mem; 
  end
  
  %----------
  % Momentum
  %----------

  for c=AIR:H2O
      
    % Gravity term for momentum
    g_v = - G * avg(Y, prop(c).rho) * min(ts/100,1);
  
    [u(c),v(c),w(c)] = calc_uvw(u(c),v(c),w(c), u(c),v(c),w(c),     ...
                                prop(c).rho, prop(c).mu,            ...
                                p_tot(c).val,                       ...
                                zeros(r(c).u), g_v, zeros(r(c).w),  ...
                                dt, ts, del(c).x,del(c).y,del(c).z, obst);
  end
  
  %----------
  % Pressure
  %----------
  for c=AIR:H2O
      
    p(c) = calc_p(p(c), u(c),v(c),w(c),               ... 
                  prop(c).rho,                        ... 
                  dt, del(c).x,del(c).y,del(c).z, obst);      
            
    p_tot(c).val = p_tot(c).val + p(c).val;
  end
  
  %---------------------
  % Velocity correction
  %---------------------
  for c=AIR:H2O
      
    [u(c),v(c),w(c)] = corr_uvw(u(c),v(c),w(c), p(c),  ...
                                prop(c).rho,           ... 
                                dt, del(c).x,del(c).y,del(c).z, obst);
 
    % Compute the source for pressure again to check mass error
    err = vol_balance(u(c),v(c),w(c), del(c).x,del(c).y,del(c).z, obst);
    disp( sprintf('maximum mass error after correction  = %12.5e', ...
          max(max(max(abs(err))))));
                                 
    cfl = cfl_max(u(c),v(c),w(c), dt, del(c).x,del(c).y,del(c).z);
    disp( sprintf('cfl max = %5.3f', cfl ) );
 
  end  
  
%==========================================================================
%
% Visualisation
%
%==========================================================================
  
  if mod(ts,20) == 0

    %-------------------------
    % Plot to Matlab's figure
    %-------------------------
    kp = 3;
    
    % Catenate coordinates
    xp = avg(x(AIR).nod);     
    yp = cat(2, avg(y(AIR).nod), avg(y(H2O).nod));

    for c=AIR:H2O
      up(c).c = avg(X, cat(X, u(c).bnd(W).val, u(c).val, u(c).bnd(E).val));
      vp(c).c = avg(Y, cat(Y, v(c).bnd(S).val, v(c).val, v(c).bnd(N).val));
    end  
    
    up_c = cat(2, up(AIR).c, up(H2O).c );
    vp_c = cat(2, vp(AIR).c, vp(H2O).c );
    
    tp = cat(2,     t(AIR).val,     t(H2O).val);
    pp = cat(2, p_tot(AIR).val, p_tot(H2O).val);
    
    subplot(2,1,1)
      [cnt,hnd] = contourf(xp, yp, tp(:,:,kp)', linspace(49.999,80.001,10) ); 
      hold on;
%     set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2);
      quiver(xp, yp, up_c(:,:,kp)', vp_c(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Temperature at time step %d ', ts),   ...
             sprintf('out of %d\n',                  ndt),  ...
             sprintf('CFL = %5.3f\n',                cfl)]);
    subplot(2,1,2)
      [cnt,hnd] = contourf(xp, yp, pp(:,:,kp)'); 
      hold on;
      set(hnd,'ShowText','on','TextStep',get(hnd,'LevelStep')*2);
      quiver(xp, yp, up_c(:,:,kp)', vp_c(:,:,kp)', 2, 'k-'); 
      hold off;
      axis equal
      title([sprintf('Pressure correction at time step %d ', ts),   ...
             sprintf('out of %d\n',                          ndt),  ...
             sprintf('CFL = %5.3f\n',                        cfl)]);
    drawnow;
    
    %--------------------------
    % Export to Tecplot format
    %--------------------------
%     container.var(1:4) = [u(AIR), v(AIR), w(AIR), p(AIR)];
%     file_name = sprintf('results-%6.6d.dat', n);
%     export_tecplot(file_name, container, r(AIR).c, [dx(AIR) dy(AIR) dz(AIR)]);

  end
end


