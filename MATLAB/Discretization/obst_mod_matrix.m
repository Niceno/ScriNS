%==========================================================================
function [c] = obst_matrix(phi, c, obst, obc)
%--------------------------------------------------------------------------
% Adjusts the system matrix for obstacles and cell centered varaibles 
% (such as pressure)
% 
% phi  - variable
% c    - coefficients in system matrix
% obst - obstacle array
% obc  - obstacles's boundary condition, ('n' - Neumann, 'd' - Dirichlet)
%--------------------------------------------------------------------------

pos = phi.pos;

%--------------------------
%
% For collocated variables
%
%--------------------------
if pos == C

  %------------------------------------
  % Neumann's boundary on the obstacle
  %------------------------------------
  if obc == 'n'
      
    % Correct west and east
    sol_x = dif(X, obst);  % will be +1 east of obst, -1 west of obst
    corr = 1 - (sol_x < 0);  c.W(2:end,  :,:) = c.W(2:end,  :,:) .* corr;
    corr = 1 - (sol_x > 0);  c.E(1:end-1,:,:) = c.E(1:end-1,:,:) .* corr;
 
    % Correct south and north
    sol_y = dif(Y, obst);  % will be +1 north of obst, -1 south of obst  
    corr = 1 - (sol_y < 0);  c.S(:,2:end,  :) = c.S(:,2:end,  :) .* corr;
    corr = 1 - (sol_y > 0);  c.N(:,1:end-1,:) = c.N(:,1:end-1,:) .* corr;
  
    % Correct bottom and top
    sol_z = dif(Z, obst);  % will be +1 north of obst, -1 south of obst  
    corr = 1 - (sol_z < 0);  c.B(:,:,2:end  ) = c.B(:,:,2:end  ) .* corr;
    corr = 1 - (sol_z > 0);  c.T(:,:,1:end-1) = c.T(:,:,1:end-1) .* corr;

  %--------------------------------------
  % Dirichlet's boundary on the obstacle
  %--------------------------------------
  elseif obc == 'd' 
  
    % Set central coefficient to 1 in obst, unchanged elsewhere
    c.P = c.P .* not(obst) + obst; 

    % Set neighbour coefficients to zero in obst  
    c.W = c.W .* not(obst);
    c.E = c.E .* not(obst);
    c.S = c.S .* not(obst);
    c.N = c.N .* not(obst);
    c.B = c.B .* not(obst);
    c.T = c.T .* not(obst);
      
    % Increase coefficients close to obst (makes sense for momentum)
    sol_x = dif(X, obst);  % will be +1 east of obst, -1 west of obst
    corr = 1 + (sol_x > 0);  c.E(1:end-1,:,:) = c.E(1:end-1,:,:) .* corr;  
    corr = 1 + (sol_x < 0);  c.W(2:end,  :,:) = c.W(2:end,  :,:) .* corr;  
      
    sol_y = dif(Y, obst);  % will be +1 north of obst, -1 south of obst
    corr = 1 + (sol_y > 0);  c.N(:,1:end-1,:) = c.N(:,1:end-1,:) .* corr;  
    corr = 1 + (sol_y < 0);  c.S(:,2:end,  :) = c.S(:,2:end,  :) .* corr;  
      
    sol_z = dif(Z, obst);  % will be +1 top of obst, -1 bottom of obst
    corr = 1 + (sol_z > 0);  c.T(:,:,1:end-1) = c.T(:,:,1:end-1) .* corr;  
    corr = 1 + (sol_z < 0);  c.B(:,:,2:end  ) = c.B(:,:,2:end  ) .* corr;  
  end

%-------------------------
%
% For staggered variables
%
%-------------------------
elseif pos==X 
    
  % Set central coefficient to 1 in obst, unchanged elsewhere
  obst_x = max(obst(1:end-1,:,:), obst(2:end,:,:));
  c.P = c.P .* not(obst_x) + obst_x; 

  % Set neighbour coefficients to zero in obst  
  c.W = c.W .* not(obst_x);
  c.E = c.E .* not(obst_x);
  c.S = c.S .* not(obst_x);
  c.N = c.N .* not(obst_x);
  c.B = c.B .* not(obst_x);
  c.T = c.T .* not(obst_x);
      
  % Increase coefficients close to obst (makes sense for momentum)
  sol_y = dif(Y, obst_x);  % will be +1 north of obst, -1 south of obst
  corr = 1 + (sol_y > 0);  c.N(:,1:end-1,:) = c.N(:,1:end-1,:) .* corr;  
  corr = 1 + (sol_y < 0);  c.S(:,2:end,  :) = c.S(:,2:end,  :) .* corr;  
     
  sol_z = dif(Z, obst_x);  % will be +1 top of obst, -1 bottom of obst
  corr = 1 + (sol_z > 0);  c.T(:,:,1:end-1) = c.T(:,:,1:end-1) .* corr;  
  corr = 1 + (sol_z < 0);  c.B(:,:,2:end  ) = c.B(:,:,2:end  ) .* corr;  
      
elseif pos==Y  
        
  % Set central coefficient to 1 in obst, unchanged elsewhere
  obst_y = max(obst(:,1:end-1,:), obst(:,2:end,:));
  c.P = c.P .* not(obst_y) + obst_y; 

  % Set neighbour coefficients to zero in obst  
  c.W = c.W .* not(obst_y);
  c.E = c.E .* not(obst_y);
  c.S = c.S .* not(obst_y);
  c.N = c.N .* not(obst_y);  
  c.B = c.B .* not(obst_y);
  c.T = c.T .* not(obst_y);  
      
  % Increase coefficients close to obst (makes sense for momentum)
  sol_x = dif(X, obst_y);  % will be +1 north of obst, -1 south of obst
  corr = 1 + (sol_x > 0);  c.E(1:end-1,:,:) = c.E(1:end-1,:,:) .* corr;  
  corr = 1 + (sol_x < 0);  c.W(2:end,  :,:) = c.W(2:end,  :,:) .* corr;   
      
  sol_z = dif(Z, obst_y);  % will be +1 north of obst, -1 south of obst
  corr = 1 + (sol_z > 0);  c.T(:,:,1:end-1) = c.T(:,:,1:end-1) .* corr;  
  corr = 1 + (sol_z < 0);  c.B(:,:,2:end  ) = c.B(:,:,2:end  ) .* corr;   
      
elseif pos==Z  
        
  % Set central coefficient to 1 in obst, unchanged elsewhere
  obst_z = max(obst(:,:,1:end-1), obst(:,:,2:end));
  c.P = c.P .* not(obst_z) + obst_z; 

  % Set neighbour coefficients to zero in obst  
  c.W = c.W .* not(obst_z);
  c.E = c.E .* not(obst_z);
  c.S = c.S .* not(obst_z);
  c.N = c.N .* not(obst_z);  
  c.B = c.B .* not(obst_z);
  c.T = c.T .* not(obst_z);  
      
  % Increase coefficients close to obst (makes sense for momentum)
  sol_x = dif(X, obst_z);  % will be +1 north of obst, -1 south of obst
  corr = 1 + (sol_x > 0);  c.E(1:end-1,:,:) = c.E(1:end-1,:,:) .* corr;  
  corr = 1 + (sol_x < 0);  c.W(2:end,  :,:) = c.W(2:end,  :,:) .* corr;   
      
  sol_y = dif(Y, obst_z);  % will be +1 north of obst, -1 south of obst
  corr = 1 + (sol_y > 0);  c.N(:,1:end-1,:) = c.N(:,1:end-1,:) .* corr;  
  corr = 1 + (sol_y < 0);  c.S(:,2:end,  :) = c.S(:,2:end,  :) .* corr;   
end      
          
end