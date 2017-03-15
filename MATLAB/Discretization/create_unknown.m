%==========================================================================
function [phi] = create_unkonwn( name, pos, res, def_bc )
%--------------------------------------------------------------------------
% This function creates a new unkown; helping to shorten the main program.
%--------------------------------------------------------------------------
% Input parameters are:
%
% name   - string holding the name of the variable 
%          it should be used for post-processing
% pos    - character specifying if the variable is cell centered (C),
%          staggered in x direction (X) or in y direction (Y)
% res    - vector specifying resolutions in x and y direction
% def_bc - character specifying if the default boundary condition is of 
%          dirichlet type ('d'), neumann type ('n') or outlet ('o')
%--------------------------------------------------------------------------

%---------------
% Save the name
%---------------
phi.name = name;

%------------------
% Set the position
%------------------
if pos ~= C && pos ~= X && pos ~= Y && pos ~= Z 
  error('Variable must be defined at positions ' + ...
        'C, X, Y or Z');
end
phi.pos = pos;

%-----------------------------
% Dimension and initial value
%-----------------------------
phi.val = zeros(res);
phi.old = zeros(res);

%---------------------------------
% Set default boundary conditions
%---------------------------------
if def_bc ~= 'd' && ...
   def_bc ~= 'n' && ...
   def_bc ~= 'o' 
  error('Variable must be defined with boundary conditions ' + ...
        '''d'', ''n'' or ''o''');
end
phi.bnd(W).type(1, 1:res(2), 1:res(3)) = def_bc;
phi.bnd(E).type(1, 1:res(2), 1:res(3)) = def_bc; 
phi.bnd(S).type(1:res(1), 1, 1:res(3)) = def_bc;
phi.bnd(N).type(1:res(1), 1, 1:res(3)) = def_bc;
phi.bnd(B).type(1:res(1), 1:res(2), 1) = def_bc;
phi.bnd(T).type(1:res(1), 1:res(2), 1) = def_bc;

phi.bnd(W).val(1, 1:res(2), 1:res(3)) = 0;
phi.bnd(E).val(1, 1:res(2), 1:res(3)) = 0; 
phi.bnd(S).val(1:res(1), 1, 1:res(3)) = 0;
phi.bnd(N).val(1:res(1), 1, 1:res(3)) = 0;
phi.bnd(B).val(1:res(1), 1:res(2), 1) = 0;
phi.bnd(T).val(1:res(1), 1:res(2), 1) = 0;

disp( sprintf('Created variable %s', name) );

end

