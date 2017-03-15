%==========================================================================
function [src] = vol_balance(uf, vf, wf, dx, dy, dz, obst)
%--------------------------------------------------------------------------
% Computes the volume balance, which is essentially a right hand side in
% the Poisson's equation for pressure.
%
% Note that "obst" is an optional parameter.  If it is not sent, the source
% won't be zeroed inside the obstacle.  That is important for calculation
% of pressure, see function "calc_p".
%--------------------------------------------------------------------------

% Compute it throughout the domain
src = - dif(X, cat(X, uf.bnd(W).val, uf.val, uf.bnd(E).val)).*dy.*dz ...
      - dif(Y, cat(Y, vf.bnd(S).val, vf.val, vf.bnd(N).val)).*dx.*dz ...   
      - dif(Z, cat(Z, wf.bnd(B).val, wf.val, wf.bnd(T).val)).*dx.*dy;   

% Zero it inside obstacles, if obstacle is sent as parameter
if nargin == 7
  if size(obst,1) ~= 0
    src = obst_zero_val(C, src, obst);
  end
end    
  
end  
  