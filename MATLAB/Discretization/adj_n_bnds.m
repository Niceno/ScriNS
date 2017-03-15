%==========================================================================
function phi = adj_n_bnds(phi)
%--------------------------------------------------------------------------
% Copies last domain cell values to Neumann boundary condition values.
%--------------------------------------------------------------------------

% These arrays will hold values true (1) in cells with boundary ... 
% ... condition of Neumann type, and false (0) otherwise
if_w_n = ( phi.bnd(W).type(1,:,:) == 'n' );  % 1 if west is Neumann type
if_e_n = ( phi.bnd(E).type(1,:,:) == 'n' );  % 1 if east is Neumann type
if_s_n = ( phi.bnd(S).type(:,1,:) == 'n' );  % 1 if south is Neumann type
if_n_n = ( phi.bnd(N).type(:,1,:) == 'n' );  % 1 if north is Neumann type
if_b_n = ( phi.bnd(B).type(:,:,1) == 'n' );  % 1 if bottom is Neumann type
if_t_n = ( phi.bnd(T).type(:,:,1) == 'n' );  % 1 if top is Neumann type

% In what follows, a linear combination of true (1) and false (0) 
% will copy the values of variable phi to the boundaries.
phi.bnd(W).val(1,:,:) = phi.bnd(W).val(1,:,:) .* ( not(if_w_n) ) + ...
                        phi.val(1,:,:)        .*       if_w_n;
phi.bnd(E).val(1,:,:) = phi.bnd(E).val(1,:,:) .* ( not(if_e_n) ) + ...
                        phi.val(end,:,:)      .*       if_e_n;

phi.bnd(S).val(:,1,:) = phi.bnd(S).val(:,1,:) .* ( not(if_s_n) ) + ...
                        phi.val(:,1,:)        .*       if_s_n;
phi.bnd(N).val(:,1,:) = phi.bnd(N).val(:,1,:) .* ( not(if_n_n) ) + ...
                        phi.val(:,end,:)      .*       if_n_n;

phi.bnd(B).val(:,:,1) = phi.bnd(B).val(:,:,1) .* ( not(if_b_n) ) + ...
                        phi.val(:,:,1)        .*       if_b_n;
phi.bnd(T).val(:,:,1) = phi.bnd(T).val(:,:,1) .* ( not(if_t_n) ) + ...
                        phi.val(:,:,end)      .*       if_t_n;

end  