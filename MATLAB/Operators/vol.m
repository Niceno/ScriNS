%==========================================================================
function [ dv ] = vol( dx, dy, dz )
%--------------------------------------------------------------------------
% Returns a 3D array of computational volumes
%--------------------------------------------------------------------------

nx = max( size(dx) );
ny = max( size(dy) );
nz = max( size(dz) );

dv = repmat(reshape(dx, nx, 1,  1),  1, ny, nz) ...
  .* repmat(reshape(dy, 1,  ny, 1),  nx, 1, nz) ...
  .* repmat(reshape(dz, 1,  1,  nz), nx, ny, 1);

end

