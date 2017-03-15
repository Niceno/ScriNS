%==========================================================================
function [dx, dy, dz, nx, ny, nz] = cartesian_grid(xn, yn, zn)
%--------------------------------------------------------------------------
% Spawns a 3D Cartesian grid from three arrays with node coordinates.
%--------------------------------------------------------------------------

dx = dif(xn);
dy = dif(yn);
dz = dif(zn);

nx = max( size(dx) );
ny = max( size(dy) );
nz = max( size(dz) );

dx = repmat(reshape(dx, nx, 1,  1),  1, ny, nz);
dy = repmat(reshape(dy, 1,  ny, 1),  nx, 1, nz);
dz = repmat(reshape(dz, 1,  1,  nz), nx, ny, 1);

end

