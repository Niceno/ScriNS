%==========================================================================
function [ rho, mu, cp, lambda ] = water_properties( rx, ry, rz )
%--------------------------------------------------------------------------
% Returns physical properties of water in given ranges "rx", "ry" and "rz"
%
% For 60 deg from: 
%   http://www.engineeringtoolbox.com/water-properties-d_1508.html 
%--------------------------------------------------------------------------

rho   (rx, ry, rz) =  983;         % density              [kg/m^3]
mu    (rx, ry, rz) =    0.466E-3;  % viscosity            [Pa s]
cp    (rx, ry, rz) = 4185;         % thermal capacity     [J/kg/K]
lambda(rx, ry, rz) =    0.654;     % thermal conductivity [W/m/K]
  
end

