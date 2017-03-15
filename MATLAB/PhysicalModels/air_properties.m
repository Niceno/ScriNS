%==========================================================================
function [ rho, mu, cp, lambda ] = air_properties( rx, ry, rz )
%--------------------------------------------------------------------------
% Returns physical properties of air in given ranges "rx", "ry" and "rz"
%
% For 60 deg from: 
%   http://www.engineeringtoolbox.com/air-properties-d_156.html 
%--------------------------------------------------------------------------

rho   (rx, ry, rz) =    1.067;     % density              [kg/m^3]
mu    (rx, ry, rz) =   20.17E-06;  % viscosity            [Pa s]
cp    (rx, ry, rz) = 1009;         % thermal capacity     [J/kg/K]
lambda(rx, ry, rz) =    0.0285;    % thermal conductivity [W/m/K]
  
end

