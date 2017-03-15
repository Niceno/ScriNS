%==========================================================================
function [x] = nodes(s, e, n, del_1, del_e)
%--------------------------------------------------------------------------
% Creates node coordinates in one dimension
%--------------------------------------------------------------------------

%---------------------------------------------
%
% Just a simple constant-spacing distribution
% 
%---------------------------------------------
if nargin == 3
  x = linspace(s, e, n+1);

else
%------------------------------------------------------------------------
%
% Grid distribution is determined from the following function:
%   x = a + b y + c y^2 + d y^3
%
% One should imagine this function as taking integer arguments.
%
% Boundary conditions:
%
%  1. x(1)   = s          =>  a + b + c + d = s
%
%  2. x(2)   = s + del_1  =>  a + b*2 + c*4 + d*8 = s + del_1
%
%  3. x(n)   = e - del_e  =>  a + b*n + c*n^2 + d*n^3 = e - del_e
%
%  4. x(n+1) = e          =>  a + b*(n+1) + c*(n+1)^2 + d*(n+1)^3 = e
%
%  =>
%
%  |1   1      1        1     |  |a|  |s      |
%  |1   2      4        8     |  |b|  |s+del_1|
%  |1   n      n^2      n^3   |  |c|  |e-del_e|
%  |1  (n+1)  (n+1)^2  (n+1)^3|  |d|  |e      |
%
%------------------------------------------------------------------------
  A = [1  1     1       1;    ...
       1  2     4       8;    ...
       1  n     n^2     n^3;  ...
       1 (n+1) (n+1)^2 (n+1)^3];
 
  f = [s;        ...
       s+del_1;  ...
       e-del_e;  ...
       e];  
   
  abcd = A \ f;
  
  for i=1:n+1
    x(i) = abcd(1) + abcd(2) * i + abcd(3) * i^2 + abcd(4) * i^3; 
  end    
  
end

