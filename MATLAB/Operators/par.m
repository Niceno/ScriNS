%==========================================================================
function [ y ] = par( mean_val, n )

max_val = mean_val * 3/2;

x = linspace(-1,1,n);

y = (1-x.*x) * max_val;


end

