%==========================================================================
function b = dif(d, a)

if nargin == 1
  b = (d(2:end) - d(1:end-1)); 
else
  if     d == X  b = (a(2:end,:,:) - a(1:end-1,:,:)); 
  elseif d == Y  b = (a(:,2:end,:) - a(:,1:end-1,:));      
  elseif d == Z  b = (a(:,:,2:end) - a(:,:,1:end-1));   
  end    
end    


end