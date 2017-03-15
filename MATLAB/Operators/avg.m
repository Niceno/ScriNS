%==========================================================================
function b = avg(d, a)

if nargin == 1
  b = (d(2:end) + d(1:end-1))/2; 
else
  if     d == X  b = (a(2:end,:,:) + a(1:end-1,:,:))*0.5; 
  elseif d == Y  b = (a(:,2:end,:) + a(:,1:end-1,:))*0.5;      
  elseif d == Z  b = (a(:,:,2:end) + a(:,:,1:end-1))*0.5;   
  else           b = a;    
  end    
end    
    
end