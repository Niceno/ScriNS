%==========================================================================
function [val] = obst_zero_source(pos, val, obst)
%--------------------------------------------------------------------------
% Set value to zero inside obstacle.
% 
% pos  - position of the variable, C, X, Y or Z
% val  - value to be zeroed in obstacle
% obst - matrix holding positions of obstacle
%--------------------------------------------------------------------------

if pos==C  
  val = val .* not(obst); 
  
elseif pos==X  
  obst_x = max(obst(1:end-1,:,:), obst(2:end,:,:));
  val = val .* not(obst_x); 
  
elseif pos==Y  
  obst_y = max(obst(:,1:end-1,:), obst(:,2:end,:));
  val = val .* not(obst_y); 
  
elseif pos==Z  
  obst_z = max(obst(:,:,1:end-1), obst(:,:,2:end));
  val = val .* not(obst_z); 
end

end