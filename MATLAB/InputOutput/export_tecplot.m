%==========================================================================
function export_tecplot(file_name, container, resolution, xn, yn, zn)
%--------------------------------------------------------------------------

nx = resolution(1);
ny = resolution(2);
nz = resolution(3);

%-------------- 
% write header
%--------------
file_id = fopen(file_name, 'w');
fprintf(file_id, 'title = "MatNaSt-FV Output"\n');
fprintf(file_id, 'variables = "x" "y" "z" ');
for n=1:size(container.var, 2)
  fprintf(file_id, '"%s" ', container.var(n).name);
end
fprintf(file_id, '\n');
fprintf(file_id, 'zone i=%0.0d, j=%0.0d k=%0.0d\n', nx+1, ny+1, nz+1);
fprintf(file_id, 'datapacking=block\n');
fprintf(file_id, 'varlocation=([%0.0d-%0.0d]=cellcentered)\n', ... 
                  4, size(container.var, 2) + 3);

%------------------------
% write node coordinates
%------------------------ 
fprintf(file_id, '# x coordinates\n');
c = 0;
for k=1:size(zn)
  for j=1:size(yn)
    for i=1:size(xn)
      fprintf(file_id, '%12.5e ', xn(i));
      c = c + 1;
      if mod(c,4)==0
        fprintf(file_id, '\n');
      end     
    end
  end
end  
fprintf(file_id, '\n');

fprintf(file_id, '# y-coordinates\n');
c = 0;
for k=1:size(zn)
  for j=1:size(yn)
    for i=1:size(xn)
      fprintf(file_id, '%12.5e ', yn(j));
      c = c + 1;
      if mod(c,4)==0
        fprintf(file_id, '\n');
      end
    end
  end
end  
fprintf(file_id, '\n');

fprintf(file_id, '# z-coordinates\n');
c = 0;
for k=1:size(zn)
  for j=1:size(yn)
    for i=1:size(xn)
      fprintf(file_id, '%12.5e ', zn(k));
      c = c + 1;
      if mod(c,4)==0
        fprintf(file_id, '\n');
      end
    end
  end
end  
fprintf(file_id, '\n');

%---------------------
% write all variables
%---------------------
for n=1:size(container.var, 2)
  
  % print a sub-header
  fprintf(file_id, '# %s\n', container.var(n).name);
  
  % arrange the variables according to their positions
  r = container.var(n);  
  values = r.val;
  if container.var(n).pos == X
    values = avg(X,cat(1, r.bnd(W).val, r.val, r.bnd(E).val));  
  end
  if container.var(n).pos == Y
    values = avg(Y,cat(2, r.bnd(S).val, r.val, r.bnd(N).val));  
  end
  if container.var(n).pos == Z
    values = avg(Z,cat(3, r.bnd(B).val, r.val, r.bnd(T).val));  
  end

  % print to file
  c = 0;
  for k=1:nz
    for j=1:ny
      for i=1:nx
        fprintf(file_id, '%12.5e ', values(i, j, k));
        c = c + 1;
        if mod(c,4)==0
          fprintf(file_id, '\n');
        end     
      end
    end
  end  
end

fclose(file_id);

end