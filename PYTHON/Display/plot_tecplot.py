# Standard Python modules
from standard import *

# ScriNS modules
from Constants.all      import *
from Operators.all      import *

#==============================================================================
def plot_tecplot(file_name, xyzn, variables):
#------------------------------------------------------------------------------
# Exports results to Tecplot (TM) format.
#------------------------------------------------------------------------------
  
  # Unpack tuples
  xn, yn, zn = xyzn    
  
  # Compute cell resolutions (remember: cell resolutions)
  nx = len(xn)-1
  ny = len(yn)-1
  nz = len(zn)-1
  
  file_id = open(file_name, 'w')

  #---------------------------
  # Write the file header out
  #---------------------------
  file_id.write("# File header \n")
  file_id.write("title=\"ScriNS Output\"\n")
  file_id.write("variables=\"x\" \"y\" \"z\" ")
  for v in variables:
    file_id.write("\"%s\" " % v.name)
  file_id.write("\n")
  file_id.write("zone i=%d" % (nx+1) + " j=%d" % (ny+1) + " k=%d\n" % (nz+1))
  file_id.write("datapacking = block\n")  
  file_id.write("varlocation=([1-3]=nodal ")
  file_id.write("[4-%d]=cellcentered)" % (3+len(variables)) )
  
  #--------------------------------------------------------------------
  # Write the coordinates out (remember - those are nodal coordinates)
  #--------------------------------------------------------------------  
  file_id.write("\n# X coordinates\n")
  c = 0                                    # column counter
  for k in range(0, nz+1):
    for j in range(0, ny+1):
      for i in range(0, nx+1):
        file_id.write("%12.5e " % xn[i])
        c = c + 1        
        if c % 4 == 0:                     # go to new line after 4th column
          file_id.write("\n")
          
  file_id.write("\n# Y coordinates\n")
  c = 0                                    # column counter
  for k in range(0, nz+1):
    for j in range(0, ny+1):
      for i in range(0, nx+1):
        file_id.write("%12.5e " % yn[j])
        c = c + 1        
        if c % 4 == 0:                     # go to new line after 4th column
          file_id.write("\n")
          
  file_id.write("\n# Z coordinates\n")
  c = 0                                    # column counter
  for k in range(0, nz+1):
    for j in range(0, ny+1):
      for i in range(0, nx+1):
        file_id.write("%12.5e " % zn[k])
        c = c + 1        
        if c % 4 == 0:                     # go to new line after 4th column
          file_id.write("\n")

  #-------------------------
  # Write the variables out 
  #-------------------------
  
  # Average values to be written for staggered variables
  for v in variables:
    if v.pos == C:
      val = v.val
    elif v.pos == X:
      val = avg(X,cat(X,(v.bnd[W].val[:1,:,:], v.val, v.bnd[E].val[:1,:,:])))
    elif v.pos == Y:
      val = avg(Y,cat(Y,(v.bnd[S].val[:,:1,:], v.val, v.bnd[N].val[:,:1,:])))
    elif v.pos == Z:
      val = avg(Z,cat(Z,(v.bnd[B].val[:,:,:1], v.val, v.bnd[T].val[:,:,:1])))
      
    file_id.write("\n# %s \n" % v.name)
    c = 0
    for k in range(0, nz):
      for j in range(0, ny):
        for i in range(0, nx):
          file_id.write("%12.5e " % val[i,j,k])
          c = c + 1        
          if c % 4 == 0:
            file_id.write("\n")

  file_id.close()

  return  # end of function