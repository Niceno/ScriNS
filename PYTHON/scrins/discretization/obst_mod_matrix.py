
from scipy import logical_not as lnot
from scipy import maximum as mx

from scrins.constants.coordinates import X, Y, Z
from scrins.constants.compass import W, E, S, N, B, T, C
from scrins.operators.dif import dif

#==========================================================================
def obst_mod_matrix(phi, c, obst, obc):
#--------------------------------------------------------------------------
# Adjusts the system matrix for obstacles and cell centered varaibles 
# (such as pressure)
# 
# phi  - variable
# c    - coefficients in system matrix
# obst - obstacle array
# obc  - obstacles's boundary condition, ('n' - Neumann, 'd' - Dirichlet)
#--------------------------------------------------------------------------

  pos = phi.pos
  
  #--------------------------
  #
  # For collocated variables
  #
  #--------------------------
  if pos == C:
  
    #------------------------------------
    # Neumann's boundary on the obstacle
    #------------------------------------
    if obc == 'n':
        
      # Correct west and east
      sol_x = dif(X, obst)  # will be +1 east of obst, -1 west of obst
      corr = 1 - (sol_x < 0)  
      c.W[1:,:,:] = c.W[1:,:,:] * corr
      corr = 1 - (sol_x > 0)  
      c.E[:-1,:,:] = c.E[:-1,:,:] * corr
   
      # Correct south and north
      sol_y = dif(Y, obst)  # will be +1 north of obst, -1 south of obst  
      corr = 1 - (sol_y < 0)  
      c.S[:,1:,:] = c.S[:,1:,:] * corr
      corr = 1 - (sol_y > 0)  
      c.N[:,:-1,:] = c.N[:,:-1,:] * corr
    
      # Correct bottom and top
      sol_z = dif(Z, obst)  # will be +1 north of obst, -1 south of obst  
      corr = 1 - (sol_z < 0)  
      c.B[:,:,1:] = c.B[:,:,1:] * corr
      corr = 1 - (sol_z > 0)  
      c.T[:,:,:-1] = c.T[:,:,:-1] * corr
  
    #--------------------------------------
    # Dirichlet's boundary on the obstacle
    #--------------------------------------
    elif obc == 'd':
    
      # Set central coefficient to 1 in obst, unchanged elsewhere
      c.P[:] = c.P[:] * lnot(obst) + obst 
  
      # Set neighbour coefficients to zero in obst  
      c.W[:] = c.W[:] * lnot(obst)
      c.E[:] = c.E[:] * lnot(obst)
      c.S[:] = c.S[:] * lnot(obst)
      c.N[:] = c.N[:] * lnot(obst)
      c.B[:] = c.B[:] * lnot(obst)
      c.T[:] = c.T[:] * lnot(obst)
        
      # Increase coefficients close to obst (makes sense for momentum)
      sol_x = dif(X, obst)  # will be +1 east of obst, -1 west of obst
      corr = 1 + (sol_x > 0)  
      c.E[:-1,:,:] = c.E[:-1,:,:] * corr  
      corr = 1 + (sol_x < 0)  
      c.W[1:,:,:] = c.W[1:,:,:] * corr  
        
      sol_y = dif(Y, obst)  # will be +1 north of obst, -1 south of obst
      corr = 1 + (sol_y > 0)  
      c.N[:,:-1,:] = c.N[:,:-1,:] * corr  
      corr = 1 + (sol_y < 0)  
      c.S[:,1:,:] = c.S[:,1:,:] * corr  
        
      sol_z = dif(Z, obst)  # will be +1 top of obst, -1 bottom of obst
      corr = 1 + (sol_z > 0)  
      c.T[:,:,:-1] = c.T[:,:,:-1] * corr  
      corr = 1 + (sol_z < 0)  
      c.B[:,:,1:] = c.B[:,:,1:] * corr  
  
  #-------------------------
  #
  # For staggered variables
  #
  #-------------------------
  elif pos == X: 
      
    # Set central coefficient to 1 in obst, unchanged elsewhere
    obst_x = mx(obst[:-1,:,:], obst[1:,:,:])
    c.P[:] = c.P[:] * lnot(obst_x) + obst_x 
  
    # Set neighbour coefficients to zero in obst  
    c.W[:] = c.W[:] * lnot(obst_x)
    c.E[:] = c.E[:] * lnot(obst_x)
    c.S[:] = c.S[:] * lnot(obst_x)
    c.N[:] = c.N[:] * lnot(obst_x)
    c.B[:] = c.B[:] * lnot(obst_x)
    c.T[:] = c.T[:] * lnot(obst_x)
        
    # Increase coefficients close to obst (makes sense for momentum)
    sol_y = dif(Y, obst_x)  # will be +1 north of obst, -1 south of obst
    corr = 1 + (sol_y > 0)  
    c.N[:,:-1,:] = c.N[:,:-1,:] * corr  
    corr = 1 + (sol_y < 0)  
    c.S[:,1:,:] = c.S[:,1:,:] * corr  
       
    sol_z = dif(Z, obst_x)  # will be +1 top of obst, -1 bottom of obst
    corr = 1 + (sol_z > 0)  
    c.T[:,:,:-1] = c.T[:,:,:-1] * corr  
    corr = 1 + (sol_z < 0)  
    c.B[:,:,1:] = c.B[:,:,1:] * corr  
        
  elif pos == Y:  
          
    # Set central coefficient to 1 in obst, unchanged elsewhere
    obst_y = mx(obst[:,:-1,:], obst[:,1:,:])
    c.P[:] = c.P[:] * lnot(obst_y) + obst_y 
  
    # Set neighbour coefficients to zero in obst  
    c.W[:] = c.W[:] * lnot(obst_y)
    c.E[:] = c.E[:] * lnot(obst_y)
    c.S[:] = c.S[:] * lnot(obst_y)
    c.N[:] = c.N[:] * lnot(obst_y)  
    c.B[:] = c.B[:] * lnot(obst_y)
    c.T[:] = c.T[:] * lnot(obst_y)  
        
    # Increase coefficients close to obst (makes sense for momentum)
    sol_x = dif(X, obst_y)  # will be +1 north of obst, -1 south of obst
    corr = 1 + (sol_x > 0)  
    c.E[:-1,:,:] = c.E[:-1,:,:] * corr  
    corr = 1 + (sol_x < 0)  
    c.W[1:,:,:] = c.W[1:,:,:] * corr   
        
    sol_z = dif(Z, obst_y)  # will be +1 north of obst, -1 south of obst
    corr = 1 + (sol_z > 0)  
    c.T[:,:,:-1] = c.T[:,:,:-1] * corr  
    corr = 1 + (sol_z < 0)  
    c.B[:,:,1:] = c.B[:,:,1:] * corr   
        
  elif pos == Z:
          
    # Set central coefficient to 1 in obst, unchanged elsewhere
    obst_z = mx(obst[:,:,:-1], obst[:,:,1:])
    c.P[:] = c.P[:] * lnot(obst_z) + obst_z 
  
    # Set neighbour coefficients to zero in obst  
    c.W[:] = c.W[:] * lnot(obst_z)
    c.E[:] = c.E[:] * lnot(obst_z)
    c.S[:] = c.S[:] * lnot(obst_z)
    c.N[:] = c.N[:] * lnot(obst_z)  
    c.B[:] = c.B[:] * lnot(obst_z)
    c.T[:] = c.T[:] * lnot(obst_z)  
        
    # Increase coefficients close to obst (makes sense for momentum)
    sol_x = dif(X, obst_z)  # will be +1 north of obst, -1 south of obst
    corr = 1 + (sol_x > 0)  
    c.E[:-1,:,:] = c.E[:-1,:,:] * corr  
    corr = 1 + (sol_x < 0)  
    c.W[1:,:,:] = c.W[1:,:,:] * corr   
        
    sol_y = dif(Y, obst_z)  # will be +1 north of obst, -1 south of obst
    corr = 1 + (sol_y > 0)  
    c.N[:,:-1,:] = c.N[:,:-1,:] * corr  
    corr = 1 + (sol_y < 0)  
    c.S[:,1:,:] = c.S[:,1:,:] * corr   

  return c  # end of function