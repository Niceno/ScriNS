from numpy import concatenate

#==========================================================================
def cat(d, tup):
#--------------------------------------------------------------------------
# An interface to the standard NumPy's function concatenate 
#--------------------------------------------------------------------------

  return concatenate(tup, d)  # end of function