#==========================================================================
def print_matrix(*args):
#--------------------------------------------------------------------------

  if len((args)) == 1:
    a      = args[0]
    format = "%7.2f"
  else:  
    a      = args[0]
    format = args[1]

  print('Matrix['+('%d' %a.shape[0])+']['+('%d' %a.shape[1])+']')
  rows = a.shape[0]
  cols = a.shape[1]

  for i in range(0,rows):
    for j in range(0,cols):
      print((format %a[i,j]), end='')
    print('')
  print('')      

  return  # end of function