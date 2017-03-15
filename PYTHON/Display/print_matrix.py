#==========================================================================
def print_matrix(a):
#--------------------------------------------------------------------------

  print('Matrix['+('%d' %a.shape[0])+']['+('%d' %a.shape[1])+']')
  rows = a.shape[0]
  cols = a.shape[1]

  for i in range(0,rows):
    for j in range(0,cols):
      if a[i,j] == 0:
        print('    .  ', end='')
      else:
        print(('%7.2f' %a[i,j]), end='')
    print('')
  print('')      

  return  # end of function