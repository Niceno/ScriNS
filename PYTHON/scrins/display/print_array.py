#==========================================================================
def print_array(*args):
#--------------------------------------------------------------------------

  if len((args)) == 1:
    a      = args[0]
    format = "%7.2f"
  else:  
    a      = args[0]
    format = args[1]

  print('Array ['+('%d' %a.shape[0])+']')
  rows = a.shape[0]
  for i in range(0,rows):
    print(format %a[i])

  return  # end of function