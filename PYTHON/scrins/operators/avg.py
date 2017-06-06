"""
Docstring.
"""
from scrins.constants.coordinates import X, Y, Z
#==========================================================================
def avg(*args):
    """
    Docstring.
    """

    # Only one argument is sent - perform
    # averaging of one-dimensional array
    if len((args)) == 1:
        a = args[0]
        return (a[1:] + a[:-1]) * 0.5

    # Two arguments are sent - perform averaging of a
    # three-dimensional array, with a given direction
    elif len((args)) == 2:
        d = args[0]  # direction
        a = args[1]  # array

        if d == X:
            return (a[1:, :, : ] + a[:-1, :, :]) * 0.5
        elif d == Y:
            return (a[:, 1:, : ] + a[:, :-1, :]) * 0.5
        elif d == Z:
            return (a[:, :, 1:] + a[:, :, :-1]) * 0.5

    # Some error message might
    # be printed if number of
    # arguments is wrong

    return a  # end of function
