# Standard modules used
from collections         import namedtuple
from math                import ceil,   \
                                floor,  \
                                sqrt
from matplotlib          import pyplot      as plt
from matplotlib          import cm
from numpy               import logical_not as lnot
from numpy               import minimum     as mn
from numpy               import maximum     as mx
from numpy               import array,        \
                                concatenate,  \
                                copy,         \
                                empty,        \
                                linspace,     \
                                matrix,       \
                                meshgrid,     \
                                ndarray,      \
                                ones,         \
                                outer,        \
                                prod,         \
                                reshape,      \
                                transpose,    \
                                zeros
from numpy.linalg        import solve
from random              import random
from scipy.sparse        import spdiags
from scipy.sparse.linalg import bicgstab
