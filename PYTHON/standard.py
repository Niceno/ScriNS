# Standard modules used
from collections         import namedtuple
from math                import ceil,        \
                                floor,       \
                                factorial,   \
                                pi,          \
                                sqrt
from matplotlib          import pyplot      as plt
from matplotlib          import cm
from random              import random
from scipy               import logical_not as lnot
from scipy               import minimum     as mn
from scipy               import maximum     as mx
from scipy               import array,        \
                                append,       \
                                concatenate,  \
                                copy,         \
                                cos,          \
                                delete,       \
                                dot,          \
                                empty,        \
                                linspace,     \
                                log,          \
                                log2,         \
                                log10,        \
                                matrix,       \
                                meshgrid,     \
                                ndarray,      \
                                ones,         \
                                outer,        \
                                prod,         \
                                sin,          \
                                reshape,      \
                                tile,         \
                                transpose,    \
                                zeros
from scipy.linalg        import solve
from scipy.sparse        import spdiags
from scipy.sparse.linalg import bicgstab,     \
                                spsolve