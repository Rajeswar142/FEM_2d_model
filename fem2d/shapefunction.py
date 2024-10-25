import numpy as np


def shape(xi):
    x, y = tuple(xi)
    n = [(1.0 - x)*(1.0 - y), (1.0 + x)*(1.0 - y), (1.0 + x)*(1.0 + y), (1.0 - x)*(1.0 + y)]
    return 0.25*np.array(n)


def grad_shape(xi):
    x, y = tuple(xi)
    dn = [[-(1.0-y), (1.0-y), (1.0+y), -(1.0+y)],
          [-(1.0-x), -(1.0+x), (1.0+x), (1.0-x)]]
    return 0.25*np.array(dn)
