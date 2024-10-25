import fem2d.constants as c
import numpy as np

# compliance matrix 2D

C = c.E/(1.0 + c.V)/(1.0 - 2.0*c.V) * np.array([[1.0-c.V,      c.V,   0.0],
                                                [c.V,  1.0-c.V,       0.0],
                                                [0.0,    0.0,     0.5-c.V]])
