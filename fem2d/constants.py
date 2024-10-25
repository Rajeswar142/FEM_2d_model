import numpy as np
import math
# Constants that are further going to be used
# to create mesh
MESH_EX = 9             # number of elements in x direction
MESH_EY = 49            # number of elements in y direction
MESH_LX = 10            # length of shape in x direction
MESH_LY = 50            # length of shape in y direction

# derived
MESH_NX = MESH_EX + 1               # number of nodes in x direction
MESH_NY = MESH_EY + 1               # number of nodes in y direction
NUM_NODES = MESH_NX * MESH_NY       # total number of nodes in geometry
NUM_ELEMENTS = MESH_EX * MESH_EY    # total number of elements in geometry
MESH_HX = MESH_LX / MESH_EX         # element length in x direction
MESH_HY = MESH_LY / MESH_EY         # element length in y direction

# Material constants
E = 100
V = 0.48


