import fem2d.constants as c
import numpy as np
import math
import fem2d.shapefunction as shape
import fem2d.material as mat
from matplotlib import pyplot as plt

nodes = []
# nodes in a geometry
for y in np.linspace(0.0, c.MESH_LY, c.MESH_NY):
    for x in np.linspace(0.0, c.MESH_LX, c.MESH_NX):
        nodes.append([x, y])
nodes = np.array(nodes)
print(nodes)

conn = []
for j in range(c.MESH_EY):
    for i in range(c.MESH_EX):
        n0 = i + j*c.MESH_NX
        conn.append([n0, n0+1, n0+1+c.MESH_NX, n0+c.MESH_NX])

# Create global stiffness matrix
K = np.zeros((2*c.NUM_NODES, 2*c.NUM_NODES))
q4 = [[x/math.sqrt(3.0), y/math.sqrt(3.0)] for y in [-1.0, 1.0] for x in [-1.0, 1.0]]
B = np.zeros((3, 8))
u_loc = np.zeros((8, 1))
stress = np.zeros((c.NUM_NODES, 3))
strain = np.zeros((c.NUM_NODES, 3))

for co in conn:
    xIe = nodes[co, :]
    ke = np.zeros((8, 8))
    for q in q4:
        dN = shape.grad_shape(q)
        J = np.dot(dN, xIe).T
        dN = np.dot(np.linalg.inv(J), dN)
        B[0, 0::2] = dN[0, :]
        B[1, 1::2] = dN[1, :]
        B[2, 0::2] = dN[1, :]
        B[2, 1::2] = dN[0, :]
        ke += np.dot(np.dot(B.T, mat.C), B) * np.linalg.det(J)
    for i, I in enumerate(co):
        for j, J in enumerate(co):
            K[2*I, 2*J] += ke[2*i, 2*j]
            K[2*I+1, 2*J] += ke[2*i+1, 2*j]
            K[2*I+1, 2*J+1] += ke[2*i+1, 2*j+1]
            K[2*I, 2*J+1] += ke[2*i, 2*j+1]

print('assign nodal forces and boundary conditions')
f = np.zeros((2*c.NUM_NODES))
for i in range(c.NUM_NODES):
    if nodes[i, 1] == 0.0:
        K[2*i, :] = 0.0
        K[2*i+1, :] = 0.0
        K[2*i, 2*i] = 1.0
        K[2*i+1, 2*i+1] = 1.0
    if nodes[i, 1] == c.MESH_LY:
        x = nodes[i, 0]
        f[2*i+1] = 20.0
        if x == 0.0 or x == c.MESH_LX:
            f[2*i+1] *= 0.5

print('Solving linear system')
u = np.linalg.solve(K, f)
print('max u=', max(u))

print('Stress calculation')
prop_xx = []
prop_yy = []
prop_xy = []

for co in conn:
    xIe = nodes[co, :]
    ke = np.zeros((8, 8))
    for q in q4:
        dN = shape.grad_shape(q)
        J = np.dot(dN, xIe).T
        dN = np.dot(np.linalg.inv(J), dN)
        B[0, 0::2] = dN[0, :]
        B[1, 1::2] = dN[1, :]
        B[2, 0::2] = dN[1, :]
        B[2, 1::2] = dN[0, :]

    for i, I in enumerate(co):
        u_loc[2*i] = u[2*I]
        u_loc[2*i+1] = u[2*I+1]

    strain_loc = np.dot(B, u_loc)
    stress_loc = np.dot(mat.C, strain_loc)
    prop_xx.append([strain_loc[0], stress_loc[0]])
    prop_yy.append([strain_loc[1], stress_loc[1]])
    prop_xy.append([strain_loc[2], stress_loc[2]])

print('plotting displacement')
ux = np.reshape(u[0::2], (c.MESH_NY, c.MESH_NX))
uy = np.reshape(u[1::2], (c.MESH_NY, c.MESH_NX))
x_vec = []
y_vec = []
res = []
for i in range(c.MESH_NX):
    for j in range(c.MESH_NY):
        x_vec.append(i*c.MESH_HX + ux[j, i])
        y_vec.append(j*c.MESH_HY + uy[j, i])
        res.append(uy[j, i])
t = plt.tricontourf(x_vec, y_vec, res, levels=14, cmap=plt.cm.jet)
plt.scatter(x_vec, y_vec, marker='o', c='b', s=2)
plt.grid()
plt.colorbar(t)
plt.axis('equal')
plt.show()

# print('plotting stress - xx')
# sigma_x = np.reshape(prop_xx[0::2], (c.MESH_NY, c.MESH_NX))
# uy = np.reshape(u[1::2], (c.MESH_NY, c.MESH_NX))
# x_vec = []
# y_vec = []
# res = []
# for i in range(c.MESH_NX):
#     for j in range(c.MESH_NY):
#         x_vec.append(i*c.MESH_HX + ux[j, i])
#         y_vec.append(j*c.MESH_HY + uy[j, i])
#         res.append(uy[j, i])
# t = plt.tricontourf(x_vec, y_vec, res, levels=14, cmap=plt.cm.jet)
# plt.scatter(x_vec, y_vec, marker='o', c='b', s=2)
# plt.grid()
# plt.colorbar(t)
# plt.axis('equal')
# plt.show()
