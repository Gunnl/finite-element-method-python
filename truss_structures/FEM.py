"""
@author: Mad Dreamer
"""

# NoN = Number of Nodes
# NoE = Number of Elements
# NPE = Nodes per Element
# ENL = Extended Node List
# PD = Problem Dimension
# NL = Node List
# EL = Element List

import numpy as np
from functions import *

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

NL = np.array([[0,0],
              [1,0],
              [0.5,1]]
              )

EL = np.array([[1,2],
              [2,3],
              [3,1]]
              )

# Boundary conditions
# -1 = direchelet (dynamic)
# 1 = neuman (fixed)
# node 1 fixed in both x and y directions
DorN = np.array([[-1,-1],
                 [1,-1],
                 [1,1]]
                 )

# External forces
# force applied to node 3
Fu = np.array([[0,0],
               [0,0],
               [0,-20]])

U_u = np.array([[0,0],
                [0,0],
                [0,0]]
                )

# youngs modulus
E = 10**6
# cross section area
A = 0.01

# problem dimension
PD = np.size(NL,1)
# Number of Nodes
NoN = np.size(NL,0)
# ENL = Extended Node List
ENL = np.zeros([NoN,6*PD])
# setup first two columns (to Node List)
ENL[:,0:PD] = NL[:,:]
# setup columns 2 and 3 (to Node List)
ENL[:,PD:2*PD] = DorN[:,:]

(ENL, DOFs, DOCs) = assign_BCs(NL, ENL)

K = assemble_stiffness(ENL, EL, NL, E, A)

ENL[:,4*PD:5*PD] = U_u[:,:]
ENL[:,5*PD:6*PD] = Fu[:,:]

#print(U_u)
U_u = U_u.flatten()
#print(U_u)
#print(Fu)
Fu = Fu.flatten()
#print(Fu)

# calculate forces
Fp = assemble_forces(ENL, NL)

# calculate displacements
Up = assemble_displacements(ENL, NL)

#print(Up)
#print(Fp)

K_UU = K[0:DOFs, 0:DOFs]
K_UP = K[0:DOFs, DOFs:DOFs+DOCs]
K_PU = K[DOFs:DOFs+DOCs, 0:DOFs]
K_PP = K[DOFs:DOFs+DOCs, DOFs:DOFs+DOCs]

F = Fp - np.matmul(K_UP, Up)
U_u = np.matmul(np.linalg.inv(K_UU),F)
Fu = np.matmul(K_PU, U_u) + np.matmul(K_PP,Up)

# print(F)
#print(U_u)
#print(Fu)

#print(ENL)
ENL = update_nodes(ENL, U_u, NL, Fu)
#print(ENL)

scale = 100 # Exageration
coord = []
dispx_array = []

for i in range(np.size(NL,0)):
    dispx = ENL[i,8]
    dispy = ENL[i,9]

    x = ENL[i,0] + dispx * scale
    y = ENL[i,1] + dispy * scale

    dispx_array.append(dispx)
    coord.append(np.array([x,y]))

coord = np.vstack(coord)
dispx_array = np.vstack(dispx_array)

x_scatter = []
y_scatter = []
color_x = []

for i in range(0,np.size(EL,0)):
    x1 = coord[EL[i,0]-1,0]
    x2 = coord[EL[i,1]-1,0]
    y1 = coord[EL[i,0]-1,1]
    y2 = coord[EL[i,1]-1,1]

    dispx_EL = np.array([dispx_array[EL[i,0]-1]])
    dispx_EL = np.array([dispx_array[EL[i,0]-1],dispx_array[EL[i,1]-1]])

    if x1 == x2:
        x = np.linspace(x1,x2,200)
        y = np.linspace(y1,y2,200)
    else:
        m = (y2-y1)/(x2-x1)
        x = np.linspace(x1,x2,200)
        y = m*(x-x1)+y1
    
    x_scatter.append(x)
    y_scatter.append(y)

    color_x.append(np.linspace(np.abs(dispx_EL[0]),np.abs(dispx_EL[1]),200))

x_scatter = np.vstack([x_scatter]).flatten()
y_scatter = np.vstack([y_scatter]).flatten()
color_x = np.vstack([color_x]).flatten()

dispFigure = plt.figure(1)
ax_dispx = dispFigure.add_subplot(111)

cmap = plt.get_cmap('jet')
ax_dispx.scatter(x_scatter, y_scatter, c = color_x, cmap = cmap, s=10, edgecolor='none')

norm_x = Normalize(np.abs(dispx_array.min()),np.abs(dispx_array.max()))

dispFigure.colorbar(ScalarMappable(norm = norm_x, cmap = cmap))
plt.show()
