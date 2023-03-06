"""
@author: Mad Dreamer
"""
import numpy as np
#import matplotlib.pyplot as plt
#from mesh_structures.uniform_mesh import *
from functions import *
from void_mesh import *
from post_process import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from datetime import datetime
startTime = datetime.now()

d1 = 1
d2 = 1
p = 9
m = 9
R = 0.2
element_type = 'D2QU4N' # 2 dimension quadrilateral 4 node element
defV = 0.1

NL, EL = void_mesh(d1,d2,p,m,R,element_type)

BC_flag = 'extension'

(ENL, DOFs, DOCs) = assign_BCs(NL, BC_flag, defV)

K = assemble_stiffness(ENL, EL, NL)

Fp = assemble_forces(ENL, NL)
Up = assemble_displacements(ENL, NL)

#print(ENL)

# stiffness matrix
K_reduced = K[0:DOFs, 0:DOFs] # K_UU in 1D code
K_UP = K[0:DOFs, DOFs:DOCs+DOFs]
K_PU = K[DOFs:DOCs+DOFs, 0:DOFs]
K_PP = K[DOFs:DOCs+DOFs, DOFs:DOCs+DOFs]

F = Fp - (K_UP @ Up) # np.matmul()
Uu = np.linalg.solve(K_reduced,F)
Fu = (K_PU @ Uu) + (K_PP @ Up)

ENL = update_nodes(ENL, Uu, Fu, NL)

print(datetime.now() - startTime)

################################################

(stress_xx, stress_xy, stress_yx, stress_yy, strain_xx, strain_xy, strain_yx, strain_yy, disp_x, disp_y, X, Y) = post_process(NL,EL,ENL)

# Normalize color values for the colormap
stress_xxNormalized = (stress_xx - stress_xx.min()) / (stress_xx.max() - stress_xx.min())
disp_xNormalized = (disp_x - disp_x.min()) / (disp_x.max() - disp_x.min())

# Plot for each element

# import matplotlib.colors as mcolor
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval,maxval,n))
    )
    return new_cmap

fig_1 = plt.figure(1)
plt.title('Stress_XX')
axstress_xx = fig_1.add_subplot(111)

for i in range(np.size(EL,0)):
    x = X[:,i]
    y = Y[:,i]
    c = stress_xxNormalized[:,i]

    # Setting colormap boundaries
    cmap = truncate_colormap(plt.get_cmap('jet'), c.min(), c.max())

    # Plot the colors
    t = axstress_xx.tripcolor(x,y,c,cmap=cmap,shading='gouraud')
    # t = axstress_xx.tricontourf(x,y,c,cmap=cmap,levels=10)

    # plot the black lines
    p = axstress_xx.plot(x,y,'k-',linewidth=0.5)

fig_2 = plt.figure(2)
plt.title('Displacement X')
axdisp_x = fig_2.add_subplot(111)

for i in range(np.size(EL,0)):
    x = X[:,i]
    y = Y[:,i]
    c = disp_xNormalized[:,i]

    # Setting colormap boundaries
    cmap = truncate_colormap(plt.get_cmap('jet'), c.min(), c.max())

    # Plot the colors
    t = axdisp_x.tripcolor(x,y,c,cmap=cmap,shading='gouraud')
    # t = axstress_xx.tricontourf(x,y,c,cmap=cmap,levels=10)

    # plot the black lines
    p = axdisp_x.plot(x,y,'k-',linewidth=0.5)

print(datetime.now() - startTime)
plt.show()
