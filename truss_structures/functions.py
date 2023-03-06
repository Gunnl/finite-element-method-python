"""
@author: Mad Dreamer
"""
import numpy as np
import math

"""
Assign boundary conditions, local Degrees of Freedom and global Degrees of Freedom, for truss structures
"""
def assign_BCs(NL, ENL):
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    DOFs = 0
    DOCs = 0

    # update local DoF
    for i in range(0,NoN):
        for j in range(0,PD):
            if ENL[i,PD+j] == -1:
                # it is a dirichelet node (-1)
                DOCs -= 1
                ENL[i,2*PD+j] = DOCs
            else:
                # it is a neuman node (1)
                DOFs += 1
                ENL[i,2*PD+j] = DOFs
    
    # update global DoF
    for i in range(0,NoN):
        for j in range(0,PD):
            if ENL[i,2*PD+j] < 0:
                # it is a dirichelet node (-1)
                ENL[i,3*PD+j] = abs(ENL[i,2*PD+j]) + DOFs
            else:
                # it is a neuman node (1)
                ENL[i,3*PD+j] = abs(ENL[i,2*PD+j])

    DOCs = abs(DOCs)

    return (ENL, DOFs, DOCs)

"""
Assemble stiffness, for truss structures
"""
def assemble_stiffness(ENL, EL, NL, E, A):
    # number of elements
    NoE = np.size(EL,0)
    # nodes per element
    NPE = np.size(EL,1)
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    K = np.zeros([NoN*PD,NoN*PD])

    for i in range(0,NoE):
        nl = EL[i,0:NPE]
        # element stiffness
        k = element_stiffness(nl,ENL,E,A)
        # insert element stiffness in global stiffness matrix
        for r in range(0, NPE):
            for p in range(0, PD):
                for q in range(0, NPE):
                    for s in range(0, PD):
                        row = ENL[nl[r]-1, p+3*PD]
                        column = ENL[nl[q]-1, s+3*PD]
                        value = k[r*PD+p, q*PD+s]
                        K[int(row)-1,int(column)-1] += value
    return K

"""
Element stiffness, for truss structures
"""
def element_stiffness(nl,ENL,E,A):
    X1 = ENL[nl[0]-1,0]
    Y1 = ENL[nl[0]-1,1]
    X2 = ENL[nl[1]-1,0]
    Y2 = ENL[nl[1]-1,1]

    # element lenght
    L = math.sqrt((X1-X2)**2+(Y1-Y2)**2)
    # cosine
    C = (X2-X1)/L
    # sine
    S = (Y2-Y1)/L

    #calculate element stiffness
    k = (E*A)/L * np.array([[C**2, C*S, -C**2, -C*S],
                            [C*S, S**2, -C*S, -S**2],
                            [-C**2, -C*S, C**2, C*S],
                            [-C*S, -S**2, C*S, S**2]])
    return k

"""
Assemble forces, for truss structures
"""
def assemble_forces(ENL, NL):
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    DOF = 0

    Fp = []

    for i in range(0, NoN):
        for j in range(0, PD):
            if ENL[i,PD+j] == 1:
                # it is a neuman node (1)
                DOF += 1
                Fp.append(ENL[i,5*PD+j])
    
    Fp = np.vstack([Fp]).reshape(-1,1)

    return Fp

"""
Assemble displacements, for truss structures
"""
def assemble_displacements(ENL, NL):
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    DOC = 0

    Up = []

    for i in range(0, NoN):
        for j in range(0, PD):
            if ENL[i,PD+j] == -1: # was 1 corrected? to -1
                # it is a neuman node (1)
                DOC += 1
                Up.append(ENL[i,4*PD+j])

    Up = np.vstack([Up]).reshape(-1,1)

    return Up


"""
Assemble displacements, for truss structures
"""
def update_nodes(ENL, U_u, NL, Fu):
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    DOFs = 0
    DOCs = 0
    
    for i in range(0,NoN):
        for j in range(0,PD):
            if ENL[i,PD+j] == 1:
                # it is a neuman node (1)
                DOFs += 1
                ENL[i,4*PD+j] = U_u[DOFs-1]
            else:
                # it is a dirichelet node (-1)
                DOCs += 1
                ENL[i,5*PD+j] = Fu[DOCs-1]

    return ENL
