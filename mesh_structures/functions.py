"""
@author: Mad Dreamer
"""
import numpy as np
import math

"""
Assign boundary conditions, local Degrees of Freedom and global Degrees of Freedom, for mesh structures
"""
def assign_BCs(NL, BC_flag, defV):
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    ENL = np.zeros([NoN,6*PD])

    ENL[:,0:PD] = NL

    if BC_flag == 'extension':
        for i in range(0,NoN):
            if ENL[i,0] == 0:
                # dirichelet node boundary condition
                ENL[i,2] = -1
                ENL[i,3] = -1
                # apply displacements
                ENL[i,8] = -defV # x axis
                ENL[i,9] = 0 # y axis
            elif ENL[i,0] == 1:
                # dirichelet node boundary condition
                ENL[i,2] = -1
                ENL[i,3] = -1
                # apply displacements
                ENL[i,8] = defV # x axis
                ENL[i,9] = 0 # y axis
            else:
                # neuman node boundary condition
                ENL[i,2] = 1
                ENL[i,3] = 1
                # apply forces
                ENL[i,10] = 0 # x axis
                ENL[i,11] = 0 # y axis

    elif BC_flag == 'expansion':
        for i in range(0,NoN):
            if ENL[i,0] == 0 or ENL[i,0] == 1 or ENL[i,1] == 0 or ENL[i,1] == 1:
                # dirichelet node boundary condition
                ENL[i,2] = -1
                ENL[i,3] = -1
                # apply displacements
                ENL[i,8] = defV * ENL[i,0] # x axis
                ENL[i,9] = defV * ENL[i,1] # y axis
            else:
                # neuman node boundary condition
                ENL[i,2] = 1
                ENL[i,3] = 1
                # apply forces
                ENL[i,10] = 0 # x axis
                ENL[i,11] = 0 # y axis
    elif BC_flag == 'shear':
        for i in range(0,NoN):
            if ENL[i,1] == 0:
                # dirichelet node boundary condition
                ENL[i,2] = -1
                ENL[i,3] = -1
                # apply displacements
                ENL[i,8] = 0 # x axis
                ENL[i,9] = 0 # y axis
            if ENL[i,1] == 1:
                # dirichelet node boundary condition
                ENL[i,2] = -1
                ENL[i,3] = -1
                # apply displacements
                ENL[i,8] = defV # x axis
                ENL[i,9] = 0 # y axis
            else:
                # neuman node boundary condition
                ENL[i,2] = 1
                ENL[i,3] = 1
                # apply forces
                ENL[i,10] = 0 # x axis
                ENL[i,11] = 0 # y axis

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
Assemble stiffness, for mesh structures
"""
def assemble_stiffness(ENL, EL, NL):
    # number of elements
    NoE = np.size(EL,0)
    # nodes per element
    NPE = np.size(EL,1)
    # problem dimension
    PD = np.size(NL,1)
    # Number of Nodes
    NoN = np.size(NL,0)

    K = np.zeros([NoN*PD,NoN*PD])

    for i in range(1,NoE+1):
        nl = EL[i-1,0:NPE]
        # element stiffness
        k = element_stiffness(nl,NL)
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
    # Number of Nodes
    NoN = np.size(NL,0)
    # problem dimension
    PD = np.size(NL,1)

    DOC = 0

    Up = []

    for i in range(0, NoN):
        for j in range(0, PD):
            if ENL[i,PD+j] == -1:
                # it is a neuman node (1)
                DOC += 1
                Up.append(ENL[i,4*PD+j])

    Up = np.vstack([Up]).reshape(-1,1)

    return Up

"""
Element stiffness, for truss structures
"""
def element_stiffness(nl,NL):
    # nodes per element
    NPE = np.size(nl,0)
    # problem dimension
    PD = np.size(NL,1)

    x = np.zeros([NPE,PD])
    x[0:NPE, 0:PD] = NL[nl[0:NPE]-1,0:PD]

    K = np.zeros([NPE*PD,NPE*PD])

    coord = x.T

    if NPE == 3:
        GPE = 1
    elif NPE == 4:
        GPE = 4
    else:
        NotImplemented()

    for i in range(1, NPE+1):
        for j in range(1, NPE+1):
            
            k = np.zeros([PD, PD])

            for gp in range(1, GPE+1):
                
                J = np.zeros([PD, PD])

                grad = np.zeros([PD, NPE])

                (xi, eta, alpha) = gauss_point(NPE, GPE, gp)

                grad_nat = grad_N_nat(NPE, xi, eta)

                # jacobian
                J = coord @ grad_nat.T # np.matmul(,)

                grad = np.linalg.inv(J).T @ grad_nat

                for a in range(1, PD+1):
                    for c in range(1, PD+1):
                        for b in range(1, PD+1):
                            for d in range(1, PD+1):
                                if NPE == 4:
                                    k[a-1, c-1] = k[a-1, c-1] + grad[b-1,i-1] * constitutive(
                                        a,b,c,d) * grad[d-1,j-1] * np.linalg.det(J) * alpha
                                elif NPE == 3:
                                    k[a-1, c-1] = k[a-1, c-1] + grad[b-1,i-1] * constitutive(
                                        a,b,c,d) * grad[d-1,j-1] * np.linalg.det(J) * alpha * 1/2
                                else:
                                    NotImplemented()
                
                K[((i-1)*PD+1)-1:i*PD , ((j-1)*PD+1)-1:j*PD] = k
    
    return K

"""
Calculate stiffness, for mesh structures
"""
def stiffness(x, GPE):
    # nodes per element
    NPE = np.size(x,0)
    # problem dimension
    PD = np.size(x,1)

    K = np.zeros([NPE*PD,NPE*PD])

    coord = x.T

    for i in range(1, NPE+1):
        for j in range(1, NPE+1):
            
            k = np.zeros([PD, PD])

            for gp in range(1, GPE+1):
                
                J = np.zeros([PD, PD])

                grad = np.zeros([PD, NPE])

                (xi, eta, alpha) = gauss_point(NPE, GPE, gp)

                grad_nat = grad_N_nat(NPE, xi, eta)

                # jacobian
                J = coord @ grad_nat.T # np.matmul(,)

                grad = np.linalg.inv(J).T @ grad_nat

                for a in range(1, PD+1):
                    for c in range(1, PD+1):
                        for b in range(1, PD+1):
                            for d in range(1, PD+1):
                                if NPE == 4:
                                    k[a-1, c-1] = k[a-1, c-1] + grad[b-1,i-1] * constitutive(
                                        a,b,c,d) * grad[d-1,j-1] * np.linalg.det(J) * alpha
                                elif NPE == 3:
                                    k[a-1, c-1] = k[a-1, c-1] + grad[b-1,i-1] * constitutive(
                                        a,b,c,d) * grad[d-1,j-1] * np.linalg.det(J) * alpha * 1/2
                
                K[((i-1)*PD+1)-1:i*PD , ((j-1)*PD+1)-1:j*PD] = k
    
    return K

"""
Calculate gaus point, for mesh structures
"""
def gauss_point(NPE, GPE, gp):
    if NPE == 4:
        # D2QU4N case
        if GPE == 1:
            if gp == 1:
                xi = 0
                eta = 0
                alpha = 4
        elif GPE == 4:
            if gp == 1:
                xi = -1 / math.sqrt(3)
                eta = -1 / math.sqrt(3)
                alpha = 1
            elif gp == 2:
                xi = 1 / math.sqrt(3)
                eta = -1 / math.sqrt(3)
                alpha = 1
            elif gp == 3:
                xi = 1 / math.sqrt(3)
                eta = 1 / math.sqrt(3)
                alpha = 1
            elif gp == 4:
                xi = -1 / math.sqrt(3)
                eta = 1 / math.sqrt(3)
                alpha = 1
    elif NPE == 3:
        # D2TR3N case
        if (gp == 1):
            xi = 1/3
            eta = 1/3
            alpha = 1
    
    return (xi, eta, alpha)

"""
Calculate grad N nat, for mesh structures
"""
def grad_N_nat(NPE, xi, eta):
    PD = 2
    result = np.zeros([PD,NPE])

    if NPE == 3:
        # D2TR3N case
        result[0,0] = 1
        result[0,1] = 0
        result[0,2] = -1

        result[1,0] = 0
        result[1,1] = 1
        result[1,2] = -1
    elif NPE == 4:
        # D2QU4N case
        result[0,0] = -1/4*(1-eta)
        result[0,1] = 1/4*(1-eta)
        result[0,2] = 1/4*(1+eta)
        result[0,3] = -1/4*(1+eta)

        result[1,0] = -1/4*(1-xi)
        result[1,1] = -1/4*(1+xi)
        result[1,2] = 1/4*(1+xi)
        result[1,3] = 1/4*(1-xi)
    
    return result

"""
Calculate constitutive, for mesh structures
"""
def constitutive(i,j,k,l):
    E = 8/3
    nu = 1/3

    C = (E/(2*(1+nu))) * (delta(i,l)*delta(j,k) + delta(i,k)*delta(j,l)) + (E*nu)/(1-nu**2) * delta(i,j)*delta(k,l)

    return C

"""
Calculate delta, for mesh structures
"""
def delta(i,j):
    if i == j:
        delta = 1
    else:
        delta = 0
    return delta

"""
Assemble displacements, for truss structures
"""
def update_nodes(ENL, U_u, Fu, NL):
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
