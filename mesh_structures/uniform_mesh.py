"""
@author: Mad Dreamer
"""
import numpy as np

def uniform_mesh(d1,d2,p,m,element_type):
    PD = 2

    q = np.array([[0,0],[d1,0],[0,d2],[d1,d2]]) # 4 corners

    # Number of Nodes
    NoN = (p+1)*(m+1)
    # Numer of Elements
    NoE = p*m
    # Nodes per Element
    NPE = 4

    ## NODES ##

    NL = np.zeros([NoN, PD])

    # increment in horizontal direction
    a = (q[1,0]-q[0,0])/p
    # increment in vertical direction
    b = (q[2,1]-q[0,1])/m

    n = 0 # This will iterate the rows

    for i in range(1,m+2):
        for j in range(1,p+2):
            NL[n,0] = q[0,0] + (j-1)*a
            NL[n,1] = q[0,1] + (i-1)*b
            n += 1
    
    ## ELEMENTS ##

    EL = np.zeros([NoE,NPE], dtype=int)

    for i in range(1,m+1):
        for j in range(1,p+1):
            if j == 1:
                EL[(i-1)*p+j-1,0] = (i-1)*(p+1) + j
                EL[(i-1)*p+j-1,1] = EL[(i-1)*p+j-1,0] + 1
                EL[(i-1)*p+j-1,3] = EL[(i-1)*p+j-1,0] + (p+1)
                EL[(i-1)*p+j-1,2] = EL[(i-1)*p+j-1,3] + 1
            else:
                EL[(i-1)*p+j-1,0] = EL[(i-1)*p+j-2,1]
                EL[(i-1)*p+j-1,3] = EL[(i-1)*p+j-2,2]
                EL[(i-1)*p+j-1,1] = EL[(i-1)*p+j-1,0] + 1
                EL[(i-1)*p+j-1,2] = EL[(i-1)*p+j-1,3] + 1

    if element_type == 'D2TR3N':
        NPE_new = 3
        NoE_new = 2*NoE

        EL_new = np.zeros([NoE_new,NPE_new], dtype=int)

        for i in range(1, NoE+1):
            # first triangular element
            EL_new[2*(i-1), 0] = EL[i-1,0]
            EL_new[2*(i-1), 1] = EL[i-1,1]
            EL_new[2*(i-1), 2] = EL[i-1,2]
            # second triangular element
            EL_new[2*(i-1)+1, 0] = EL[i-1,0]
            EL_new[2*(i-1)+1, 1] = EL[i-1,2]
            EL_new[2*(i-1)+1, 2] = EL[i-1,3]

        EL = EL_new

    return NL, EL