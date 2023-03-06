"""
@author: Mad Dreamer
"""
import numpy as np
import math 

def void_mesh(d1,d2,p,m,R,element_type):
    PD = 2

    q = np.array([[0,0],[d1,0],[0,d2],[d1,d2]]) # 4 corners

    # Number of Nodes
    NoN = 2*(p+1)*(m+1)+2*(p-1)*(m+1)
    # Numer of Elements
    NoE = 4*p*m
    # Nodes per Element
    NPE = 4

    ## NODES ##

    NL = np.zeros([NoN, PD])

    # increment in horizontal direction
    a = (q[1,0]-q[0,0])/p
    # increment in vertical direction
    b = (q[2,1]-q[0,1])/p

    n = 0 # This will iterate the rows

    ## region 1
    coord11 = np.zeros([(p+1)*(m+1),PD])

    for i in range(1,p+2):
        coord11[i-1,0] = q[0,0] + (i-1)*a
        coord11[i-1,1] = q[0,1]
    
    for i in range(1,p+2):
        coord11[m*(p+1)+i-1,0] = R*np.cos( (5*math.pi/4) + (i-1)*((math.pi/2)/p)) + d1/2
        coord11[m*(p+1)+i-1,1] = R*np.sin( (5*math.pi/4) + (i-1)*((math.pi/2)/p)) + d2/2
    
    for i in range(1, m):
        for j in range(1, p+2):
            dx = (coord11[m*(p+1)+j-1,0] - coord11[j-1,0])/m
            dy = (coord11[m*(p+1)+j-1,1] - coord11[j-1,1])/m

            coord11[i*(p+1)+j-1,0] = coord11[(i-1)*(p+1)+j-1,0] + dx
            coord11[i*(p+1)+j-1,1] = coord11[(i-1)*(p+1)+j-1,1] + dy
    
    ## region 2
    coord22 = np.zeros([(p+1)*(m+1),PD])

    for i in range(1,p+2):
        coord22[i-1,0] = q[2,0] + (i-1)*a
        coord22[i-1,1] = q[2,1]
    
    for i in range(1,p+2):
        coord22[m*(p+1)+i-1,0] = R*np.cos( (3*math.pi/4) - (i-1)*((math.pi/2)/p)) + d1/2
        coord22[m*(p+1)+i-1,1] = R*np.sin( (3*math.pi/4) - (i-1)*((math.pi/2)/p)) + d2/2
    
    for i in range(1, m):
        for j in range(1, p+2):
            dx = (coord22[m*(p+1)+j-1,0] - coord22[j-1,0])/m
            dy = (coord22[m*(p+1)+j-1,1] - coord22[j-1,1])/m

            coord22[i*(p+1)+j-1,0] = coord22[(i-1)*(p+1)+j-1,0] + dx
            coord22[i*(p+1)+j-1,1] = coord22[(i-1)*(p+1)+j-1,1] + dy
    
    ## region 3
    coord33 = np.zeros([(p-1)*(m+1),PD])

    for i in range(1,p):
        coord33[i-1,0] = q[0,0]
        coord33[i-1,1] = q[0,1] + i*b
    
    for i in range(1,p):
        coord33[m*(p-1)+i-1,0] = R*np.cos( (5*math.pi/4) - (i)*((math.pi/2)/p)) + d1/2
        coord33[m*(p-1)+i-1,1] = R*np.sin( (5*math.pi/4) - (i)*((math.pi/2)/p)) + d2/2
    
    for i in range(1, m):
        for j in range(1, p):
            dx = (coord33[m*(p-1)+j-1,0] - coord33[j-1,0])/m
            dy = (coord33[m*(p-1)+j-1,1] - coord33[j-1,1])/m

            coord33[i*(p-1)+j-1,0] = coord33[(i-1)*(p-1)+j-1,0] + dx
            coord33[i*(p-1)+j-1,1] = coord33[(i-1)*(p-1)+j-1,1] + dy
    
    ## region 4
    coord44 = np.zeros([(p-1)*(m+1),PD])

    for i in range(1,p):
        coord44[i-1,0] = q[1,0]
        coord44[i-1,1] = q[1,1] + i*b
    
    for i in range(1,p):
        coord44[m*(p-1)+i-1,0] = R*np.cos( (7*math.pi/4) + (i)*((math.pi/2)/p)) + d1/2
        coord44[m*(p-1)+i-1,1] = R*np.sin( (7*math.pi/4) + (i)*((math.pi/2)/p)) + d2/2
    
    for i in range(1, m):
        for j in range(1, p):
            dx = (coord44[m*(p-1)+j-1,0] - coord44[j-1,0])/m
            dy = (coord44[m*(p-1)+j-1,1] - coord44[j-1,1])/m

            coord44[i*(p-1)+j-1,0] = coord44[(i-1)*(p-1)+j-1,0] + dx
            coord44[i*(p-1)+j-1,1] = coord44[(i-1)*(p-1)+j-1,1] + dy
    
    # reordering the nodes

    for i in range(1,m+2):
        NL[(i-1)*4*p:i*4*p,:] = np.vstack([coord11[(i-1)*(p+1):(i)*(p+1),:],
                                           coord44[(i-1)*(p-1):(i)*(p-1),:],
                                           np.flipud(coord22[(i-1)*(p+1):(i)*(p+1),:]),
                                           np.flipud(coord33[(i-1)*(p-1):(i)*(p-1),:])])

    ## ELEMENTS ##

    EL = np.zeros([NoE,NPE], dtype=int)

    for i in range(1,m+1):
        for j in range(1,4*p+1):
            if j == 1:
                EL[(i-1)*(4*p)+j-1,0] = (i-1)*(4*p) + j
                EL[(i-1)*(4*p)+j-1,1] = EL[(i-1)*(4*p)+j-1,0] + 1
                EL[(i-1)*(4*p)+j-1,3] = EL[(i-1)*(4*p)+j-1,0] + (4*p)
                EL[(i-1)*(4*p)+j-1,2] = EL[(i-1)*(4*p)+j-1,3] + 1
            elif j == 4*p:
                EL[(i-1)*(4*p)+j-1,0] = i*(4*p)
                EL[(i-1)*(4*p)+j-1,1] = (i-1)*(4*p)+1
                EL[(i-1)*(4*p)+j-1,2] = EL[(i-1)*(4*p)+j-1,0] + 1
                EL[(i-1)*(4*p)+j-1,3] = EL[(i-1)*(4*p)+j-1,0] + (4*p)
            else:
                EL[(i-1)*(4*p)+j-1,0] = EL[(i-1)*(4*p)+j-2,1]
                EL[(i-1)*(4*p)+j-1,3] = EL[(i-1)*(4*p)+j-2,2]
                EL[(i-1)*(4*p)+j-1,2] = EL[(i-1)*(4*p)+j-1,3] + 1
                EL[(i-1)*(4*p)+j-1,1] = EL[(i-1)*(4*p)+j-1,0] + 1

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