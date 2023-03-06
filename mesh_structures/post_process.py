"""
@author: Mad Dreamer
"""
import numpy as np
import math
from functions import gauss_point, grad_N_nat, constitutive

"""
Post process matrix, for mesh structures
"""
def post_process(NL,EL,ENL):
    PD = np.size(NL,1)
    NoE = np.size(EL,0)
    NPE = np.size(EL,1)

    scale = 1 # exagerate deflection

    disp, stress, strain = element_post_process(NL,EL,ENL)

    stress_xx = np.zeros([NPE,NoE])
    stress_xy = np.zeros([NPE,NoE])
    stress_yx = np.zeros([NPE,NoE])
    stress_yy = np.zeros([NPE,NoE])

    strain_xx = np.zeros([NPE,NoE])
    strain_xy = np.zeros([NPE,NoE])
    strain_yx = np.zeros([NPE,NoE])
    strain_yy = np.zeros([NPE,NoE])

    disp_x = np.zeros([NPE,NoE])
    disp_y = np.zeros([NPE,NoE])

    X = np.zeros([NPE,NoE])
    Y = np.zeros([NPE,NoE])

    if NPE in [3,4]: #D2QU4N and D2TR3N
        # calculate X and Y coordinates of the nodes
        X = ENL[EL-1,0] + scale*ENL[EL-1,4*PD]
        Y = ENL[EL-1,1] + scale*ENL[EL-1,4*PD+1]

        X = X.T
        Y = Y.T

        #[[xx,xy],
        # [yx,yy]]

        #[:,:] is to compensate for D2TR3N since it only has 1 gauss point. All the r

        stress_xx[:,:] = stress[:,:,0,0].T
        stress_xy[:,:] = stress[:,:,0,1].T
        stress_yx[:,:] = stress[:,:,1,0].T
        stress_yy[:,:] = stress[:,:,1,1].T

        strain_xx[:,:] = strain[:,:,0,0].T
        strain_xy[:,:] = strain[:,:,0,1].T
        strain_yx[:,:] = strain[:,:,1,0].T
        strain_yy[:,:] = strain[:,:,1,1].T

        disp_x = disp[:,:,0,0].T
        disp_y = disp[:,:,1,0].T

    # you may need to write other if blocks for D2QU8N, D2QU9N and D2TR6N due to their
    
    return (stress_xx, stress_xy, stress_yx, stress_yy, strain_xx, strain_xy, strain_yx, strain_yy, disp_x, disp_y, X, Y)

"""
Post process element, for mesh structures
"""
def element_post_process(NL,EL,ENL):
    PD = np.size(NL,1)
    NoE = np.size(EL,0)
    NPE = np.size(EL,1)

    if NPE == 3:
        GPE = 1
    elif NPE == 4:
        GPE = 4
    else:
        NotImplemented()

    # NoE: specified the element I'm looking at
    # NPE: specified the node on the element I'm looking at
    # PD: Displacement is a vector that has x and y coordinates
    # Specified the direction I'm looking at
    #1: Added to make it similar to stress and strain matrices

    disp = np.zeros([NoE, NPE, PD, 1]) # Disp is written on the nodes

    # Stress and strain are normally supposed to mapped on the Gauss Points
    # but because of the shortcomings of Python and Matlab, we will cheat by
    # mapping these to the corners as well

    # In short, stress and stain are calculated on the Gauss Points but mapped
    # on to the corners.

    stress = np.zeros([NoE, GPE, PD, PD])
    strain = np.zeros([NoE, GPE, PD, PD])

    for e in range(1, NoE+1): # First, find the displacements. By using disps,
        nl = EL[e-1,0:NPE]      #plot streaa and strain

        # Assignning displacements to the corresponding nodes
        for i in range(1,NPE+1):
            for j in range(1,PD+1):
                disp[e-1,i-1,j-1,0] = ENL[nl[i-1]-1,4*PD+j-1]

        # Specify the corners of the elements (just like in the stiffness calculation)
        x = np.zeros([NPE,PD])
        x[0:NPE, 0:PD] = NL[nl[0:NPE]-1,0:PD]

        # specify the displacements for these corners (required for strain)
        u = np.zeros([PD,NPE])
        for i in range(1,NPE+1):
            for j in range(1,PD+1):
                u[j-1,i-1] = ENL[nl[i-1]-1,4*PD+j-1]


        # coordinates of the corners transposed (just like in the stiffness calculation)
        coord = x.T

        # going over gauss points since every gauss point will have their own
        for gp in range(1, GPE+1):

            # strain for each gauss point (2x2 matrix)
            epsilon = np.zeros([PD,PD])

            # going over the nodes in the element
            for i in range(1, NPE+1):
                # The process for the shape functions is the same as in stiffness
                J = np.zeros([PD,PD])

                grad = np.zeros([PD,NPE])

                (xi,eta,alpha) = gauss_point(NPE,GPE,gp)

                grad_nat = grad_N_nat(NPE,xi,eta)

                J = coord @ grad_nat.T

                grad = np.linalg.inv(J).T @ grad_nat

                # calculate strain
                # define dyadic in another function
                epsilon = epsilon + 1/2 * (dyad(grad[:,i-1],u[:,i-1]
                                                ) + dyad(u[:,i-1],grad[:,i-1]))
        
        # initialize stress as a 2x2 matrix
        sigma = np.zeros([PD,PD])

        # The same logic as in the stiffness calculation (4 directions)
        # sigma = E * epsilon
        for a in range(1, PD+1):
            for b in range(1, PD+1):
                for c in range(1, PD+1):
                    for d in range(1, PD+1):
                        sigma[a-1, b-1] = sigma[a-1, b-1] + constitutive(
                                a,b,c,d) * epsilon[c-1,d-1]
        
        # Compile the results. Remember the sizes of stress and strain matrices:
        for a in range(1,PD+1):
            for b in range(1,PD+1):
                strain[e-1,gp-1,a-1,b-1] = epsilon[a-1,b-1]
                stress[e-1,gp-1,a-1,b-1] = sigma[a-1,b-1]

    return disp, stress, strain

"""
Dyad, for mesh structures
"""
def dyad(u,v):
    # Takes two matrixes
    # shapes them into row matrixes
    u = u.reshape(len(v),1)
    v = v.reshape(len(v),1)
    PD = 2
    # Matrix multiplication
    A = u @ v.T
    return A
