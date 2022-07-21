import numpy as np
import math as m

def generate_xyz_mtrx_pcs(coords=None):
    #------------------------------------------------------------------- 
    #   Calculates the A matrix from the input co-ordinate set
    #   In this current ver. coords must be an array, not a list
    #   Note that element operation is default
    #---------------------------------------------------------------------
    #eps=0;
    
    X = np.copy(coords[:,0])
    Y = np.copy(coords[:,1])
    Z = np.copy(coords[:,2])
    R2 = (X ** 2 + Y ** 2 + Z ** 2)
    R = np.sqrt(R2)
    R5 = R * R2 ** 2 + np.finfo(float).eps #added eps to avoid zero division, eps requires numpy

    A = np.full((X.size, 5), float('nan'))
#DF: just define a matrix of this size/

    A[:, 0] = (Y ** 2 - X ** 2) / R5
    A[:, 1]= (Z ** 2 - X ** 2) / R5
    A[:, 2] = (2 * X * Y) / R5
    A[:, 3] = (2 * X * Z) / R5
    A[:, 4] = (2 * Y * Z) / R5
    A = A / 4 / m.pi * 10000

    return A

    #=======================================================================
