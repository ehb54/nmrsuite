import numpy as np
import math as m
def generate_xyz_mtrx_rdc(coords=None):
    #------------------------------------------------------------------- 
    #   Calculates the XYZ matrix for RDCs from the input coordinates set
    #---------------------------------------------------------------------
    #eps=0;
    
    X = np.copy(coords[:,0])
    Y = np.copy(coords[:,1])
    Z = np.copy(coords[:,2])
    R2 = (X ** 2 + Y ** 2 + Z ** 2)
    R = np.sqrt(R2)

    X /= R
    Y /= R
    Z /= R

    A = np.ones((np.shape(X)[0], 5))

    A[:, 0] = (Y ** 2 - X ** 2)
    A[:, 1]= (Z ** 2 - X ** 2) 
    A[:, 2] = (2 * X * Y)
    A[:, 3] = (2 * X * Z)
    A[:, 4] = (2 * Y * Z)

    return A * -3

    #=======================================================================
