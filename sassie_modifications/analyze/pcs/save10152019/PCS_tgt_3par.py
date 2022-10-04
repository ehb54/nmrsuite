import numpy as np

#import invert_A_new as inv_new
#import generate_xyz_mtrx_pcs as gen
#import x_tensor3 as xt

import sassie.analyze.pcs.invert_A_new as inv_new
import sassie.analyze.pcs.generate_xyz_mtrx_pcs as gen
import sassie.analyze.pcs.x_tensor3 as xt

def PCS_tgt_3par(guess=None, data_PCS=None, coords=None): #coord_XL = vertices from gen, data_PCS = PCS_data.dat, coords = coordinates_H.dat
    #--------------------------------------------------------------------
    #
    #-------------------------------------------------------------------- 
    
    #First 3 values of guess are x,y,z coordinates of spin-label location
    coord_XL = np.copy(np.asarray(guess[0:3]))
    #Create copy of coordinates to prevent constantly overwriting the original file
    coord_mat = np.copy(coords)

    for i in range(3):
        coord_mat[:, i] -= coord_XL[i]

    xyz_mtrx = gen.generate_xyz_mtrx_pcs(coord_mat)
    A_inv = np.linalg.pinv(xyz_mtrx)

    tensor_vector = np.matmul(A_inv, data_PCS)  
    backcalc_PCS = np.matmul(xyz_mtrx, tensor_vector)

    diff = (data_PCS - backcalc_PCS)
    try:
        err = np.copy(data_PCS[:, 2]) #Third column in PCS data should be error values
        Chi2 = np.array([np.sum((diff / err) ** 2)])
    except:
        Chi2 = np.array([np.sum(diff ** 2)])

    #Calculate target function
    return Chi2
    #=====================================================
