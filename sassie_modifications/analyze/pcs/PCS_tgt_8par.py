import numpy as np

#import generate_xyz_mtrx_pcs as gxmpfile
import sassie.analyze.pcs.generate_xyz_mtrx_pcs as gxmpfile

def PCS_tgt_8par(guess=None, data_PCS=None, coords=None):
    #--------------------------------------------------------------------
    #
    #--------------------------------------------------------------------

    #subtract the XL coordinates to get the atom--XL vector

    #First 3 values of guess are x,y,z coordinates of spin-label location
    coord_XL = np.copy(np.asarray(guess[0:3]))
    #Create copy of coordinates to prevent constantly overwriting the original file
    coord_mat = np.copy(coords)

    for i in range(3):
        coord_mat[:, i] -= coord_XL[i]

    tensor_vect = np.copy(guess[3:8])

    xyz_mtrx = gxmpfile.generate_xyz_mtrx_pcs(coord_mat)
    
    if np.shape(xyz_mtrx)[1] == np.shape(tensor_vect)[0]: calc_PCS = np.matmul(xyz_mtrx, np.transpose(tensor_vect))      
    else: calc_PCS = np.matmul(xyz_mtrx, tensor_vect)

    diff = np.subtract(data_PCS, calc_PCS)
    try:
        err = np.copy(data_PCS[:, 2]) #Third column in PCS data should be error values
        Chi2 = np.array([np.sum((diff / err) ** 2)])
    except:
        Chi2 = np.array([np.sum(diff ** 2)])

    return Chi2
    #=====================================================
