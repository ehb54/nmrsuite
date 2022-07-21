import numpy as np

from predict_scripts.generate_xyz_mtrx_pcs import generate_xyz_mtrx_pcs
#from generate_xyz_mtrx_pcs import generate_xyz_mtrx_pcs

def tensor2pcs (coords, tensorvector, SLposition):
#----------------------------------------------------
#   df-oct-18 original matlab
#   mc-jun-20 translated to python
#   given the location of the lanthanide, and the 
#   susceptibility tensor in column-vector form
#   calculate the PCS values for each atom
#   defined in the coordinates list
#   INPUT:
#       coords = [res# x y z]  (nat x 4) array
#----------------------------------------------------

    #calculate atom coordinates relative to the SL

    coor = coords[:, 1:4]
    rlist = np.array(coords[:, 0])
    X_vect = tensorvector.flatten("F") # flatten in column major orders

    #return coor, SLposition

    coor[:,0] -= SLposition[0]    #shift in X
    coor[:,1] -= SLposition[1]    #shift in X
    coor[:,2] -= SLposition[2]    #shift in X

    #generate the XYZ matrix
    xyz_mtrx = generate_xyz_mtrx_pcs(coor)
    calc_PCS=xyz_mtrx@X_vect

    z = np.column_stack((rlist,calc_PCS))

    return z

#====================================================
