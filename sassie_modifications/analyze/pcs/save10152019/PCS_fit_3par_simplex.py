import numpy as np

#import pcs_simplex_3par as pcss3
#import generate_vertices as genfile
#import generate_xyz_mtrx_pcs as genmat

import sassie.analyze.pcs.pcs_simplex_3par as pcss3
import sassie.analyze.pcs.generate_vertices as genfile
import sassie.analyze.pcs.generate_xyz_mtrx_pcs as genmat

def PCS_fit_3par_simplex(self, PCS=None, coord=None, guess=None, tolerance = None):
    #-------------------------------------------------------------
    #   df-aug-18
    #   given PCS data and atom coordinates 
    #   determine the position of the lantanide (through simplex)
    #   and the DeltaX tensor through SVD  
    #   essentially a 3-parameter fit (lanthanide coordinates)
    #-------------------------------------------------------------

    #need a few lines that make sure the residue numbers in PCS and Coordinates match

    
    [position, Chi2] = pcss3.pcs_simplex_3par(self, guess, PCS, coord, tolerance)#simplex minimization

    #position
    #Chi2
    coord_s = np.copy(coord)
    for i in range(3):
        coord_s[:, i] -= position[i]#shift in X,Y,Z

    xyz_mtrx = genmat.generate_xyz_mtrx_pcs(coord_s)
    cond_number = np.linalg.cond(xyz_mtrx)

    return (position, Chi2, cond_number)
