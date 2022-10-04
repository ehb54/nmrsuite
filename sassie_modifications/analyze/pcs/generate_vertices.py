import numpy as np
import sys

def generate_vertices(coord=None, k_face=None):
    #--------------------------------------------------------------
    #   Given protein coordinates generate vertices of a cube around the protein
    #   INPUT: 
    #       1) coord = [x,y,z]  nx3 array, n = number of residues
    #       2) k_face = 0 (default) only generate 8 corners
    #                   1 add face-points (shifted 0.7 away from the face)
    #--------------------------------------------------------------
    if k_face is None:
        k_face = 0
    #default: only corners

    if (coord.shape[1] >= 4):
        coord = np.delete(coord, 0, axis=1) #If coordinate input contains more than 3 columns, assume the first column is residue number
    
    coord = np.array(coord) # converts standard array to numpy array
    
    mx = np.max(coord, axis=0)
    mn = np.min(coord, axis=0)
    mc = np.mean(coord, axis=0)
    
    vertices = np.array([mn, [mn[0], mn[1], mx[2]], [mn[0], mx[1], mn[2]], [mx[0], mn[1], mn[2]], [mn[0], mx[1], mx[2]], [mx[0], mn[1], mx[2]], [mx[0], mx[1], mn[2]], mx])

    if k_face == 1:
        faces = np.array([[mn[0] - (mc[0] - mn[0]) * 0.7, mc[1], mc[2]], [mx[0] + (mx[0] - mc[0]) * 0.7, mc[1], mc[2]], [mc[0], mn[1] - (mc[1] - mn[1]) * 0.7, mc[2]], [mc[0], mx[1] + (mx[1] - mc[1]) * 0.7, mc[2]], [mc[0], mc[1], mn[2] - (mc[2] - mn[2]) * 0.7], [mc[0], mc[1], mx[2] + (mx[2] - mc[2]) * 0.7]]) 
        vertices = np.concatenate([vertices, faces])


    return vertices

#=====================================================

