import numpy as np

'''
This program is used to subtract column vectors from matrices
Performs the subtraction as matlab would
'''
def col_sub(mat = None, vec = None):
    if mat.shape[1] == vec.shape[0]:
        for i in range(mat.shape[1]):
            mat[:, i] -= vec[:]
    else: 
        print('Error in dimensions of input')
        return
    return mat
