import numpy as np

def vector2tensor(X_vector=None): 
    #-----------------------------------------------------
    #   THIS SCRIPT REQUIRES CAREFUL CHECKING OF EIG FUNCTION
    #-----------------------------------------------------
    #--------Reconstruct the tensor-matrix---------------

#    S = float('nan') * np.ones((3,3))
    S = np.zeros( (3,3))
    S[0, 0] = -X_vector[0] - X_vector[1]
    S[1, 1] = X_vector[0]
    S[2, 2] = X_vector[1]
    S[0, 1] = X_vector[2]
    S[0, 2] = X_vector[3]
    S[1, 2] = X_vector[4]
    S[2, 0] = S[0, 2]
    S[2, 1] = S[1, 2]
    S[1, 0] = S[0, 1]

    tensor = S

    '''
    #---------Diagonalize the tensor-matrix--------------

    [d, v] = np.linalg.eig(S) #v and d are arrays, check if eig function returns same as matlab
    v[:,1] = -v[:,1] #temporary fix to direction of 2nd colmn vector

    #------Sort the eigenvalues by absolute value--------
    rotmat = np.real(v) #combine eigenvectors into rotation matrix
    d = np.real(d)
    ii = np.argsort(abs(d)) #np.sort actually doesn't return anything, need to write logically equivalent 
    d = d[ii-1]
    d = d.reshape(len(d),1) #only need this line if we care about whether d is a row or column vector, though np doesnt seem to differentiate the two
    v = v[:, ii-1]
 
    eigenvalues = [d[0], d[1], d[2]]

    print (eigenvalues)
    print (rotmat)
    '''

    return tensor
'''
inputvector = np.array(([1, 2, 3, 4, 5, 6, 7, 8]))
tensor = vector2tensor(inputvector)
print (tensor)
'''
