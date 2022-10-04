#from sassie import matlabeig

import numpy as np
import matlabeig as me
def tensor2eigen_sort (S, ksort = 0):
    #--------------------------------------------------------
    # df-jul-19, df-jun-19 
    #   for input tensor/matrix S, 
    #   determine its eigenvalues and eigenvectors 
    #   sort the eigensystem based on ksort options
    #   check and restore (if necessary) right-handedness of eigenvectors
    #   and determnine the rotation matrix
    #   
    #  INPUT: ksort=1 sort by absolute value
    #         ksort= -1 sort by value
    #         ksort not equal 1 or -1  do nothing (default)
    #  OUTPUT eigenval -- row of eigenvalues
    #         eigenvect -- 3x3 matrix where columns are eigenvectors
    #         rotmat -- 3x3 rotation matrix rotmat=transpose(eigenvect)
    #            rotmat defined as in Arfken, i.e. passive rotation
    #(c) Copyright by David Fushman, U.Maryland
    #--------------------------------------------------------

    #---------Diagonalize the S matrix-----------------------
    [v,d] = me.matlabeig(S)
    dia = np.diag(np.real(d))

    #---------Sort the eigensystem------------------
    if ksort == 1:           #sort by absolute value
        dia_s, ind = np.sort(abs(dia)), np.argsort(abs(dia))   
    elif ksort == -1:   
        dia_s, ind = np.sort(dia), np.argsort(dia) #sort by value not by abs value
    else:
        dia_s=dia
        ind = [0, 1, 2]
    eigenval = dia[ind] #make it a row

    #rearrange/permute the eigenvectors, if necessary
    eigenvect=v[:,ind]           #
    
    #check and restore (if necessary) right-handedness of eigenvectors
    rhand_check=np.transpose(np.dot(np.cross(eigenvect[:, 0], eigenvect[:, 1]), eigenvect[:, 2]))
    if rhand_check < 0:
        eigenvect[:, 1] = -eigenvect[:, 1]
        print('Eigenvectors: changed sign of 2nd eigenvector to keep right-handedness')
    #ideally, need to change sign when odd number of permutations, but I found
    #that sometimes eig returns left-handed eigenvectors, so it's safer just
    #to check for and restore right-handedness regardless of permutations

    #determine rotmat as in Arfken, i.e. passive rotation
    rotmat = np.transpose(eigenvect)       

    return eigenval, eigenvect, rotmat

'''
test = np.array([[0.8166, -1.3445, -3.3450], [-1.3445, 0.9318, -3.1860], [-3.3450, -3.1860, -1.7483]])
eigenval, eigenvect, rotmat = tensor2eigen_sort(test)
print (eigenval)
print (eigenvect)
print (rotmat)
'''