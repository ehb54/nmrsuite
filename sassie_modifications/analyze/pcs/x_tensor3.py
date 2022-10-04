import numpy as np

def x_tensor3(a_inv=None, d_eff=None):

    q_v = a_inv.dot(d_eff)
    S = np.empty((3,3)) #should give an empty matrix, testing in shell gives matrix with values, hopefully not an issue

    #--------Reconstruct the Q-matrix---------------

    S[0, 0] = -q_v[0] - q_v[1]
    S[1, 1] = q_v[0]
    S[2, 2] = q_v[1]
    S[0, 1] = q_v[2]
    S[0, 2] = q_v[3]
    S[1, 2] = q_v[4]
    S[2, 0] = S[0, 2]
    S[2, 1] = S[1, 2]
    S[1, 0] = S[0, 1]

    #---------Diagonalize the Q matrix--------------

    [v, d] = np.linalg.eig(S)

    #---------Sort the eigensystem------------------

    v = np.real(v)

    d = np.diag(np.real(d),k=0)

    i = np.argsort(abs(d))

    d = d[i]


    rot = v[:, i]


    S_val = [d[0], d[1], d[2]]


    return [rot, S_val,S]

    #====================================================
