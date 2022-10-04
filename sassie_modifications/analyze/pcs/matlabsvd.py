import numpy as np
def matlabsvd (A, full_matrices=True, compute_uv=True):
    [u, w, v] = np.linalg.svd(A, full_matrices=True, compute_uv=True)
    w = np.diag(w ** -1)
    v = np.transpose(v)
    return (u, w, v)
'''
a = np.array([[1, 2], [3, 4]])
[u, w, v] = matlabsvd(a)
print ('u=',u)
print ('w=',w)
print ('v=',v)
print (np.matmul(np.matmul(u, w), v))
inv = np.matmul(np.matmul(u,w), v)
print ('inv=',inv)
def arrayOfArrayToMatlabString(array):
    return '[' + "\n ".join(" ".join("%6g" % val for val in line) for line in array) 
'''