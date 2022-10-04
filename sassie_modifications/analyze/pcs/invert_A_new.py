import numpy as np
def invert(A):
    u, w, v = np.linalg.svd(A)
    inverse = np.dot (v.transpose(), np.dot(np.diag(w ** -1), u.transpose()))
    return (inverse)

'''
# Check function
n = 24
A = np.random.randn(n, n)
Ainv = invert(A)
print (np.allclose(Ainv @ A, np.eye((n))))
'''