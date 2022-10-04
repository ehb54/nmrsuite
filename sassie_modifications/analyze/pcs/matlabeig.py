import numpy as np
def matlabeig (matrix):
    w, v = np.linalg.eig(matrix)
    d = np.diag(w)
    return v, d
def testmatlabeig (matrix):
    v, d = matlabeig(matrix)
    print (v)
    print (d)
'''
a = np.array([[0.8166, -1.3445, -3.3450], [-1.3445, 0.9318, -3.1860], [-3.3450, -3.1860, -1.7483]])
testmatlabeig(a)
b = np.array([[0.74899068, 3.50879255, 3.92726324], [0.96852007, -2.25680279, -2.19297766], [-2.36688388, 3.49953363, 2.46581692]])
testmatlabeig(b)
'''