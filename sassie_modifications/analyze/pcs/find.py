import numpy as np
def mfind (x, target, negate = False) { # if negate is true, it means we are looking for elements of x that are not equal to target
    x = np.transpose(np.nonzero(x))
    numrows = np.shape(x, 0)
    numcols = np.shape(x, 1)
    output = np.array([])
    i = 1
    for col in range(numcols):
        for row in range(numrows):
            value = x[row][col]
            if ((value == target) != negate)):
                np.append(output, i)
            i += 1
    return output
}