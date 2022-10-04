import numpy as np
import math
import matplotlib.pyplot as plt
import time

def SLfit_ss(X, inten, at_coord, factor, R2dia, Htime, kvis = 0):

    a = at_coord[:, 1] - X[0]
    b = at_coord[:, 2] - X[1]
    c = at_coord[:, 3] - X[2]

    dist = np.sqrt(a ** 2 + b ** 2 + c ** 2)

    R2para = np.divide(factor, dist ** 6)

    ratio_sim = np.divide(R2dia * np.exp(-R2para * Htime), (R2para + R2dia))

    Chi2 = np.sum((ratio_sim - inten[:, 1]) ** 2)

    if kvis == 1:
        plt.plot(inten[:, 0], inten[:, 1])
        plt.plot(inten[:, 0], ratio_sim)
        plt.pause(0.0001)
        plt.draw()
        plt.clf()

    return Chi2