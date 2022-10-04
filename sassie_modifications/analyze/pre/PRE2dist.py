import numpy as np
import math
import scipy

def PRE2dist(ratio, T2dia, Htime, freq, TAUc):


    omega = freq * 2 * math.pi * 1e6
    d2 = 1.23e-44
    d2 = d2 * (1e10 ** 6)
    tauC = TAUc * 1e-9
    R2dia = 1/float(T2dia)
    factor = d2 * (4 * tauC+3 * tauC/float(1+(omega * tauC)^2))

    nres = ratio.shape[0]
    z = np.full((nres, 3), float('nan'))
    z[:, 0] = ratio[:, 0]
    for i in range(nres):
        if ~np.isnan(ratio[i, 1]):
            if ratio[i, 1] == 0:
                ratio[i, 1] = 0.00000001   
            if ratio[i, 1] > 1:
                ratio[i, 1] = 0.99999999 
            X0 = R2dia * (1 - ratio[i, 1]) / float(ratio[i, 1] + Htime * R2dia)
            X = scipy.optimize.fsolve(lambda x: ratio[i, 1] * (x + R2dia) - R2dia * np.exp(-x * Htime), 2)
            z[i, 1] = X

    z[:, 2] = np.divide(factor, z[:, 1] ** (1/6.0))

    return z
