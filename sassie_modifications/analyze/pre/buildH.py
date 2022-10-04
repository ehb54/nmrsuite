import numpy as np
import math

from sassie.analyze.pre.rotate_df import rotate_df

def buildH(vCaN, vNC, angle  =  121.9947):
    xx = np.zeros(3,)
    yy = np.zeros(3,)
    zz = np.zeros(3,)
    nr1 = vCaN.shape[0]
    nr2 = vNC.shape[0]
    vHN = np.zeros((nr1,3))
    for i in range(nr1):
        Zax = vCaN[i, :]
        REF = vNC[i, :]
        norm = math.sqrt(np.sum(Zax ** 2))
        zz = Zax/norm
        xx[0] = zz[2] * REF[1] - REF[2] * zz[1]
        xx[1] = zz[0] * REF[2] - REF[0] * zz[2]
        xx[2] = zz[1] * REF[0] - REF[1] * zz[0]
        norm = math.sqrt(np.sum(xx ** 2))
        xx = xx/norm
        yy[0] = zz[1] * xx[2] - xx[1] * zz[2]
        yy[1] = zz[2] * xx[0] - xx[2] * zz[0]
        yy[2] = zz[0] * xx[1] - xx[0] * zz[1]
        rrr = rotate_df(np.vstack(xx, yy, zz), 90 - angle, 1)
        vHN[i, :] = rrr[1,:]
    return vHN
