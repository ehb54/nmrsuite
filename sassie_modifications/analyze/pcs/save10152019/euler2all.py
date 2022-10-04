import numpy as np
def euler2all (alpha, beta = None, gamma = None):
    #-----------------------------------------------------
    #   For a given set of euler angles (in degrees), 
    #   produce all other possible angles
    #   based on the symmetry properties of
    #   the diffusion tensor
    #   {a,b,g}={a,b,g+180o}={a+180o,180o-b,180o-g}={a+180o,180o-b,360o-g}
    #   and select the set that has all three angles in the [0 180]-cube
    #------------------------------------------------------    
    if (beta == None or gamma == None):
        if (beta == None and gamma == None):
            if (np.size(np.asarray(alpha)) == 3):
                beta = np.copy(alpha[1])
                gamma = np.copy(alpha[2])
                alpha = alpha[0]
            else:
                print ("wrong input!")
        else:
            print ("wrong input!")
    z = np.full((4, 3), float('nan'))
    zsel = float('nan')
    z[0, :] = [alpha,beta,gamma]
    z[1, :] = [alpha,beta,gamma+180]
    z[2, :] = [alpha+180,180-beta,180-gamma]
    z[3, :] = [alpha+180,180-beta,360-gamma]

    #convert angles to "normal" range
    z = np.where(z >= 270, z - 360, z) # changes to a - 360 if a >= 270 and keeps as a else
    z = np.where(z <= -270, z + 360, z)
    
    #select
    indsel = np.array([])
    for i in range(4):
        alpha = z[i, 0]
        beta = z[i, 1]
        gamma = z[i, 2]
        indomain = lambda angle : angle >= 0 and angle <= 180
        if (indomain(alpha) and indomain(beta) and indomain(gamma)):
            indsel = np.append(indsel, i)
    if (indsel.size > 0):
        zsel = z[indsel.astype(int), :]
    return z, zsel[0]
    #======================================================
