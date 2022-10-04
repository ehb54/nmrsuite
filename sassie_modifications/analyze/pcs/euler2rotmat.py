import numpy as np
import math

def euler2rotmat(a, b, g):
    #----------------------------------------------------
    #  df-ovt-2017 (comments added)
    #  Calculates the rotation matrix for a given set
    #  of Euler angles alpha, beta and gamma.
    #
    #  NOTE:This matrix corresponds to rotation of coordinates of a given 
    #  vectors, i.e. R*x, not rotation of the coordinate frames.
    #  Rotation of the coordinate frames is described by 
    #  transpose(R), see Arfken book
    #
    #  The angles are in degrees!   
    #-----------------------------------------------------
    
    z = np.zeros ((3, 3)) # sets up an array of size 3 x 3 and fills it with zeros
    
    def sin (theta): return math.sin(math.radians(theta)) # sine of an angle in degree form
    def cos (theta): return math.cos(math.radians(theta)) # cosine of an angle in degree form
    
    z[0][0], z[0][1], z[0][2] = cos(g) * cos(b) * cos(a) - sin(g) * sin(a), -sin(g) * cos(b) * cos(a) - cos(g) * sin(a), sin(b) * cos(a)
    z[1][0], z[1][1], z[1][2] = cos(g) * cos(b) * sin(a) + sin(g) * cos(a), -sin(g) * cos(b) * sin(a) + cos(g) * cos(a), sin(b) * sin(a)
    z[2][0], z[2][1], z[2][2] = -cos(g) * sin(b), sin(g) * sin(b), cos(b)

    return 

    #=====================================================
   
