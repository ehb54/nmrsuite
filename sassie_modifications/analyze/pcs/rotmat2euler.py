import math as m
import numpy as np

def rotmat2eulerZYZ(rotmat): # backcalculates the euler angles (in ZYZ convention) of a rotation matrix
    rotmat = np.asarray(rotmat)
    beta = m.acos(rotmat[2, 2])
    alpha = m.atan2(rotmat[2, 1], rotmat[2, 0])
    gamma = m.atan2(rotmat[1, 2], -rotmat[0, 2])
    return [m.degrees(alpha), m.degrees(beta), m.degrees(gamma)]

def rotmat2eulerZYX(rotmat): # backcalculates the euler angles (in ZYX convention) of a rotation matrix
    rotmat = np.asarray(rotmat)
    beta = m.asin(-rotmat[0, 2])
    alpha = m.atan2(rotmat[0, 1], rotmat[0, 0])
    gamma = m.atan2(rotmat[1, 2], rotmat[2, 2])
    return [m.degrees(alpha), m.degrees(beta), m.degrees(gamma)]

def rotmat2eulerXYZ(rotmat): # backcalculates the euler angles (in XYZ convention) of a rotation matrix
    rotmat = np.asarray(rotmat)
    beta = m.acos(rotmat[2, 0])
    alpha = m.atan2(-rotmat[2, 1], rotmat[2, 2])
    gamma = m.atan2(-rotmat[1, 0], rotmat[0, 0])
    return [m.degrees(alpha), m.degrees(beta), m.degrees(gamma)]

round_digits = 4

def unpack (angles): # turns an array containing the three euler angles into a user-friendly output
#    return "Alpha =", round(angles[0], round_digits), "Beta =", round(angles[1], round_digits), "Gamma =", round(angles[2], round_digits)
    out = ""
    out += "          Alpha = " + str(round(angles[0], round_digits)) + "\n"
    out += "          Beta  = " + str(round(angles[1], round_digits)) + "\n"
    out += "          Gamma = " + str(round(angles[2], round_digits)) + "\n"

    return out 

def printeuler (rotmat, stdprint = True): # computes, unpacks, and identifies the euler angles for each of the three rotation conventions
    output = "\n"
    output += "    1) ZYZ Convention: \n" 
#    output += "".join(str(unpack(rotmat2eulerZYZ(rotmat))))
    output += unpack(rotmat2eulerZYZ(rotmat))
    output += "    2) ZYX Convention: \n"
    output += unpack(rotmat2eulerZYX(rotmat))
#    output += "".join(str(unpack(rotmat2eulerZYX(rotmat))))
    output += "    3) XYZ Convention: \n"
    output += unpack(rotmat2eulerXYZ(rotmat))
#    output += "".join(str(unpack(rotmat2eulerXYZ(rotmat))))
    return output
