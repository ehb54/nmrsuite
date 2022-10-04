import numpy as np
import math

def rotate_df(orts, angle, axis):
    cs = math.cos(math.radians(angle))
    si = math.sin(math.radians(angle))
    r = orts
    if axis == 1: 		
        rot = np.array([[1, 0, 0], [0, cs, si], [0, si, -cs]])
    elif axis == 2:
        rot = np.array([[cs, 0, -si], [0, 1, 0], [si, 0, cs]])
    else:
        rot = np.array([[cs, si, 0], [-si, cs, 0], [0, 0, 1]])    
    orts1 = orts * r.H * rot * r
    return orts1
