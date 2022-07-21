import numpy as np
import math as m

def tensor_axrh2xyz(Aax, Arh):
    #-------------------------------------
    #   for a traceless tensor, convert
    #   the axial and rhombic values into
    #   the Axx, Ayy, Azz components
    #-------------------------------------
    #Arh=Azz-Ayy Aax=Azz-(Axx+Ayy)/2=Azz*3/2

    Axx = -Aax / 3 + Arh / 2
    Ayy = -(Aax / 3 + Arh / 2)
    #Azz=-(Axx+Ayy)
    Azz = Aax * 2 / 3
    z = np.array([Axx, Ayy, Azz])

    return z
    #=====================================


def euler2rotmat (a, b, g):
    #----------------------------------------------------
    #  Calculates the rotation matrix for a given set
    #  of Euler angles alpha, beta and gamma.

    #  NOTE:This matrix corresponds to rotation of coordinates of a given 
    #  vectors, i.e. R*x, not rotation of the coordinate frames.
    #  Rotation of the coordinate frames is described by 
    #  transpose(R), see Arfken book
    #
    #  The angles are in degrees!   
    #-----------------------------------------------------
    
    z = np.zeros ((3, 3)) # sets up an array of size 3 x 3 and fills it with zeros

    a = m.radians(a)
    b = m.radians(b)
    g = m.radians(g)

    z[0,0] = m.cos(a) * m.cos(b) * m.cos(g) - m.sin(a) * m.sin(g)
    z[0,1] = m.sin(a) * m.cos(b) * m.cos(g) + m.cos(a) * m.sin(g)
    z[0,2] = -m.sin(b) * m.cos(g)
    z[1,0] = -m.cos(a) * m.cos(b) * m.sin(g) - m.sin(a) * m.cos(g)
    z[1,1] = -m.sin(a) * m.cos(b) * m.sin(g) + m.cos(a) * m.cos(g)
    z[1,2] = m.sin(b) * m.sin(g)
    z[2,0] = m.cos(a) * m.sin(b)
    z[2,1] = m.sin(a) * m.sin(b)
    z[2,2] = m.cos(b)

    return z

    #=====================================================

def tensor_diag2full (tensor_diag, alpha, beta, gamma):
    # Converts the tensor containing Axx, Ayy, Azz and the Euler angles into a full 3x3 tensor

    R = euler2rotmat(alpha, beta, gamma)
    tensor = R.conj().T @ np.diag(tensor_diag) @ R # @ is matrix multiplication
    return tensor
    

def tensor2vector (tensor):
    # Converts full tensor into column vector
    return np.array([[tensor[1,1], tensor[2, 2], tensor[0,1], tensor[0,2], tensor[1,2]]])   


def components2tensor (Aax, Arh, a, b, g):
    z = tensor_axrh2xyz (Aax, Arh)
    tensor = tensor_diag2full(z, a, b, g)
    vector = tensor2vector(tensor)
    return vector

Aax = -0.70593
Arh = -0.41931
a = 120.9908
b = 115.357
g = 92.958

print (components2tensor(Aax, Arh, a, b, g))

#-0.20733     0.27613     0.14986    -0.16548     0.31343
