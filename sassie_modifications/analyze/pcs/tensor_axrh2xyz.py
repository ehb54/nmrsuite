import numpy as np
def tensor_axrh2xyz(Aax=None, Arh=None):
    #-------------------------------------
    #   for a traceless tensor, convert
    #   the axial and rhombic values into
    #   the Axx, Ayy, Azz components
    #-------------------------------------
    #Arh=Azz-Ayy; Aax=Azz-(Axx+Ayy)/2=Azz*3/2;

    Axx = -Aax / 3.0 + Arh / 2.0
    Ayy = -(Aax / 3.0 + Arh / 2.0)
    #Azz=-(Axx+Ayy);
    Azz = Aax * 2.0 / 3.0
    z = np.array([Axx, Ayy, Azz])

    return z
    #=====================================
