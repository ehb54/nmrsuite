import numpy as np
import math as m

from predict_scripts.generate_xyz_mtrx_rdc import generate_xyz_mtrx_rdc
#from generate_xyz_mtrx_rdc import generate_xyz_mtrx_rdc

def tensor2rdc(coords, tensorvector, Sorder, freq = 850, T = 298):
    #----------------------------------------------------
    #   Python: mc-jun-2020
    #   Matlab: df-jun-2020  df-oct-18
    #   given the  susceptibility tensor in column-vector form
    #   calculate the RDC values for each atom-pair
    #   defined in the coordinates list
    #   INPUT:
    #       coords = [res# x y z]  (nat x 4) array
    #       tensorvector = [Xyy,Xzz,Xxy,Xxz,Xyz] column vector (in 10^-32 units)
    #       Sorder
    #       freq = spectrometer 1H Larmor frequency (in MHz units)
    #       T = temperature (in K) (optional)
    #----------------------------------------------------

    gH=267.522*1e6
    gN=-27.126*1e6         #It looks like we need the sign here
    Bo=freq/gH*2*m.pi*1e6

    #S=1    
    h_=1.0536e-034         #Planck's const/2pi
    kB= 1.3800e-023        #Boltzmann const
    r0=1.04e-10            #length of NH vector, to keep the order of magnitude right
    tensor_units=1e-32     #X-tensor units, to keep the order of magnitude right
    mu0_4pi=1e-7           #not used, for now

    Dmax = -1*Bo**2*gH*gN*h_*Sorder/(120*kB*T)/m.pi/m.pi/r0**3*tensor_units   
    #Note: the factor of 3 that comes from conversion of the tensor matrix equation into 
    #a tensor-vector equation is already included in the generated xyz matrix

    #need to make sure the output of generate_matrx_pcs doesn't have extra scalings
    #return

    #N-H-vector coordinates
    coor = coords[:, 1:4]
    rlist = coords[:, 0]
    X_vect = tensorvector.flatten("F") # flatten in column major order

    #generate the XYZ matrix
    xyz_mtrx = generate_xyz_mtrx_rdc(coor)

    #calculate the RDCs
    calc_RDC = xyz_mtrx@X_vect*Dmax

    #z = np.column_stack((coords, tensorvector))

    z = np.column_stack((rlist,calc_RDC))

    #z = Dmax

    return z

    #====================================================

    #S=0.9, B0=850 MHz, gH=267.522 106 rad.s-1.T-1, gN=-27.126 106 rad.s-1.T-1 and T=298 K. 
