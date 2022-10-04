#from helperFunctions import getResidueNumbers
#from Predict.Predict.readpdb import readpdb
from predict_scripts.readpdb import readpdb
from predict_scripts.pdb2nhcoor import pdb2nhcoor

import numpy as np
import math as m

""" def predictPRE (atom_coordinates, atom_types, XYZ_coordinates, frequency, T2dia, OneHT2, TauC):
    residue_nums = getResidueNumbers(3) # Change parameter to modify input type

    for atom_type in atom_types:
        print("Results for Atom Type " + str(atom_type))
        predictPREHelper(atom_coordinates, residue_nums, atom_type, XYZ_coordinates, frequency, T2dia, OneHT2, TauC)
        print("------------------------")
 """
#def predictPREHelper (atom_coordinates, residue_nums, atom_type, XYZ_coordinates, freq, T2dia, OneHT2, TauC):


def predictPRE (reslist, T2dia, TAUc, Htime, SL_position, freq, pdb_filename, spin = 1/2, gammaRatio = 1, pdb_model = 0, chainID = "A", atom_type = "H"):
    # DO PCS stuff
    omega = freq * gammaRatio * 2 * m.pi * 1e6   
    S = spin
    d2 = (gammaRatio  ** 2) * S * (S+1) * 1.6414e-44
    d2 = d2 * (1e10 ** 6)
    tauC = TAUc * 1e-9
    R2dia = 1/float(T2dia)
    factor = d2 * (4 * tauC + 3 * tauC / float(1 + (omega * tauC) ** 2))

    readpdb_output = readpdb(fname = pdb_filename, reslst = reslist, atlst = np.array([atom_type]), model_num = 0, chainID = chainID)
    NH_coord = readpdb_output[0]
    reslst = readpdb_output[2]
    #pdb2nhcoor(fname = pdb_filename, model = pdb_model, chainID = chainID)
    #return (NH_coord, NH_coord)
    #a,b = pdb2nhcoor(fname = pdb_filename, model = pdb_model, chainID = chainID)
    #return np.column_stack((a, b))

    # translation for pymol function slfit_ss -> copy/paste into here
    """ Hvect = NH_coord[i, 4:7] - SL_position
    dist = np.sqrt(np.dot(Hvect, np.conj(Hvect.reshape(Hvect.size, 1))))
    R2para = factor/(dist ** 6)
    ratio = R2dia * np.exp(-R2para * Htime) / float(R2para + R2dia) """

    """ for i in range(nres):
        Hvect = NH_coord[i, 4:7] - position
        dist = np.sqrt(np.dot(Hvect, np.conj(Hvect.reshape(Hvect.size, 1))))
        R2para = factor/(dist ** 6)
        ratio_sim = R2dia * math.exp(-R2para * Htime) / float(R2para + R2dia)
 """
    
    at_coord = NH_coord#[:, [0, 4, 5, 6]]

    a = at_coord[:, 0] - SL_position[0]
    b = at_coord[:, 1] - SL_position[1]
    c = at_coord[:, 2] - SL_position[2] 

    dist = np.sqrt(a ** 2 + b ** 2 + c ** 2)

    R2para = np.divide(factor, dist ** 6)

    ratio_sim = np.divide(R2dia * np.exp(-R2para * Htime), (R2para + R2dia))

    #np.savetxt("Predict/ratio_sim.txt", ratio_sim)
    #np.savetxt("ratio_sim.txt", ratio_sim) 
    
    return np.column_stack((reslst, ratio_sim))

""" T2dia = 50e-3
TAUc = 4.5
Htime = 5e-3
SL_position = np.array([48.88335461474581, -67.33450318575808, -8.346363336662463])
freq = 600.13
pdb_filename = "./1d3z_mod1.pdb"
spin = .5

print (predictPRE (T2dia, TAUc, Htime, SL_position, freq, pdb_filename, spin)) """

# ratio = 1
""" 
[' X = ', '48.88335461474581', ' A']
[' Y = ', '-67.33450318575808', ' A']
[' Z = ', '-8.346363336662463', ' A']


['Corr Coeff = ', '0.98733']
['qR_factor = ', '0.11853']
['electron spin S = ', '0.5'] """

"""
0.97251046 0.93718155 0.8526038  0.80143928 0.4942296  0.57206116
 0.23573874 0.47472837 0.56230937 0.70690369 0.87596817 0.813073
 0.93169002 0.92533446 0.95913044 0.96257553 0.96049093 0.97079898
 0.95420255 0.93195836 0.80866631 0.79528525 0.87421968 0.84524815
 0.74901942 0.84113179 0.88810193 0.83839116 0.8412932  0.91665811
 0.92252845 0.88812247 0.89539337 0.81320067 0.62822768 0.42654379
 0.24108574 0.02245113 0.03695149 0.03860665 0.08617862 0.31228882
 0.09210954 0.0855569  0.00659599 0.05784676 0.37574196 0.46758333
 0.6933994  0.74689526 0.87528875 0.89260829 0.9261074  0.90420119
 0.84405851 0.8901388  0.87842582 0.90671325 0.96374844 0.95197559
 0.91569049 0.78399516 0.63852387 0.17865229 0.15490812 0.01119569
 0.00354158 0.03738071 0.02871575 0.34178434 0.49127205 0.65577518
 """
 
 # Make sure everything lines up in the plots, match, on y=x, Chi squared close to 0
