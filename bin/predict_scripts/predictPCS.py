# Inputs: 
# Atom Coordinates
# List of Residues
# Type of Atoms
# XYZ Coordinates of the Paramagnetic Center
# Magnetic Susceptibility Tensor

# Output: PCS values for each atom/residue (numerical and also plotted as a function of residue number)

import numpy as np
import math as m
from predict_scripts.components2tensor import components2tensor
from predict_scripts.readpdb import readpdb
from predict_scripts.tensor2pcs import tensor2pcs


""" def predictPCS (atom_coordinates, atom_types, XYZ_coordinates):
    residue_nums = getResidueNumbers(3) # Change parameter to modify input type
    susceptibility_tensor = getSusceptibilityTensor(1) # Change parameter to modify input type

    for atom_type in atom_types:
        print("Results for Atom Type " + str(atom_type))
        predictPCSHelper(atom_coordinates, residue_nums, atom_type, XYZ_coordinates, susceptibility_tensor)
        print("------------------------") """

def predictPCS (susceptibility_tensor, reslist, pdb_filename, SLposition, atom_type = "H", pdb_model = 0, chainID = "A"):

    #reslist = np.array(range(1, 76))

    #convert output of PCS analysis into SL coordinates and the tensor

    #---------------------------------
    #---------------------------------
    result_8par = np.array([]) # MODIFY THIS LINE to specify the output you want to use
    #---------------------------------
    #---------------------------------

    #result_Xpar = np.copy(result_8par)
    #result_Xpar = result_3par

    #tensorvector = np.transpose(result_Xpar[3:8]) # leaving this to test
    tensorvector = susceptibility_tensor

    #read in PDB and extract proton coords for selected residues
    #[coorH,atnamH,at_resH] = readpdb(pdb_filename,reslist,'H ',pdb_model)
    #coordH = [at_resH(:,2),coorH]
    #read coordinates of all amide H and N atoms 
    #[coorH,atnamHall,at_resH] = readpdb(pdb_filename,np.array([2, 3]),np.array(["H"]),model_num=0,chainID="A")
    [coorH,atnamHall,at_resH] = readpdb(pdb_filename,reslist,np.array([atom_type]),pdb_model,chainID)
    #print (readpdb("mnodesS00001.pdb", np.array([2, 4]), np.array(["H", "N"]), 0, "B"))

    
    #return reslist, at_resH

    """ print (coorHall)
    print (atnamHall)
    print (at_resHall) """
    #read in only those amide Ns that have Hs attached

    #[coorNall,atnamNall,at_resNall] = readpdb(pdb_filename,at_resHall[:, 1],'N ',pdb_model)

    #now predict PCSs using this information
    #PCS_pred = tensor2pcs(coordH,tensorvector,SLposition)

    #susceptibility tensor = result_8par[3:8]

    """ print ("Tensor vector: {}".format(tensorvector))
    print ("SL Position: {}".format(SLposition))
    print ("Shapes")
    print (at_resHall[:, 1].shape)
    print (coorHall.shape) """

    coordH = np.column_stack((at_resH, coorH))

    PCS_H_predicted = tensor2pcs(coordH, tensorvector, SLposition)

    return PCS_H_predicted

    '''
    PCS_N_predicted = tensor2pcs([at_r  esNall[:, 1], coorNall], tensorvector, SLposition )

    PCS_HN_predicted = [PCS_H_predicted,PCS_N_predicted[:, 1]]
    '''

#reslist = getResidueNumbers(3) # Change parameter to modify input type
#reslist = np.array(range(1, 76))
""" reslist = np.array(range(0, 80))

susceptibility_tensor = np.array([-0.2073, 0.2761, 0.1499, -0.1655, 0.3134]) # Change parameter to modify input type
SLposition = np.array([57.7920, -86.4502, -2.8023])

pdb_filename = "./1d3z_mod1.pdb"

PCS_H_predicted = predictPCS(susceptibility_tensor, reslist, pdb_filename, SLposition)
print ("SLposition: {}".format(SLposition))
print ("Susceptibility Tensor (tensorvector): {}".format(SLposition))
print ("PCS H Predicted:")
print (PCS_H_predicted) """


"""
---> 3-parameter fit trials results:
[[ 6.37638459e+01 -9.16307520e+01 -1.14260122e+01  1.67616351e-02
   6.55275792e+00]
 [ 6.68584016e+01 -9.24745747e+01 -1.45053396e+01  1.07256452e-01
   1.10579060e+01]
 [ 9.19704404e+01 -6.38456631e+01 -4.11570981e+01  9.52314503e-01
   9.36516783e+01]
 [ 6.37638459e+01 -9.16307521e+01 -1.14260122e+01  1.67616351e-02
   6.55275792e+00]
 [ 1.59791104e+05  2.29536481e+07  3.73760070e+06  1.35108901e+00
   1.09325036e+13]
 [ 6.53649035e+01 -1.03196576e+02  5.72374476e+00  8.72400757e-01
   1.39315323e+01]
 [ 6.37638459e+01 -9.16307520e+01 -1.14260122e+01  1.67616351e-02
   6.55275791e+00]
 [-4.01750394e+06  2.31675317e+07 -1.06957090e+06  1.20244264e+00
   1.02496241e+13]
 [ 6.54797737e+01 -9.90886256e+01 -5.84618519e-01  5.23140938e-01
   1.08533831e+01]
 [ 6.37638459e+01 -9.16307520e+01 -1.14260122e+01  1.67616351e-02
   6.55275790e+00]
 [ 6.37638459e+01 -9.16307521e+01 -1.14260122e+01  1.67616351e-02
   6.55275793e+00]
 [-1.17457467e+06  3.19582210e+07 -5.31025697e+05  1.20504838e+00
   1.96533284e+13]
 [ 7.08414088e+01 -9.52770212e+01 -2.62896070e+01  4.61578950e-01
   2.96063313e+01]
 [ 6.54795393e+01 -2.01729036e+02  2.86883659e+01  1.06945246e+00
   2.94514692e+02]]

---> 3-Parameter Fit Trials Best Solution:
Best Solution:  [ 63.7638459  -91.630752   -11.42601224]
Chi2min =  0.01676163510744273
Condition number =  6.552757915291881
Eigenvectors: changed sign of 2nd eigenvector to keep right-handedness
Position =  [ 63.7638459  -91.630752   -11.42601224]
X_vector =  [ 0.92716251 -1.73912909 -1.35935075 -3.3503403  -3.20032881]
ΔX  Tensor Eigenvalues =  [ 2.23011875  3.56002264 -5.7901414 ]
Axial component ( ΔX a) =  -8.685212093359606
Rhombic component ( ΔX r) =  -1.3299038923264934  ( x10^-32 m^3)
Rhombicity ( ΔX r/ ΔX a)=  0.1531227882556016
Euler Angles (ZYZ Convention): Alpha = 43.72164 Beta = 41.16442 Gamma = -90.17728
Euler Angles (ZYX Convention): Alpha = -46.41182 Beta = -0.11669 Gamma = -41.16428
Euler Angles (XYZ Convention): Alpha = -31.1447 Beta = 61.59499 Gamma = -38.38674
Elapsed Time: 3.799872398376465 seconds


['---> 8-Parameter Fit']
Direct Hit: Chi2 =  0.01676163510744275
Condition Number =  6.552757916321632
Position =  [ 63.76384587 -91.63075201 -11.42601225]
X-vector =  [ 0.9271625  -1.73912907 -1.35935075 -3.35034029 -3.20032879]
Eigenvectors: changed sign of 2nd eigenvector to keep right-handedness
ΔX  Tensor Eigenvalues =  [ 2.23011874  3.56002263 -5.79014137]
Axial component ( ΔX a) =  -8.685212059180161  Rhombic component ( ΔX r) =  -1.3299038844430209  ( x10^-32 m^3)
Rhombicity ( ΔX r/ ΔX a)=  0.15312278795050593
Euler Angles (ZYZ Convention): Alpha = 43.72164 Beta = 41.16442 Gamma = -90.17728
Euler Angles (ZYX Convention): Alpha = -46.41182 Beta = -0.11669 Gamma = -41.16428
Euler Angles (XYZ Convention): Alpha = -31.1447 Beta = 61.59499 Gamma = -38.38674
Elapsed Time: 16.32589101791382 seconds
"""

#XYZ: [[ 0.19429785  0.69044479 -3.51941577  3.06871542  1.77046453]]
#ZYX: [[ 2.23248317 -5.79010521 -0.05621226 -0.01425687 -0.01066854]]
#ZYZ: [[-0.88772155 -2.04006977 -0.48829711 -0.45051906  4.28767003]] 

