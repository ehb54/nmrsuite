# Inputs: 
# Atom Coordinates
# List of Residues
# Type of Atoms
# Magnetic Susceptibility Tensor
# Frequency
# Order Parameter

# Output: PCS values for each atom/residue (numerical and also plotted as a function of residue number)

import numpy as np
import math as m
from predict_scripts.components2tensor import components2tensor
from predict_scripts.readpdb import readpdb
from predict_scripts.tensor2rdc import tensor2rdc
'''
from components2tensor import components2tensor
from readpdb import readpdb
from tensor2rdc import tensor2rdc
'''

#from helperFunctions import getResidueNumbers, getSusceptibilityTensor

def predictRDC (susceptibility_tensor, reslist, pdb_filename, freq, atom1 = "H", atom2 = "N", Sorder = 0.9, pdb_model = 0, chainID = "A", T = 298):
    #this is a script to use the results of PCS analysis
    #in order to predict RDCs for user-defined residues
    #modified df-jun-2020
    #modified mc-jun-2020
    #modified mc-jul-2021


    #pdb_filename = 'D:\Programs\PDB\1d3z_mod1.pdb'
    #pdb_filename = 'D:\Dropbox\MyMatlab\altens\PCS\PDBs\fixA_Bmnode_1.pdb'
    #pdb_filename = 'D:\Dropbox\MyMatlab\altens\PCS\PDBs\fixB_Amnode_1.pdb'
    #pdb_filename = 'D:\Dropbox\MyMatlab\altens\PCS\PDBs\1D3Z_on_1AAR.pdb'
    #pdb_filename = 'D:\Dropbox\MyMatlab\altens\PCS\PDBs\1AAR_hydr_on_1D3Z_byA_ZonB.pdb'
    #pdb_filename = 'Z:\Dropbox\MyMatlab\PCS\1d3z_mod1.pdb'


    #reslist = np.array(np.range(1, 76))

    #reslist = rdc_800(:,1)
    #reslist = rlist_ub2b
    #reslist = rlist_ub2a
    #reslist = []

    #convert output of PCS analysis into SL coordinates and the tensor

    #---------------------------------
    #---------------------------------
    result_8par = np.array([]) # MODIFY THIS LINE to specify the output you want to use
    #---------------------------------
    #---------------------------------

    result_Xpar = np.copy(result_8par)

    #SLposition = result_Xpar[0:3]      #  we don't care about this for RDCs
    
    #tensorvector = np.transpose(result_Xpar[3:8]) # leaving this to test
    tensorvector = susceptibility_tensor

    #read in PDB and extract proton coords for selected residues
    #vNH = pdb2nh(pdb_filename,reslist,pdb_model)

    #read coordinates of all amide H-N atoms 
    #reslist = np.array([2, 3])
    [coorHall,atnamHall,at_resHall] = readpdb(pdb_filename,reslist,np.array([atom1]),pdb_model,chainID)
    [coorNall,atnamNall,at_resNall] = readpdb(pdb_filename,reslist,np.array([atom2]),pdb_model,chainID)
    coorH_size = np.size(coorHall)
    coorN_size = np.size(coorNall)
    if (coorH_size > coorN_size):
        [coorHall,atnamHall,at_resHall] = readpdb(pdb_filename,at_resNall,np.array([atom1]),pdb_model,chainID)
    elif (coorN_size > coorH_size):
        [coorNall,atnamNall,at_resNall] = readpdb(pdb_filename,at_resHall,np.array([atom2]),pdb_model,chainID)
    #read in only those amide Ns that have Hs attached
    #at_resHall = np.array([2, 3])
    #[coorNall,atnamNall,at_resNall] = readpdb(pdb_filename,at_resHall,np.array([atom2]),pdb_model,chainID)
    
    #build NH vectors
    try:
        vectNH = np.subtract(coorHall,coorNall)
    except:
        raise Exception ("The residue range provided goes out of bounds. Please check it.")
    #normalize the BN vectors and combine with the residue numbers
    for i in range (np.shape(vectNH)[0]):
        vectNH[i] = vectNH[i]/m.sqrt(np.sum(vectNH[i] ** 2))

    #normNH = sqrt(diag(vectNH*transpose(vectNH)))
    vNH = np.column_stack((at_resHall, vectNH))


    #now predict PCSs using this information
    #RDC_pred = tensor2rdc(vNH,tensorvector) 
    return (tensor2rdc(vNH,tensorvector,Sorder,freq,T))


    # Check the slfit.py code for column vector 

""" reslist = np.array(range(6, 11))
susceptibility_tensor = np.array([-0.2073, 0.2761, 0.1499, -0.1655, 0.3134]) # Change parameter to modify input type
pdb_model = 1                  #read in first model
freq = 600   #1H frequency
pdb_filename = "./1d3z_mod1.pdb"

print ("Predicted RDCs: \n{}".format(predictRDC (susceptibility_tensor, reslist, pdb_filename, freq))) """
