#!/opt/miniconda3/bin/python
from Bio import PDB
import numpy as np

# New version of readpdb coded 
# Written by MC (along with all other Python programs in "Predict", "MaxEnt", and "SES")

def readpdb(fname, reslst = np.array([]), atlst = np.array(["H"]), model_num = 0, chainID = "A"): 
    parser = PDB.PDBParser()
    struct = parser.get_structure('1abz',fname)
    reslstIsNone = np.size(reslst) == 0
    coordinates = np.array([])
    atnam = np.array([])
    at_res = np.array([])
    coordinatesIsEmpty = True # Variable that helps construct the coordinates array
    for model in struct:
        if (model.id == model_num):
            for chain in model:
                if (chain.id == chainID):
                    for residue in chain:
                        if reslstIsNone or np.isin(residue.id[1], reslst):
                            for atom in residue:
                                if (np.isin(atom.id, atlst)):
                                    XYZ = atom.get_coord()
                                    if coordinatesIsEmpty:
                                        coordinates = np.reshape(XYZ, (1, XYZ.size))
                                        coordinatesIsEmpty = False
                                    else:
                                        coordinates = np.vstack([coordinates, XYZ])
                                    atnam = np.append(atnam, atom.id)
                                    at_res = np.append(at_res, int(residue.id[1]))
                                    
    return [coordinates, atnam, at_res]

#print (readpdb("mnodesS00001.pdb", np.array([2, 4]), np.array(["H", "N"]), 0, "B"))
