'''
from getNHcoor import getNHcoor
from buildH import buildH
from select_at import select_at
from readpdb import readpdb
'''
from predict_scripts.getNHcoor import getNHcoor
from predict_scripts.buildH import buildH
from predict_scripts.select_at import select_at
from predict_scripts.readpdb import readpdb

import numpy as np

def pdb2nhcoor(fname, reslst = [], model = 0, kNaN = 1, chainID = "A"):
    """ if (len(chainID) > 0):
        [coor, atnam, at_res] = readpdb_chainID(fname, chainID, reslst, [], model) 
    else: """
    [coor, atnam, at_res] = readpdb(fname, reslst, [], model)

    NH_coor = getNHcoor(coor, atnam, at_res, reslst, 0)
    NH_coor = np.asarray(NH_coor)
    if NH_coor.size > 0:
        return NH_coor

    NH_coor = np.asarray(getNHcoor(coor, atnam, at_res, reslst, -1))
    if NH_coor.size > 0:
        return NH_coor

    print('no amide Hs found in PDB,  building NH vectors')

    """ if chainID.size > 0: 
        [coor, atnam, at_res] = np.readpdb_chainID(fname, chainID, [], [], model) 
    else: """ 
    [coor, atnam, at_res] = readpdb(fname, np.array([]), np.array(["H", "N"]), model, chainID)
    #readpdb(fname, reslst = np.array([]), atlst = np.array(["H"]), model_num = 0, chainID = "A"): 
    

    """     print ("ADf")
    print (atnam)
    print (at_res)
    print (reslst)
    print ("ADf") """

    """ strip = lambda x: x.strip()
    at_res = np.vectorize(strip)(at_res) """

    return at_res, reslst
    
    iN = select_at(atnam, at_res, reslst, 'N', 3, -1)
    iC = select_at(atnam, at_res, np.subtract(reslst, 1), 'C', 3, -1)
    iCa = select_at(atnam, at_res, reslst, 'CA', 3, -1)
    nN = iN.shape[0]
    NH_coor = np.full((nN, 7), float('nan'))
    NH_coor[:, 0] = at_res[iN, 1]
    angleCaNH = 121.9947
    vNHlen = 1
    for i in range(nN):
        if ("".join(atnam[iN[i], 10:13].flatten()) != "PRO"):
            indC = np.nonzero(at_res[iC, 1] == at_res[iN[i], 1] - 1)[0]
            if indC.size > 0:
                vCaN = coor[iCa[i], :] - coor[iN[i], :]
                vNC = coor[iN[i], :] - coor[iC[indC], :]
                NH_coor[i, 1:4] = coor[iN[i], :]
                NH_coor[i, 4:7] = coor[iN[i], :] + buildH(vCaN, vNC, angleCaNH) * vNHlen
    if kNaN  ==  1:
        inotNaN = np.nonzero(~np.isnan(NH_coor[:, 1]))[0]
        NH_coor = NH_coor[inotNaN, :]
    return NH_coor
