#!/opt/miniconda3/bin/python

import numpy as np
from utils import hard_replace, readasci, str2num
from at_select import at_select
'''
from predict_scripts.utils import hard_replace, readasci, str2num
from predict_scripts.at_select import at_select
'''

def readpdb(fname, reslst = np.array([]), atlst = np.array([]), model = 1, chainID = "A", off = 0):
    reslst = np.asarray(reslst)
    atlst = np.asarray(atlst)
    pdb = np.asarray(readasci(fname)) 
    print ('got the {} data set,  analyzing...'.format(fname))
    nlin = pdb.shape[0]
    select_at = np.full((nlin, 1), 0)
    isel = 1
    nmod = 1
    termflag = 0
    for i in range(nlin): 
        if ("".join(pdb[i, 0:3]) == "TER" or "".join(pdb[i, 0:3]) == "END" or "".join(pdb[i, 0]) == "."):
            if termflag == 0: 
                nmod = nmod+1
                termflag = 1
                if nmod > model:
                    break
        else:
            termflag = 0
            if nmod == model:
                if "".join(pdb[i, 0:4]) == "ATOM" and "".join(pdb[i, 21]) == chainID:
                    if reslst.size > 0:
                        number = int("".join(pdb[i, 22:26]).strip())
                        if np.nonzero(reslst == number)[0].size > 0: 
                            select_at[isel] = i
                            isel = isel + 1
                    else:
                        select_at[isel] = i
                        isel = isel + 1

    sel = list(select_at[0: isel].flatten())[1:]
    nat = len(sel)

    if nat == 0:
        raise Exception('no atoms found!!! wrong filename or model')
    print ("{} atoms read in".format(str(nat)))
    atnam = pdb[sel, 7:26]
    coor = np.full((nat, 3), 0.0)
    at_res = np.full((nat, 2), 0.0)

    #Note: This area can probably be optimized (the replace "H" with "0" method, for example, is very slow)
    new_pdb = hard_replace(pdb, "H", "0")
    at_res[:, 0] = str2num(new_pdb[sel, 6:11], integer = True)
    at_res[:, 1] = str2num(new_pdb[sel, 22:26], integer = True)
    coor[:, 0] = str2num(new_pdb[sel, 30:38])
    coor[:, 1] = str2num(new_pdb[sel, 38:46])
    coor[:, 2] = str2num(new_pdb[sel, 46:54])

    if (atlst.size > 0):
        rlist = np.array((range(int(min(at_res[:, 1])), int(max(at_res[:, 1]))+1)))
        selat = at_select(atnam,at_res,rlist,atlst,off).astype(int)
        coor = coor[selat,:]
        at_res = at_res[selat,:]
        atnam = atnam[selat,:]

    return [coor, atnam, at_res]

#readpdb("mnodesS00001.pdb")

"H     0    NH   11"

#print (readpdb("mnodesS00001.pdb", np.array([2, 4]), np.array(["H ,N "]), 1, 0, "B"))
