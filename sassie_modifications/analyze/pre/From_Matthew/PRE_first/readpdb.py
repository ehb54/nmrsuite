import numpy as np
from utils import hard_replace, readasci

def readpdb(fname, reslst = np.array([]), atlst = np.array([]), model = 1):
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
        if ("".join(pdb[i, 0:3]) == "TER" or "".join(pdb[i, 0:3]) == "END" or "".join(pdb[i, 1]) == "."):
            if termflag == 0: 
                nmod = nmod+1
                termflag = 1
                if nmod > model:
                    break
        else:
            termflag = 0
            if nmod == model:
                if 'ATOM' == "".join(pdb[i, 0:4]):
                    if reslst.size > 0:
                        if np.nonzero(reslst == int(pdb[i, 22:26]))[0].size > 0: 
                            select_at[isel] = i
                            isel = isel + 1
                    else:
                        select_at[isel] = i
                        isel = isel + 1

    sel = list(select_at[0: isel].flatten())[1:]
    print (sel[-1])
    nat = len(sel)

    if nat == 0:
        raise Exception('no atoms found!!! wrong filename or model')
    print ("{} atoms read in".format(str(nat)))
    atnam = pdb[sel, 7:26]
    coor = np.full((nat, 3), 0.0)
    at_res = np.full((nat, 2), 0.0)

    #Note: This area can probably be optimized (the replace "H" with "0" method, for example, is very slow)
    new_pdb = hard_replace(pdb, "H", "0")
    new_pdb = hard_replace(new_pdb, " ", "0")
    at_res[:, 0] = list(map(lambda x: "".join(x), new_pdb[sel, 6:11]))
    at_res[:, 1] = list(map(lambda x: "".join(x), new_pdb[sel, 22:26]))
    unstripped_list = [x.split() for x in list(map(lambda x: "".join(x), pdb[sel, 31:54]))]
    coor[:, 0:3] = [[float(y.strip()) for y in x] for x in unstripped_list]

    return [coor, atnam, at_res]