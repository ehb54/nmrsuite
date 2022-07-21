import numpy as np
from predict_scripts.utils import findstr
#from utils import findstr

def select_at(atnam, at_res, reslst, atlst, atdim, offs = 0):
    reslst = np.asarray(reslst)
    atlst = np.asarray(atlst)
    atnampos = list(range(6+offs, 7+offs+atdim-1))
    nat = atnam.shape[0]
    indsel = np.full((nat, 1), 0)
    isel = 1
    for i in range(nat): 
        if reslst.size == 0:
            if findstr(atnam[i, atnampos], atlst).size > 0:
                indsel[isel] = i
                isel = isel + 1
        else:
            if np.nonzero(at_res[i, 1] == reslst)[0].size > 0: 
                if atlst.size == 0:
                    indsel[isel] = i
                    isel = isel + 1
                else:
                    if findstr(atnam[i, atnampos], atlst).size > 0:
                        indsel[isel] = i
                        isel = isel + 1
    sel = indsel[~np.isnan(indsel)]
    return sel
