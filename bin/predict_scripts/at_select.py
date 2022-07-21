import numpy as np

from predict_scripts.utils import findstr
#from utils import findstr

def at_select(atnam, at_res, reslst, atlst, atdim = 2, offs = 0): 
    #---------------------------------------------------
    #  mc-jul-20 df-aug-03    df-dec-97
    #	select atoms according to atnam & reslst
    #   standard lists: 'all', 'heavy', 'bb', 'nh'
    #   atdim is required if atlist in not a standard one
    #---------------------------------------------------
    
    atsel_flag = 1

    lower_atlst = np.char.lower(atlst)

    if (lower_atlst == "heavy"):
        at = 'NCOS'
        atlen = 1
    elif (lower_atlst == "bb"):
        at = 'N CAC O '
        atlen = 2
    elif (lower_atlst == "nh"):
        at = 'HN'
        atlen = 2
    elif (lower_atlst == "all" or atlst.size == 0):
        #get all atoms
        #at = 'HNCOS'atlen = 1 
        atsel_flag = 0
        atlen = 1
    else:
        atlen = atdim
        at = atlst

    atnampos = list(range(6+offs, 7+offs+atlen-1))
    nat = atnam.shape[0]
    indsel = np.full((nat, 1), float('nan'))
    isel = 1

    for i in range(nat): 
        if reslst.size == 0:
            if (atsel_flag == 0):
                indsel[i] = i
            else:
                if findstr(atnam[i, atnampos], at).size > 0:
                    indsel[i] = i
        else:
            if np.nonzero(at_res[i, 1] == reslst)[0].size > 0: 
                if (atsel_flag == 0):
                    indsel[i] = i
                else:
                    if findstr(atnam[i, atnampos], at).size > 0:
                        indsel[i] = i

    sel = indsel[~np.isnan(indsel)]
    return sel
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#8 26 46 61 79 90

#i = 6 (go to 5), ii = 7
