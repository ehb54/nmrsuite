import numpy as np
from writepdb import writepdb

import platform, shutil

def SL_add2pdb(old_file, new_file, SL_coor, SL_at_res = [9999, 999]):
    [nat, dummy] = SL_coor.shape
    [nat1, dummy1] = SL_at_res.shape
    if nat != nat1:
        print ('!!!size mismatch of SL_coor and SL_at_res!!! skipped write to PDB')

    a = '/'
    b = '\\'
    
    if platform.system() == "Windows":
        old_file = old_file.replace(b, a)
        new_file = new_file.replace(b, a)
        try:
            shutil.copyfile (old_file, new_file)
            s = 0
        except Exception as error:
            s = 1
            w = error
    else:
        a = "\\"
        b = "/"
        old_file = old_file.replace(a, b)
        new_file = new_file.replace(a, b)
        try:
            shutil.copyfile (old_file, new_file)
            s = 0
        except Exception as error:
            s = 1
            w = error

    if s != 0:
        print (w)
    SL_atnam = []

    for i in range(nat):
        atnum = str(SL_at_res[i, 0])
        resnum = str(SL_at_res[i, 1])
        SL_atnam = [SL_atnam, (4 - len(atnum)) * "\n", atnum, '  S   CYS  ', (4 - len(resnum)) * "\n", resnum]
    
    fid = writepdb(SL_coor, SL_at_res, SL_atnam, 'a', new_file, 0, [], 0)
    if s == 0 and fid != -1:
        print(['data copied/appended to file ', new_file])
    return [s, w]