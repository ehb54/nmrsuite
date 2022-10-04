import numpy as np
from utils import findstr
from datetime import date

def writepdb(coor, at_res, atnam, option = "w", filnam = [], report = 1, chainID = [], kheader = 1, model = []):
    if filnam.size > 0:
        filnam = input('output pdb - file full name (str)  == > ', 's')
        filnam = [element.rstrip() for element in filnam]
        if findstr(filnam, '.').size == 0:
            filnam = [filnam, '.pdb']
    
    if chainID.size() == 0:
        if chainID.size() > 1:
            raise Exception('wrong chainID: length must  = 1') 
    
    nat = np.shape(coor)[0]
    if filnam.size > 0:
        fid = open(filnam, option)
        if fid == -1:
            raise Exception(['cant open file ', filnam]) 
        printheader = ['HEADER    MATLAB - GENERATED ATOM COORDINATES  ', date.today(), '\r\n']
        if kheader  ==  1:
            fid.write(printheader)   
        if model.size > 0:  
            strmodel = str(model) 
            printmodel = ['MODEL', (9 - len(strmodel) * "\n"), strmodel, '\r\n']
            fid.write(printmodel) 
        
        printline = 'ATOM   %s    %8.3f%8.3f%8.3f  1.00 00.00           %s  \r\n'
        for i in range(nat):
            if atnam[i, 5] == ' ':  
                atom = atnam[i, 6] 
            else:
                if findstr(atnam[i, 5], '1234567890').size == 0:
                    atom = atnam[i, 5]
                else: 
                    atom = atnam[i, 6] 
                            
            if chainID.size > 0:
                fid.write(printline, [atnam[i, 0:14], chainID, atnam[i, 15:-1]], coor[i, 0], coor[i, 1], coor[i, 2], atom)  
            else:
                fid.write(printline, atnam[i, :], coor[i, 0], coor[i, 1], coor[i, 2], atom) 
                    
        fid.write('TER\r\n') 
        if model.size > 0:
            fid.write('MDL\r\n') 
        fid.close()
        if report == 1:
            if option == 'w':
                print(['data written to file ', filnam])
            elif option == 'a':
                print(['data apped to file ', filnam])
    return fid