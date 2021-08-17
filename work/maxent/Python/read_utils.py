#!/opt/miniconda3/bin/python
import numpy as np
from utils import findstr, hard_replace, str2num

def matlabeig (matrix):
    w, v = np.linalg.eig(matrix)
    d = np.diag(w)
    return v, d

def at_select(atnam, at_res, reslst, atlst, atdim = 2, offs = 0): 
    #---------------------------------------------------
    #  mc-jul-20 df-aug-03    df-dec-97
    #	select atoms according to atnam & reslst
    #   standard lists: 'all', 'heavy', 'bb', 'nh'
    #   atdim is required if atlist in not a standard one
    #---------------------------------------------------
    
    atsel_flag = 1

    if (np.size(atlst) != 0):
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

def readasci(filename):
    with open(filename) as my_file:
        lines = my_file.readlines()

    length = max(list(map(lambda x: len(list(x)), lines)))
    lines = list(map(lambda x: x.ljust(length), lines))
    lines = list(map(lambda x: list(x), lines))
    text = np.asarray(lines, dtype = "str")
    print (text.shape)
    return text

def readpdb(fname, reslst = np.array([]), atlst = np.array([]), model = 1, off = 0, chainID = "A"): # off and chainID are new inputs that were not in previous versions of readpdb
    reslst = np.asarray(reslst)
    atlst = np.asarray(atlst)
    pdb = np.asarray(readasci(fname))
    #print ('got the {} data set,  analyzing...'.format(fname))
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
                if (("".join(pdb[i, 0:4]) == 'ATOM') and (str(pdb[i, 21]) == chainID)): # If it's a row describing an atom (Note: Chain ID part was not in previous versions of readpdb)
                    if reslst.size > 0:
                        number = int("".join(pdb[i, 22:26]).strip()) # Residue number
                        if np.nonzero(reslst == number)[0].size > 0: # If the residue number of the atom of the current row is found in the residue list
                            select_at[isel] = i
                            isel = isel + 1
                    else:
                        select_at[isel] = i
                        isel = isel + 1

    sel = list(select_at[0: isel].flatten())[1:]
    nat = len(sel)

    if (nat == 0): raise Exception('no atoms found!!! wrong filename or model')

    print ("{} atoms read in".format(str(nat)))
    atnam = pdb[sel, 7:26]
    coor = np.full((nat, 3), 0.0)
    at_res = np.full((nat, 2), 0.0)

    # Note: This area can probably be optimized (the replace "H" with "0" method, for example, is very slow)
    new_pdb = hard_replace(pdb, "H", "0")

    at_res[:, 0] = str2num(new_pdb[sel, 6:11], integer = True)
    at_res[:, 1] = str2num(new_pdb[sel, 22:26], integer = True)
    coor[:, 0] = str2num(new_pdb[sel, 30:38])
    coor[:, 1] = str2num(new_pdb[sel, 38:46])
    coor[:, 2] = str2num(new_pdb[sel, 46:54])

    if (atlst.size > 0):
        rlist = np.array((range(int(min(at_res[:, 1])), int(max(at_res[:, 1]))+1)))
        selat = at_select(atnam,at_res,rlist,atlst,2,off).astype(int) # Two new inputs added here from previous versions of readpdb (2, off)
        coor = coor[selat,:]
        at_res = at_res[selat,:]
        atnam = atnam[selat,:]

        print ("{} atoms selected".format(str(len(selat))))

    return [coor, atnam, at_res]

def readpdb1(fname,reslst,atlst,model,tag1,tag2): #NEW INPUT FLAG
    #-----------------------------------------------------------------
    #   df-nov-97 mc-aug-21
    #	read in a pdb data set 
    #       uses: readasci.m
    #(c) Copyright by David Fushman, U.Maryland
    #-----------------------------------------------------------------
    """ if nargin<4, model=1; 	end		%default: first model
    if nargin<3, atlst=[];  end      %default: all atoms
    if nargin<2, reslst=[]; end		%default: all resid. """

    selat = tag1[:,0]
    sel = tag2[:,0].astype(int)

    pdb=readasci(fname)
    print ('got the {} data set,  analyzing...'.format(fname))
    nlin=np.shape(pdb)[0]
    select_at=np.zeros((nlin,1))
    nmod=1

    nat=np.size(sel)

    if nat==0:
        raise Exception ('no atoms found!!! wrong filename or model')

    print ("{} atoms read in".format(str(nat)))

    #coor=np.zeros((nat,3))
    coor=np.zeros((max(nat, np.max(sel)),3)) # I had to do this because numpy arrays can't dynamically increase in size

    # Note: This area can probably be optimized (the replace "H" with "0" method, for example, is very slow)
    new_pdb = hard_replace(pdb, "H", "0")

    coor[sel-1, 0] = str2num(new_pdb[sel, 30:38])
    coor[sel-1, 1] = str2num(new_pdb[sel, 38:46])
    coor[sel-1, 2] = str2num(new_pdb[sel, 46:54])
    #coor(:,1:3)=str2num(pdb(sel,31:54))

    #select only required atoms 
    if (atlst.size > 0):
        coor=coor[selat,:]
        #disp([num2str(length(selat)),' atoms selected'])
    else:
        coor=coor[sel-1,:]
    
    return coor

def superimpose(coor1,atnam1,at_res1,coor2,atnam2,at_res2,reslst1,reslst2,atsel,coor3=False):
    #---------------------------------------------------
    #   mc-aug-21 (converted to Python)
    #   df-jun-13  fixed the bug in atnam positions
    #   df-dec-09 (takes care of differences in atom order)  
    #   df-aug-03
    #   given a list of residues and atom selections
    #   superimpose coor2 onto coor1, 
    #   if coor3 is given, rotate this set as well
    #   if not, rotated coor2 will be reported as coor3s
    #---------------------------------------------------

    #preliminary selection
    sel1 = at_select(atnam1,at_res1,reslst1,atsel,3).astype(int)
    #sel1 = at_select(atnam1,at_res1,reslst1,atsel,2,-1)
    sel2 = at_select(atnam2,at_res2,reslst2,atsel,3).astype(int)
    #sel2 = at_select(atnam2,at_res2,reslst2,atsel,2,-1)

    #length(sel1)
    #length(sel2)
    if np.size(sel1) != np.size(sel2):
        raise Exception('sets of atom coordinates are different!')

    #matching atom names
    sel2m = sel2
    for i in range (np.size(reslst1)):
        ind1 = np.nonzero(at_res1[sel1[:],1] == reslst1[i])
        ind2 = np.nonzero(at_res2[sel2[:],2] == reslst2[i])

        if (np.size(ind1)!=0 and np.size(ind2)!=0):
            if (np.size(ind1) != np.size(ind2)):
                    raise Exception(f'atoms in residues {reslst1[i]} and {reslst2[i]} are different!') 
            #[atnam1(sel1(ind1),:),atnam2(sel2(ind2),:)]
            ind2m = ind2
            for j in range(np.size(ind1)):
                for k in range(np.size(ind2)):
                    if atnam1[sel1[ind1[j]],5:9] == atnam2[sel2[ind2[k]], 5:9]:
                        ind2m[j] = ind2[k]

            sel2m[ind2] = sel2[ind2m]
            #disp('modif')
            #[atnam1(sel1(ind1),:),atnam2(sel2m(ind2),:)]
    
    #sel1 = at_select(atnam1,at_res1,reslst1,atsel)
    #sel2 = at_select(atnam2,at_res2,reslst2,atsel)
    #length(sel1)
    #length(sel2m)
    #[atnam1(sel1,:),atnam2(sel2,:)]
    #if length(sel1)~ = length(sel2), 
    #    error('sets of atom coordinates are different!')
    #end
    #NEED to coordinate atom names in both sets!!!
    if (not isinstance(coor3, bool)): # If it isn't provided, it's automatically set to False, so this is just checking if it wasn't provided
        [coor3s,Rotmat,coor2s,rmsd] = rotfit(coor1[sel1,:],coor2[sel2m,:],coor3)
    else:
        [coor3s,Rotmat,coor2s,rmsd] = rotfit(coor1[sel1,:],coor2[sel2m,:])
    #disp(['rmsd =  ',num2str(rmsd)])
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    return [coor3s,Rotmat,rmsd,coor2s]


def rotfit (a, b1, b2=np.array([])):
    #----------------------------------------------------
    # mc-aug-21 df-aug-03     df-nov-97		df-jul-97	df-jan-93   
    #    perform a rotation of coord.set <b1>  
    #    in order to superimpose it onto <a>
    #    and then apply this rotation to coord.set <b2>  
    #    as a first step the centers-of-mass (geom.centers)
    #    are superimposed
    #    it is assumed that <b1> is a subset of <b2>
    #
    #    for the algorithm see: 
    #    A.D.McLachlan, J.Mol.Biol.,1979,v.128,49 (appendix)
    #
    #    CALL: [b2rot,rot,b1rot] = mdrotfit(a,b1,b2)
    #          a -  coordinates (number-of-atoms x 3)
    #          b1 - coordinates (number-of-atoms x 3)
    #	   b2 - coordinates (number-of-atoms x 3)
    #    OUTPUT: b2rot, b1rot being the result of a rotation of b2
    #	     rot is a rotation matrix, superimposing 
    #	     <b1> onto <a>
    #    DOESN'T work for number-of-atoms = 1  !!!
    #    Modification of nov97: removed row2coor, mtrx2row
    #----------------------------------------------------
    sets = 2
    if (np.size(b2) == 0):
        sets = 1 #only one set to rotate
    nat = np.shape(a)[0] 
    natb1 = np.shape(b1)[0] 
    if (nat != natb1):
        raise Exception ('different number of atoms in the structures!')
    if (nat == 1):
        raise Exception ('DOESN"T work for number-of-atoms = 1 !!!')
    if (sets == 2):
        natb2 = np.shape(b2)[0]
    cma = np.mean(a, axis=0)
    cmb1 = np.mean(b1, axis=0)
    # shift to put c-o-mass in origin
    _as = a - np.ones((nat, 1))*cma # added _ to as variable because it is a reserved keyword in Python
    b1s = b1 - np.ones((nat, 1))*cmb1
    if (sets == 2):
        b2s = b2 - np.ones((natb2, 1))*cmb1
    detind = 1
    u = np.zeros((3,3))
    omega = np.zeros((6,6))
    u = _as.T.conj() @ b1s
    if (np.linalg.det(u) < 0):
        print ('Determinant negative!!!') # matlab code just had 'disp' here
        detind = 0  # negative-det(u)-records
    omega[0:3,3:6] = u				#build omega-matrix
    omega[3:6,0:3] = u.conj().T
    [V,D] = matlabeig(omega)              	#eigenvalues of omega
    rot = np.zeros((3,3)) 
    H1 = np.zeros((3,3))
    K1 = np.zeros((3,3))
    r1 = np.zeros((3,3))
    for k in range (6):
        K1 = np.ones((3,1))*V[0:3,k].conj().T
        H1 = (np.ones((3,3))*V[3:6,k]).T
        r1 = K1 * H1
        rot = rot+np.sign(D[k,k])*r1 
    
    b1rot = b1s@rot				#(rot'*b1s')'	
    #calculate RMSDs
    diff = _as-b1rot
    rmsd = np.sqrt(np.mean(np.sum((diff**2).conj().T, axis=0).conj().T))
    #shift <b>'s to COM of <a>
    b1rot = b1rot+np.ones((nat,1))*cma		# shift back to c-o-mass of <a>
    if (sets == 2): 
        b2rot = b2s@rot			#(rot'*b2s')'
        b2rot = b2rot+np.ones((natb2,1))*cma 	# shift back to c-o-mass of <a>
    else:
        b2rot = b1rot
    return [b2rot,rot,b1rot,rmsd]