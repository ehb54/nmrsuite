#!/opt/miniconda3/bin/python
#This script compares pdb files and write pairwise RMSD in a table: total_structs
#Here it's assumed that the name of the pdb files are struct# with # = 00001,00002 and so on.
#It's also assumed that the pdb files are written in a certain way like in
#the ensemble with name struct#. Functions readpdb and readpdb1 might need to have
#inputs adjusted if the pdb file is a little bit different. Like in
#ensemble mnodes. 

#Note: File is not yet complete

import numpy as np
import os

from scipy.cluster.hierarchy import dendrogram, linkage
from pathlib import Path

from read_utils import readpdb, readpdb1, superimpose


def ClusterByRMSD(run_directory):
    pdbfilenameA = r"C:\Users\Owner\Documents\Fushman Internship\MaxEntropy\Python\STRUCT_superimposed_A\struct00001.pdb"
    pdbfilenameB = r"C:\Users\Owner\Documents\Fushman Internship\MaxEntropy\Python\STRUCT_superimposed_A\struct00001.pdb"
    #pdbfilenameA = 'G:\My Drive\University of Maryland\Summer Project 2020\Ensembles\MNODES_superimposed_B\mnodesS00001.pdb'
    #pdbfilenameB = 'G:\My Drive\University of Maryland\Summer Project 2020\Ensembles\MNODES_superimposed_B\mnodesS00001.pdb'

    #struct = strfind(pdbfilenameB,'structS00001') + 11 
    #struct = strfind(pdbfilenameB,'mnodesS00001') + 11

    keep68 = np.loadtxt(os.path.join(run_directory, "cluster.txt"))

    numstruct = np.size(keep68) #TOTAL NUMBER OF STRCTURES IT'S GOING TO BE USED

    off = 1 #for STRUCT, 0 for MNODES
    #read in first domain
    [coorA1,atnamA1,at_resA1] = readpdb(pdbfilenameA,np.array([]),np.array([]),1,off,'A')
    #read in second domain
    [coorA2,atnamA2,at_resA2] = readpdb(pdbfilenameA,np.array([]),np.array([]),1,off,'B')

    [coorB1,atnamB1,at_resB1] = readpdb(pdbfilenameB,np.array([]),np.array([]),1,off,'A')
    [coorB2,atnamB2,at_resB2] = readpdb(pdbfilenameB,np.array([]),np.array([]),1,off,'B')

    total_structs = np.zeros((np.size(keep68), np.size(keep68)))

    for m in range (numstruct):
        i = keep68[m] 
        pdbfilenameA = r"C:\Users\Owner\Documents\Fushman Internship\MaxEntropy\Python\STRUCT_superimposed_A\struct00001.pdb"
        pdbfilenameA = pdbfilenameA.replace("00001", str(int(i)).zfill(5)) # zfill part adds 0s to make the string 5 characters long
                            
        coorA1 = readpdb1(pdbfilenameA,np.array([]),np.array([]),1,at_resA1,at_resA1) #FOR STRUCTS (underscores are used because the other return values are unnecessary)
        coorA2 = readpdb1(pdbfilenameA,np.array([]),np.array([]),1,at_resA2,at_resA2) #FOR STRUCTS (underscores are used because the other return values are unnecessary)
        #[coorA1,atnamA1,at_resA1] = readpdb(pdbfilenameA,[],[],1,off,'A') #FOR MNODES
        #[coorA2,atnamA2,at_resA2] = readpdb(pdbfilenameA,[],[],1,off,'B') #FOR MNODES
        
        a = 0
        for s in range (numstruct):
            j = keep68[s]
            #total_structs(m,1) = i
            
            if j != i:
                pdbfilenameB = r"C:\Users\Owner\Documents\Fushman Internship\MaxEntropy\Python\STRUCT_superimposed_A\struct00001.pdb"
                pdbfilenameB = pdbfilenameB.replace("00001", str(int(j)).zfill(5))

                #The other picture
                #coorB1 = readpdb1(pdbfilenameB,[],[],1,at_resB1,at_resB1) #FOR STRUCTS
                #coorB2 = readpdb1(pdbfilenameB,[],[],1,at_resB2,at_resB2)
                [coorB1,atnamB1,at_resB1] = readpdb(pdbfilenameB,np.array([]),np.array([]),1,off,'A')#FOR MNODES
                [coorB2,atnamB2,at_resB2] = readpdb(pdbfilenameB,np.array([]),np.array([]),1,off,'B')

                #list of residues to superimpose
                reslst = np.array([])

                [coorBrot,Rotmat,rmsd,coor2s] = superimpose(np.vstack((coorA1, coorA2)),np.vstack((atnamA1, atnamA2)),np.vstack((at_resA1, at_resA2)),np.vstack((coorB1,coorB2)),np.vstack((atnamB1, atnamB2)),np.vstack((at_resB1, at_resB2)),reslst,reslst,np.array([]),np.vstack((coorB1, coorB2)))

                total_structs[m,a] = rmsd
                a = a+1
        #save total_structs_rmsd_STRUCT_Clust_2A_maxweight.mat total_structs '-v7.3'

    #This loop is just to put the vector total_structs that has the pairwise
    #RMSD as a matrix with zeros in the diagonal. 
    for i in range (np.size(total_structs[:,0])):
        a = total_structs[i, np.nonzero(total_structs[i,:])][0]
        total_structs[i,i] = 0
        total_structs[i,i+1:] = a[i:]

    rows, columns = np.shape(total_structs) 
    v = np.array([]) #v is going to be the vector that contains all the pairwise RMSD
    for i in range(rows-1):
        v = np.append(v, total_structs[i, i+1:columns]) 

    print ("Done for now...")
    tree = linkage(v,'average') #ths function linkage does hierarchal clustering

    dendrogram(tree,0) #this plots the tree with the clusterings

    import plotly.figure_factory as ff
    fig = ff.create_dendrogram(tree)
    fig.show()

    """T = cluster(tree,'cutoff',4,'Criterion','distance') #the RMSD clustering is chose to be 4A in here
    T(:,2) = 1:length(T) 
    T(:,3) = keep68
    T = sortrows(T)

    xsol = result_of_maxent.x(:,index68)

    #Here I choose the representative of the cluster. The way I do, the
    #representative is the structure with the highest initial weight (assigned by max entropy)
    #The final weight of the representative is the summation of weights for everyone in the cluster
    for i = 1:length(unique(T(:,1)))
        aux = find(T(:,1) == i)
        weights(i) = sum(xsol(T(aux,3))) 
        structs(i) = find(xsol == max(xsol(T(aux,3)))) 
    end

    save weights.mat weights
    save structs.mat structs """


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = == 
""" def superimpose(coor1,atnam1,at_res1,coor2,atnam2,at_res2,reslst1,reslst2,atsel,coor3):
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
    sel1 = at_select(atnam1,at_res1,reslst1,atsel,3)
    #sel1 = at_select(atnam1,at_res1,reslst1,atsel,2,-1)
    sel2 = at_select(atnam2,at_res2,reslst2,atsel,3)
    #sel2 = at_select(atnam2,at_res2,reslst2,atsel,2,-1)

    #length(sel1)
    #length(sel2)
    if np.size(sel1) != np.size(sel2):
        raise Exception('sets of atom coordinates are different!')

    #matching atom names
    sel2m = sel2
    for i in range (np.size(reslst1)):
    ind1 = find(at_res1(sel1(:),2) == reslst1[i])
    ind2 = find(at_res2(sel2(:),2) == reslst2[i])
    if (np.size(ind1)!=0 and np.size(ind2)!=0):
        if (np.size(ind1) != np.size(ind2)):
                raise Exception(f'atoms in residues {reslst1[i]} and {reslst2[i]} are different!') 
        #[atnam1(sel1(ind1),:),atnam2(sel2(ind2),:)]
        ind2m = ind2
        for j in range(np.size(ind1)):
            for k in range(np.size(ind2)):
                if strcmp(atnam1(sel1(ind1(jj)),6:9),atnam2(sel2(ind2(kk)),6:9))   
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
    if exist('coor3'):
        [coor3s,Rotmat,coor2s,rmsd] = rotfit(coor1(sel1,:),coor2(sel2m,:),coor3)
    else:
        [coor3s,Rotmat,coor2s,rmsd] = rotfit(coor1(sel1,:),coor2(sel2m,:))
    #disp(['rmsd =  ',num2str(rmsd)])
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    return function [coor3s,Rotmat,rmsd,coor2s] """


""" def rotfit (a, b1, b2 = np.array([]))
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
    cma = np.mean(a)
    cmb1 = np.mean(b1)
    # shift to put c-o-mass in origin
    _as = a - np.ones((nat,) 1)*cma # added _ to as variable because it is a reserved keyword in Python
    b1s = b1 - np.ones((nat,) 1)*cmb1
    if (sets == 2):
        b2s = b2 - np.ones((natb2,) 1)*cmb1
    detind = 1
    u = np.zeros(3,3)
    omega = np.zeros(6,6)
    u = _as.H * b1s
    if (det(u) < 0):
        raise Exception ('Determinant negative!!!') # matlab code just had 'disp' here
        detind = 0  # negative-det(u)-records
    omega[0:3,3:6] = u				#build omega-matrix
    omega[3:6,0:3] = u.H
    [V,D] = np.linalg.eig(omega)              	#eigenvalues of omega
    rot = np.zeros(3,3) 
    H1 = np.zeros(3,3)
    K1 = np.zeros(3,3)
    r1 = np.zeros(3,3)
    for k = in range (6):
        K1 = ones((3,1))*V(1:3,k)'
        H1 = V(4:6,k)*ones((1,3))
        r1 = K1.*H1
        rot = rot+sign(D(k,k))*r1 
    
    b1rot = b1s*rot				#(rot'*b1s')'	
    #calculate RMSDs
    diff = _as-b1rot
    rmsd = np.sqrt(np.mean(np.sum((diff**2))))
    #shift <b>'s to COM of <a>
    b1rot = b1rot+np.ones((nat,1))*cma		# shift back to c-o-mass of <a>
    if (sets == 2): 
        b2rot = b2s*rot			#(rot'*b2s')'
        b2rot = b2rot+np.ones((natb2,1))*cma 	# shift back to c-o-mass of <a>
    else:
        b2rot = b1rot
    return [b2rot,rot,b1rot,rmsd] """

ClusterByRMSD("result_of_maxent")