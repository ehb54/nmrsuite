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

from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

from pathlib import Path

from maxent_scripts.ClusterUtils import readpdb, readpdb1, superimpose


sortrows = lambda matrix: matrix[np.lexsort(matrix.T[::-1])] # Equivalent of Matlab's sortrows

def cluster_by_rmsd(run_directory, structure_directory):
    pdbfilenameA = structure_directory
    pdbfilenameB = structure_directory
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
        pdbfilenameA = structure_directory
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
                pdbfilenameB = structure_directory
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

    tree = linkage(v,'average') #ths function linkage does hierarchal clustering

    keep68 = np.loadtxt(os.path.join(run_directory, "cluster.txt"))
    index68 = np.loadtxt(os.path.join(run_directory, "index.txt")).astype(int)

    xsol = np.loadtxt(os.path.join(run_directory, "x.txt"))[:, index68]

    T = fcluster(Z=tree, t=4, criterion="distance") # the RMSD clustering is chose to be 4A in here

    size = np.size(T)
    T_matrix = np.zeros((size, 3))
    T_matrix[:,0] = T
    T_matrix[:,1] = np.arange(size)
    T_matrix[:,2] = keep68

    #T = T_matrix[T_matrix[:,0].argsort()]
    T = sortrows(T_matrix)  # using lambda function

    num_clusters = np.size(np.unique(T[:,0]))

    weights = np.empty(num_clusters)
    structs = np.empty(num_clusters)

    # Here, the representative of the cluster is chosen
    # The representative is the structure with the highest initial weight (assigned by max entropy)
    # The final weight of the representative is the summation of weights for every point in the cluster
    for i in range (num_clusters):
        aux = np.nonzero(T[:,0] == (i+1))[0]
        x = xsol[T[aux,2].astype(int)-1]
        weights[i] = float(np.sum(x))
        structs[i] = float(np.max(x))

    np.savetxt(os.path.join(run_directory, "weights.txt"), weights)
    np.savetxt(os.path.join(run_directory, "structs.txt"), structs)


#ClusterByRMSD("result_of_maxent")