import numpy as np

from scipy.cluster.hierarchy import dendrogram, fcluster, from_mlab_linkage
import matplotlib.pyplot as plt

import os


tree = np.array([[ 3.        ,  4.        ,  5.03302994,  2.        ],
       [ 1.        ,  2.        ,  5.226462  ,  2.        ],
       [ 5.        ,  6.        ,  7.08738055,  4.        ],
       [ 0.        ,  7.        , 12.47809204,  5.        ]])


run_directory = "result_of_maxent"

keep68 = np.loadtxt(os.path.join(run_directory, "cluster.txt"))
index68 = np.loadtxt(os.path.join(run_directory, "index.txt")).astype(int)

xsol = np.loadtxt(os.path.join(run_directory, "x.txt"))[:, index68]

sortrows = lambda matrix: matrix[np.lexsort(matrix.T[::-1])] # Equivalent of Matlab's sortrows

""" 
list_to_int = lambda input_list: list(map(int, input_list))
list_to_str = lambda input_list: list(map(str, input_list))
increment_list = lambda input_list: list(map(lambda element: element+1), input_list) """
#tree = tree[:, :-1] # removing last column to match Matlab result
#tree[:, [0,1,3]] = tree[:, [0,1,3]] + 1 # adding 1 to indices to match Matlab convention

#dn = hierarchy.dendrogram(tree)

#tree = from_mlab_linkage([[4.0000, 5.0000, 5.0330], [2.0000, 3.0000, 5.2265], [6.0000, 7.0000, 7.0874], [1.0000, 8.0000, 12.4781]])
#print (tree)

T = fcluster(Z=tree, t=10, criterion="distance") # the RMSD clustering is chose to be 4A in here

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


#T[:,2] = keep68
""" x = from_mlab_linkage([[4.0000, 5.0000, 5.0330], [2.0000, 3.0000, 5.2265], [6.0000, 7.0000, 7.0874], [1.0000, 8.0000, 12.4781]])

print (x)
T = fcluster(Z=x, t=4, criterion="distance")
print (T) """

""" fig, ax = plt.subplots(figsize=(20, 10))
dn = hierarchy.dendrogram(Z=tree, p=0)

#dn["ivl"] = list(map(lambda element: str(int(element) + 1), dn["ivl"]))
#dn["leaves"] = list(map(lambda element: element + 1, dn["leaves"]))

ax = dn
print (ax) """

#dendrogram(Z, leaf_rotation=90, leaf_font_size=8, labels=df.index)


#plt.show()
""" import plotly.figure_factory as ff
fig = ff.create_dendrogram(X=tree)
fig.show() """