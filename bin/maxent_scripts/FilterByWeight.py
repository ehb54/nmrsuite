#!/opt/miniconda3/bin/python
#This script uses the definition in Jacs paper on how to decide which
#structures to keep. In the end you have a vector "keep68" that tells you
#the important structures you should be looking at. After here you can
#still cluster the structures in keep68 by rmsd using the script
#ClusterByRMSD.m

import numpy as np
import os

#from L_curveServerVersion import L_curve_for_best_lambda

def clustering(avg, sdevi, maxent):
    keep = np.array([])
    limit = avg + 2 * sdevi

    for i in range (len(maxent)):
        if (maxent[i] > limit):
            keep = np.append(keep, i)

    # Do the loop again if keep is empty
    if (np.size(keep) == 0):
        limit = 2 * sdevi
        for i in range (len(maxent)):
            if (maxent[i] > limit):
                keep = np.append(keep, i)
    #return str(keep) + "\n" + str(sdevi) + "\n" + str(maxent)
    return keep

def filter_by_weight (run_directory):
    # Load all variables from their associated files in the 'result_of_maxent' folder
    A = np.loadtxt(os.path.join(run_directory, "A.txt"))
    lambda_ = np.loadtxt(os.path.join(run_directory, "lambda.txt")) #'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
    x = np.loadtxt(os.path.join(run_directory, "weights_for_all_lambdas.txt"))
    y = np.loadtxt(os.path.join(run_directory, "data.txt"))

    with open(os.path.join(run_directory, "index.txt"), "r") as f: index68 = int(f.read()) # Load index value from its file

    std68 = np.std(x[:,index68]) 
    mean68 = np.mean(x[:,index68])

    keep68 = clustering(mean68, std68, x[:,index68])
    keep68_weights = np.transpose(x[np.asarray(keep68, int)])[index68]
    keep68_weights = np.divide(keep68_weights, np.sum(keep68_weights)) # Normalize
    
    # Corresponding weights
    np.savetxt(os.path.join(run_directory, "cluster.txt"), np.transpose(np.vstack((keep68 + 1, keep68_weights)))) # Adds 1 to match 1-indexing (structures start with 1)

#filter_by_weight("result_of_maxent43")
