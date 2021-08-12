#!/opt/miniconda3/bin/python
#This script uses the definition in Jacs paper on how to decide which
#structures to keep. In the end you have a vector "keep68" that tells you
#the important structures you should be looking at. After here you can
#still cluster the structures in keep68 by rmsd using the script
#ClusterByRMSD.m

import numpy as np
import os

def filter_by_weight (run_directory):
    # Load all variables from their associated files in the 'result_of_maxent' folder
    A = np.loadtxt(os.path.join(run_directory, "A.txt"))
    lambda_ = np.loadtxt(os.path.join(run_directory, "lambda.txt")) #'lambda' is a reserved keyword in Python, so a trailing underscore was added to the variable name per the style guide
    x = np.loadtxt(os.path.join(run_directory, "x.txt"))
    y = np.loadtxt(os.path.join(run_directory, "y.txt"))

    with open(os.path.join(run_directory, "index.txt"), "r") as f: index68 = int(f.read()) # Load index value from its file

    std68 = np.std(x[:,index68]) 
    mean68 = np.mean(x[:,index68]) 

    def clustering(avg, sdevi, maxent):
        keep = np.array([])
        limit = avg + 2 * sdevi
        
        for i in range (len(maxent)):
            if (maxent[i] > limit):
                keep = np.append(keep, i)

        return keep

    keep68 = clustering(mean68, std68, x[:,index68]) + 1
    np.savetxt(os.path.join(run_directory, "cluster.txt"), keep68)

#filter_by_weight("result_of_maxent")
