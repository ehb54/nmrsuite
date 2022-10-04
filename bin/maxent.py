#!/opt/miniconda3/bin/python

import json, sys, time, os
import numpy as np
from io import StringIO ## for Python3
from genapp3 import genapp ## python3
from maxent_scripts.RunMaxEntropy import run_max_entropy
from maxent_scripts.L_curveForBestLambda import L_curve_for_best_lambda
from maxent_scripts.FilterByWeight import filter_by_weight
from maxent_scripts.ClusterByRMSD import cluster_by_rmsd

from scipy.cluster.hierarchy import dendrogram, fcluster, linkage


# function to check if a checkbox was checked by seeing if its id exists in the JSON
is_checked = lambda json_variables, checkbox_id : checkbox_id in json_variables.keys()
sortrows = lambda matrix: matrix[np.lexsort(matrix.T[::-1])] # Equivalent of Matlab's sortrows

def printQuit (string): # Prints a desired string
    output = {}
    output["_textarea"] = str(string)
    print (json.dumps(output))
    quit()

if __name__=='__main__':

        argv_io_string = StringIO(sys.argv[1])
        json_variables = json.load(argv_io_string)

        output_str = ""
        output = {}
        run_directory = json_variables["directory_name"]

        ### initialize the genapp object
        ga = genapp( json_variables )


        """ output = {}
        output['_textarea'] = str(type(A_filename)) + "\n" + str(A_filename)# + "\n" + str(os.path.exists(A_filename)) + "\n" + str(os.path.exists(y_filename))
        print(json.dumps(output)) """

        if is_checked(json_variables, "maxent_checkbox"):
                A_filename = json_variables["matrixfile"][0]
                y_filename = json_variables["datafile"][0]
                lambda_lower = float(json_variables["lambda_lower_bound"])
                lambda_step = float(json_variables["lambda_step_size"])
                lambda_upper = float(json_variables["lambda_upper_bound"])
                sigma = float(json_variables["sigma"]) # For outliers
                run_max_entropy (A_filename, y_filename, lambda_lower, lambda_step, lambda_upper, run_directory)
                #printQuit(L_curve_for_best_lambda(run_directory, sigma))
                line_plot, histogram, scatter_plot = L_curve_for_best_lambda(run_directory, sigma)
                output_str += f"Uploaded maximum entropy result files (A.txt, index.txt, lambda.txt, weights_for_all_lambdas.txt, data.txt) to {run_directory}\n"
                output['lineplot'] = line_plot
                output['histogram'] = histogram
                output['scatterplot'] = scatter_plot
                output_str += f"Created L_curve plots\n"


        if is_checked(json_variables, "filter_by_weight_checkbox"):
                filter_by_weight(run_directory)
                #printQuit(filter_by_weight(run_directory))
                output_str += f"Uploaded filter_by_weight result vector 'keep68' to {run_directory}/cluster.txt\n"
        
        '''
        if is_checked(json_variables, "cluster_by_rmsd_checkbox"):
                structure_directory = json_variables["structure_directory"][0]
                cluster_by_rmsd(run_directory, structure_directory)
                output_str += f"Uploaded cluster_by_rmsd result matrices 'weights' and 'structs' to {run_directory}/cluster.txt\n"
        '''

        #        ga.udpmessage( { "_textarea" : "udp message to _textarea\n" } );

#        ga.udpprogress(0.5);

#        ga.tcpmessage( { "_textarea" : "tcp message to _textarea\n" } );
#        ga.tcpmessage( { "_textarea" : "JSON input to executable:\n" + json.dumps( json_variables, indent=4 ) + "\n" } );

        # output final object to the UI
        output['progress_html'] = 1.0
        output['_textarea'] = output_str

        if is_checked(json_variables, "cluster_by_rmsd_checkbox"):
                tree_filename = json_variables["tree_input"][0]
                rmsd_cut = json_variables["rmsd_cut"]

                tree = np.loadtxt(tree_filename)

                keep68 = np.loadtxt(os.path.join(run_directory, "cluster.txt"))[:,0].astype(int)
                #printQuit(keep68)
                index68 = np.loadtxt(os.path.join(run_directory, "index.txt")).astype(int)

                xsol = np.loadtxt(os.path.join(run_directory, "weights_for_all_lambdas.txt"))[:, index68]

                T = fcluster(Z=tree, t=rmsd_cut, criterion="distance") # the RMSD clustering is chose to be 4A in here

                size = np.size(keep68)
                new_T = np.zeros((size, 2))

                new_T[:,0] = keep68        
                new_T[:,1] = T[keep68-1]

                #T = T_matrix[T_matrix[:,0].argsort()]
                new_T = sortrows(new_T) # using lambda function


                # Here, the representative of the cluster is chosen
                # The representative is the structure with the highest initial weight (assigned by max entropy)
                # The final weight of the representative is the summation of weights for every point in the cluster

                aux1 = np.unique(new_T[:,1])
                num_clusters = np.size(aux1)

                weights = np.empty(num_clusters)
                structs = np.empty(num_clusters)
                
                for i in range (num_clusters):
                        aux = np.where(new_T[:,1] == aux1[i])[0]
                        x = xsol[new_T[aux,0].astype(int)-1]
                        weights[i] = np.sum(x)
                        structs[i] = np.where(xsol == np.max(x))[0]+1

                struct_sort = np.argsort(structs)
                structs = structs[struct_sort]
                weights = weights[struct_sort]
                weights = np.divide(weights, np.sum(weights)) # Normalize weights
                
                np.savetxt(os.path.join(run_directory, "result_after_cluster.txt"), np.column_stack((structs, weights)))
                output_str += f"Uploaded structure vectors 'structs' and 'weights' to {run_directory}/result_after_cluster.txt\n"

        print(json.dumps(output))

'''
{'directory_name': 'run_0', 'maxent_checkbox': 'on', '_logon': 'mcaserta', '_project': '', '_window': 'afbb1a6b-f1de-4c90-8180-e80de61f9b3f', '_uuid': 'ea502d30-f9b2-11eb-92b3-1f3dd8202a3b', '_html_maxent_checkbox-datafile_altval': '<i>Server</i>: 1_3_21_Test/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt', '_html_maxent_checkbox-matrixfile_altval': '<i>Server</i>: 1_3_21_Test/Amatrix_RDC_prox.txt', '_base_directory': 'results/users/mcaserta/no_project_specified', '_log_directory': 'results/users/mcaserta/no_project_specified/_log', '_tcphost': '10.0.0.8', '_tcpport': 30780, '_udphost': '10.0.0.8', '_udpport': 30779, 'resourcedefault': 'local', '_webroot': '/var/www/html', '_application': 'parnmr', '_menu': 'maxent_menu', '_module': 'maxent', 'lambda_lower_bound': '-0.6', 'lambda_step_size': '0.2', 'lambda_upper_bound': '7', 'datafile': ['results/users/mcaserta/1_3_21_Test/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt'], 'matrixfile': ['results/users/mcaserta/1_3_21_Test/Amatrix_RDC_prox.txt']}
'''
