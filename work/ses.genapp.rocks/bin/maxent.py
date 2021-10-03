#!/opt/miniconda3/bin/python

import json, sys, time, os
from io import StringIO ## for Python3
from genapp3 import genapp ## python3
from maxent_scripts.RunMaxEntropy import run_max_entropy
from maxent_scripts.L_curveForBestLambda import L_curve_for_best_lambda
from maxent_scripts.FilterByWeight import filter_by_weight
from maxent_scripts.ClusterByRMSD import cluster_by_rmsd


# function to check if a checkbox was checked by seeing if its id exists in the JSON
is_checked = lambda json_variables, checkbox_id : checkbox_id in json_variables.keys()


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
                run_max_entropy (A_filename, y_filename, lambda_lower, lambda_step, lambda_upper, run_directory)
                line_plot, histogram = L_curve_for_best_lambda(run_directory)
                output_str += f"Uploaded maximum entropy result files (A.txt, index.txt, lambda.txt, x.txt, y.txt) to {run_directory}\n"
                output['lineplot'] = line_plot
                output['histogram'] = histogram
                output_str += f"Created L_curve plots\n"


        if is_checked(json_variables, "filter_by_weight_checkbox"):
                filter_by_weight(run_directory)
                output_str += f"Uploaded filter_by_weight result vector 'keep68' to {run_directory}/cluster.txt\n"
        
        if is_checked(json_variables, "cluster_by_rmsd_checkbox"):
                structure_directory = json_variables["structure_directory"][0]
                cluster_by_rmsd(run_directory, structure_directory)
                output_str += f"Uploaded cluster_by_rmsd result matrices 'weights' and 'structs' to {run_directory}/cluster.txt\n"

        ga.udpprogress(0);
        ga.udpmessage( { "_textarea" : "udp message to _textarea\n" } );

        ga.udpprogress(0.5);

        ga.udpmessage( { "_textarea" : "tcp message to _textarea\n" } );
        #output += str(ga) + "\n"

        # output final object to the UI
        output['progress_html'] = 1.0
        output['_textarea'] = output_str
        print(json.dumps(output))

'''
{'directory_name': 'run_0', 'maxent_checkbox': 'on', '_logon': 'mcaserta', '_project': '', '_window': 'afbb1a6b-f1de-4c90-8180-e80de61f9b3f', '_uuid': 'ea502d30-f9b2-11eb-92b3-1f3dd8202a3b', '_html_maxent_checkbox-datafile_altval': '<i>Server</i>: 1_3_21_Test/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt', '_html_maxent_checkbox-matrixfile_altval': '<i>Server</i>: 1_3_21_Test/Amatrix_RDC_prox.txt', '_base_directory': 'results/users/mcaserta/no_project_specified', '_log_directory': 'results/users/mcaserta/no_project_specified/_log', '_tcphost': '10.0.0.8', '_tcpport': 30780, '_udphost': '10.0.0.8', '_udpport': 30779, 'resourcedefault': 'local', '_webroot': '/var/www/html', '_application': 'parnmr', '_menu': 'maxent_menu', '_module': 'maxent', 'lambda_lower_bound': '-0.6', 'lambda_step_size': '0.2', 'lambda_upper_bound': '7', 'datafile': ['results/users/mcaserta/1_3_21_Test/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt'], 'matrixfile': ['results/users/mcaserta/1_3_21_Test/Amatrix_RDC_prox.txt']}
'''

