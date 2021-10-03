#!/usr/bin/python3.6

import json
import sys
import os
import subprocess
import re
import shutil
import numpy as np
import pandas as pd
import plotly
from plotly.subplots import make_subplots
from io import StringIO                                      # <--- import Python's JSON Library
from collections import OrderedDict                          # <--- import Ordered Dictionary
import plotly.graph_objects as go

from saveL0solutions import saveAndPlot
#from parnmr.bin.saveL0solutions import saveL0Solutions

# field types = 'text', 'lrfile', 'integer', 'checkbox', 'label', 'listbox', 'float'

def blockPrint(): # prevents python module from printing
    sys.stdout = open(os.devnull, 'w')

def enablePrint(): # restores python module's printing ability
    sys.stdout = sys.__stdout__

def changeType (input_value, new_type): # function to change type of variable
    if (new_type in ['text', 'listbox']):
        return str(input_value)
    elif (new_type == 'lrfile'):
        #head, pdbfilename = os.path.split(input_value[0])
        return str(input_value[0])
    elif (new_type == 'integer'):
        return int(input_value)
    elif (new_type == 'checkbox'):
        return bool(input_value)
    elif (new_type == 'float'):
        return float(input_value)
    else:
        raise Exception ("Invalid field type '{}'.".format(new_type))

def printQuit (string):
    output = {}
    output["_textarea"] = str(string)
    print (json.dumps(output))
    quit()

def getPlot (x, y, log = False):
    fig = make_subplots(rows=1, cols=2)

    fig.append_trace(go.Scatter(x=x, y=y), row=1, col=1)
    fig.append_trace(go.Scatter(x=x, y=y), row=1, col=2)

    fig.update_yaxes(title_text="Relative Error", row=1, col=1)
    fig.update_yaxes(title_text="Log Relative Error", type="log", row=1, col=2)

    fig.update_xaxes(title_text="Ensemble Size")

    fig.update_layout(height=600, width=1000, title_text="L-curve Plots")
    fig.update_traces(marker=dict(color='Blue'))

    fig.data = [fig.data[-2], fig.data[-1]]

    return (fig.to_dict())


    # convert to dict
    #raise Exception ("Chart printed below \n\n {}. Chart variables are below \n\n".format(json.dumps(chart, indent=4)))
    


def main():
    argv_io_string = StringIO(sys.argv[1])                   # <--- Read in input JSON
    json_variables = json.load(argv_io_string)               # <--- Put input JSON into array
    
    #java -jar sesgeneral-1.1.jar -out "output" -matrix "A_data.txt" -data "y_data.txt"
    #json_variables = {"outputdir": "output", "matrixfile": "A_data.txt", "datafile": "y_data.txt"}

    from os.path import dirname, abspath

    bin_prefix = dirname(abspath(__file__)) + "/"

    module_prefix = dirname(dirname(abspath(__file__))) + "/" + "modules" + "/"

    raw_json_file = module_prefix + "ses.json"
    """ adjusted_json_file = bin_prefix + "ses_fixed.json"

    #adjusted_json = open(adjusted_json_file, "w")
    adjusted_json = ""
    for line in open(raw_json_file, "r").readlines():
        if (not line.startswith("#")):
            #adjusted_json.write(line)
            adjusted_json += line

    data = json.loads(adjusted_json) """

    with open(raw_json_file) as f:
        data = json.load(f)

    """ adjusted_json = open(adjusted_json_file, "r")
    data = json.load(adjusted_json)
    adjusted_json.close() """

    
    """

    json_file = module_prefix + "ses.json"
    json = open(json_file, "r").readlines()
    json_string = ""
    new_json = [line for line in json if not line.strip().startswith("#")]
    json_string = ""
    for line in new_json

    output = str(data)


    for key, value in data.items():
        output += "\n" + (key)

    return (output)

    with open("ses.json") as f:
        data = json.load(f) """

    #printQuit(json_variables)
    command = "java -jar sesgeneral-1.1.jar -out"

    if ("module_header" in json_variables):
        module_header = str(json_variables["module_header"])

    if ("runname" in json_variables):
        runname = str(json_variables["runname"])
            
    if ("datafile" in json_variables):
        datafile = str(json_variables["datafile"][0])
        command += f" -data {datafile}"
            
    if ("matrixfile" in json_variables):
        matrixfile = str(json_variables["matrixfile"][0])
        command += f" -matrix {matrixfile}"
            
    if ("outputdir" in json_variables):
        outputdir = str(json_variables["outputdir"])
        command += f" -out {outputdir}"
            
    if ("number_topsolution" in json_variables):
        number_topsolution = str(json_variables["number_topsolution"])
        command += f" -K {number_topsolution}"
            
    if ("nnls" in json_variables):
        nnls = str(json_variables["nnls"])
        command += f" -best {nnls}"
            
    if ("advanced_input_label" in json_variables):
        advanced_input_label = str(json_variables["advanced_input_label"])
        command += f" -best {nnls}"
            
    if ("pdb_flag" in json_variables):
        pdb_flag = str(json_variables["pdb_flag"])
    
    if ("pdbdirectory" in json_variables):
        pdbdirectory = str(json_variables["pdbdirectory"])
        command += f" -pdb {pdbdirectory}"
        
    if ("align" in json_variables):
        align = str(json_variables["align"])
        command += f" -align {align}"
        
    if ("outalign" in json_variables):
        outalign = str(json_variables["outalign"])
        command += f" -outalign {outalign}"
        
    if ("rmsd" in json_variables):
        rmsd = str(json_variables["rmsd"])
        command += f" -rmsd {rmsd}"
        
    if ("l0max" in json_variables):
        l0max = str(json_variables["l0max"])
        command += f" -l0max {l0max}"
        
    if ("maxsum" in json_variables):
        maxsum = str(json_variables["maxsum"])
        command += f" -maxsum {maxsum}"
        
    if ("precond" in json_variables):
        precond = str(int(json_variables["precond"])-1)
        command += f" -precond {precond}"
        
    if ("reltol" in json_variables):
        reltol = str(json_variables["reltol"])
        command += f" -reltol {reltol}"
        
    if ("toperror" in json_variables):
        toperror = str(json_variables["toperror"])
        command += f" -top {toperror}"

    if not os.path.exists(runname):
        os.mkdir(runname)

    os.chdir(runname)

    if not os.path.exists("SES"):
        os.mkdir("SES")

    os.chdir("SES")
    

    #command = f"""-K {number_topsolution} -align {align} -best {nnls} -data {datafile} -l0max {l0max}
    #-matrix {matrixfile} -maxsum {maxsum} -out {outputdir} -outalign {outalign} -pdb {pdbdirectory}
    #-precond {precond} -reltol {reltol} -rmsd {rmsd} -top {toperror}"""

    #printQuit(command)

    '''
    java -jar sesgeneral-1.1.jar -out -data results/users/mcasertano/no_project_specified/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt
    -matrix results/users/mcasertano/no_project_specified/Amatrix_RDC_prox.txt -out solution -K 100000 -l0max 2147483647
    -maxsum 99999999999 -precond 0 -reltol 0.0005 -top 0.005
    '''
    '''
    java -jar /opt/genapp/parnmr/bin/sesgeneral-1.1.jar -data "results/users/mcasertano/no_project_specified/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt"
    -matrix "results/users/mcasertano/no_project_specified/Amatrix_RDC_prox.txt" -out "solution" -K "100000" -l0max "2147483647" -maxsum "99999999999.0"
    -precond "1" -reltol "0.0005" -top "0.005"
    '''
    # Runs the command and captures the output
    cmd_output = str(subprocess.run(str(command) + "2>&1", stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8'))
    lines = cmd_output.split("\n")
    #printQuit(cmd_output)
    y = [x for x in lines if x.startswith("L-curve")][0]
    y = y[10:-1].split(", ")
    Lcurve_y = list(map(float, y))
    if (float(Lcurve_y[0] - 1.0) > 0.001):
        Lcurve_y.insert(0, 1.0)
    Lcurve_x = list(range(0, len(Lcurve_y) + 1))
    plot = getPlot (Lcurve_x, Lcurve_y)

    """ except:
        pass """
    """ try:
        cmd_output = subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8')
    except Exception as e:
        return ("ERROR => There was an issue in running the below command.\n\t {} \n\nThe error is printed below: \n{}\n\nThis command was based on your input. Please double check that the JSON file is correct.".format(command, e)) """

    # Separates the output into two (because it is duplicated for some reason)
    mid = int(len(cmd_output)/2)
    about_mid = mid + cmd_output[mid:].index('\n\n')
    half_cmd_output = cmd_output#[:about_mid]

    # Creates and returns a string containing the JSON input and output
    output = "JSON input to executable:\n" + json.dumps(json_variables, indent=4) + "\n\n"
    output += "JSON output from executable:\n" + half_cmd_output + "\n"
    return output, Lcurve_y, plot, runname


if __name__ == '__main__':
    output_string, Lcurve, plot, runname = main()

    f = open("{}.txt".format("output"), "w")
    f.write(output_string)
    f.close()

    f = open("Lcurve.txt", "w")
    f.write("L-curve: {}".format(str(Lcurve)))
    f.close()

    # Saves L0 solutions in their respective files
    fig = saveAndPlot(output_string.split("\n"))

    #printQuit(fig)

    # Wraps string into the text area box and returns it
    output = {}
    output['_textarea'] = output_string
    #output['_textarea'] += "\n\n" + str(plot) + "\n\n" + str(fig)
    output['lineplot'] = plot
    output['scatterplot'] = fig
    print (json.dumps(output))

    # Places standard output into runname folder
    if not os.path.exists("solution"):
        os.mkdir("solution")
    
    os.chdir("solution")
    os.rename("general_solutions.dat", "general_solutions.txt")
    """ # Archives file if it already contains stuff (only archives one iteration)
    if (os.path.exists("{}.txt".format(runname))):
        archive_file = open("{}_archive.txt".format(runname), "w")
        current_file_text = open("{}.txt".format(runname), "r").read()
        archive_file.write(current_file_text)
        archive_file.close()
 """


