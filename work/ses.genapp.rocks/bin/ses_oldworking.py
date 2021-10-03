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

    fig.data = []

    fig.append_trace(go.Scatter(x=x, y=y), row=1, col=1)
    fig.append_trace(go.Scatter(x=x, y=y), row=1, col=2)

    fig.update_yaxes(title_text="Relative Error", row=1, col=1)
    fig.update_yaxes(title_text="Log Relative Error", type="log", row=1, col=2)

    fig.update_xaxes(title_text="Ensemble Size")

    fig.update_layout(height=600, width=1000, title_text="L-curve Plots")

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

    field_ids_to_exclude = ["runname"] # field ids that won't be used for the command
    runname = str(json_variables["runname"])

    """ output = str(data) + "\n\n\n"
    for key, value in data.items():
        output += key + "\n"

    return (output, runname) """

    # Gets role id and type for each field variable, puts into 3 lists
    field_roles = [str(element['role']) for element in data['fields']]
    field_ids = [str(element['id']) for element in data['fields']]
    field_types = [str(element['type']) for element in data['fields']]    

    indices_to_remove = []

    # Removes fields with the wrong id
    for i, field_id in enumerate(field_ids):
        if (field_ids[i] in field_ids_to_exclude):
            indices_to_remove.append(i)

    for index in sorted(indices_to_remove, reverse=True):
        del field_roles[index]
        del field_ids[index]
        del field_types[index]

    
    length = len(field_roles)
    if (length != len(field_ids) or length != len(field_types)):
        raise Exception ("ERROR => Invalid module definition JSON.")

    field_ids_to_field_types = OrderedDict()
    input_length = 0

    # Creates dictionary mapping field ids to field types only for fields that have the 'input' role
    # and a type that's not 'label' since the others are irrelevant to the program
    for i in range (length):
        if (field_roles[i] == "input" and field_types[i] != "label"):
            field_ids_to_field_types[field_ids[i]] = field_types[i]
            input_length += 1

    command = "java -jar {}sesgeneral-1.1.jar".format(bin_prefix)

    # Generation of "field_flags_to_field_ids" array
    field_flags_to_field_ids = {} # For each field, contains {[command line flag]: [JSON input ID]}
    field_flags_to_field_ids["K"] = "number_topsolution"
    field_flags_to_field_ids["align"] = "align"
    field_flags_to_field_ids["best"] = "nnls"
    field_flags_to_field_ids["data"] = "datafile"
    field_flags_to_field_ids["l0max"] = "l0max"
    field_flags_to_field_ids["matrix"] = "matrixfile"
    field_flags_to_field_ids["maxsum"] = "maxsum"
    field_flags_to_field_ids["out"] = "outputdir"
    field_flags_to_field_ids["outalign"] = "outalign"
    field_flags_to_field_ids["pdb"] = "pdbdirectory"
    field_flags_to_field_ids["precond"] = "precond"
    field_flags_to_field_ids["reltol"] = "reltol"
    field_flags_to_field_ids["rmsd"] = "rmsd"
    #field_flags_to_field_ids["storejava"] = "storejava"
    field_flags_to_field_ids["top"] = "toperror"

    # Reverses the dictionary
    field_ids_to_field_flags = {v: k for k, v in field_flags_to_field_ids.items()}

    field_ids_to_field_types_list = list(field_ids_to_field_types.items())

    json_variables_keys = json_variables.keys()

    # Generation of "field_flags_to_field_inputs" array
    field_flags_to_field_inputs = {} # For each field, contains {[command line flag]: [user input for field]}
    for i in range (input_length):
        field_id, field_type = field_ids_to_field_types_list[i]
        if field_id in json_variables_keys:
            field_raw_input = json_variables[field_id]
            if (field_id not in field_ids_to_field_flags.keys()):
                raise Exception ("ERROR => There is no output flag that matches to the field id '{}'. Please change the 'field_flags_to_field_ids dictionary.".format(field_id))
            field_flag = field_ids_to_field_flags[field_id]
            try:
                field_input = changeType(field_raw_input, field_type)
            except:
                raise Exception ("ERROR => Cannot convert {} to the type '{}'".format(field_raw_input, field_type))
            field_flags_to_field_inputs[field_flag] = field_input


    # Creates the command that needs to be run in the command line
    for field_flag, field_input in field_flags_to_field_inputs.items():
        if (field_flag not in field_flags_to_field_ids.keys()):
                raise Exception ("ERROR => There is no field id set to correspond with the field flag '{}'. Please change the 'field_flags_to_field_ids dictionary.".format(field_flag))
        else:
            field_id = field_flags_to_field_ids[field_flag]
        field_type = field_ids_to_field_types[field_id]
        if field_type == 'checkbox':
            command += ' -{}'.format(str(field_flag))
        else:
            command += ' -{} "{}"'.format(str(field_flag), str(field_input))

    # Creates output directory if needed
    output_directory = field_flags_to_field_inputs["out"]
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    printQuit(command)

    # Runs the command and captures the output
    cmd_output = str(subprocess.run(str(command) + "2>&1", stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8'))
    lines = cmd_output.split("\n")
    #try:
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
    half_cmd_output = cmd_output[:about_mid]

    # Creates and returns a string containing the JSON input and output
    output = "JSON input to executable:\n" + json.dumps(json_variables, indent=4) + "\n\n"
    output += "JSON output from executable:\n" + half_cmd_output + "\n"
    return output, Lcurve_y, plot, runname


if __name__ == '__main__':
    output_string, Lcurve, plot, runname = main()

    f = open("{}.txt".format(runname), "w")
    f.write(output_string)  
    f.close()

    f = open("Lcurve.txt", "w")
    f.write("L-curve: {}".format(str(Lcurve)))
    f.close()

    # Saves L0 solutions in their respective files
    fig = saveAndPlot(output_string.split("\n"))

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

    # Archives file if it already contains stuff (only archives one iteration)
    if (os.path.exists("{}.txt".format(runname))):
        archive_file = open("{}_archive.txt".format(runname), "w")
        current_file_text = open("{}.txt".format(runname), "r").read()
        archive_file.write(current_file_text)
        archive_file.close()



