#!/opt/miniconda3/bin/python

import json
import sys
import os
import subprocess
import re
import shutil
import time
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

# Below are three file manipulation methods
# Only the second is used currently, but the others are kept in case the second ever has issues
def absoluteFilePaths(directory):
    path_list = []
    for dirpath, _, filenames in os.walk(directory):
        for f in filenames:
            path_list.append(os.path.abspath(os.path.join(dirpath, f)))

    return path_list

def joinFilePaths (f1, f2):
    return f1 + '/' + '/'.join([i for i in f2.split('/') if i not in f1.split('/')])

def combine_with_duplicate(root, rel_path):
    rs = root.split("/")
    rps = rel_path.split("/")
    popped = False
    for v in rs:
        if v == rps[0]:
            rps.pop(0)
            popped = True
        elif popped:
            break

    return "/".join(rs+rps)

def stringToList (string): # This is used to convert the column and weights strings to a list
    return string[1:-2].split(" ")

# DO NOT USE THIS FUNCTION FOR FILEPATHS
# It often just leaves out huge sections of the filepath which caused a LOT of confusion in debugging
def printQuit (string): # Prints a desired string
    output = {}
    output["_textarea"] = str(string)
    print (json.dumps(output))
    quit()

def getPlot (x, y, log = False): # This is for the top relative error plot
    fig = make_subplots(rows=1, cols=2)

    fig.append_trace(go.Scatter(x=x, y=y, marker={"color": 'Blue'}, name="Linear"), row=1, col=1)
    fig.append_trace(go.Scatter(x=x, y=y, marker={"color": 'Red'}, name="Log"), row=1, col=2)

    fig.update_yaxes(title_text="Relative Error", row=1, col=1)
    fig.update_yaxes(title_text="Log Relative Error", type="log", row=1, col=2)

    fig.update_xaxes(title_text="Ensemble Size")

    fig.update_layout(height=600, width=1000, title_text="L-curve Plots")
    #fig.update_traces(marker=dict(color='Blue'))

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

    """ with open(raw_json_file) as f:
        data = json.load(f) """

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


    command = "java -jar {}sesgeneral-1.1.jar".format(bin_prefix)

    if ("module_header" in json_variables):
        module_header = str(json_variables["module_header"])

    if ("runname" in json_variables):
        runname = str(json_variables["runname"])

    if ("datafile" in json_variables):
        datafile = str(json_variables["datafile"][0])
        #command += f' -data "{datafile}"'
        command += f' -data {datafile}'

    if ("matrixfile" in json_variables):
        matrixfile = str(json_variables["matrixfile"][0])
        #command += f' -matrix "{matrixfile}"'
        command += f' -matrix {matrixfile}'

    if ("outputdir" in json_variables):
        outputdir = str(json_variables["outputdir"])
        #command += f' -out "{outputdir}"'
        command += f' -out {outputdir}'

    if ("number_topsolution" in json_variables):
        number_topsolution = str(json_variables["number_topsolution"])
        command += f" -K {number_topsolution}"

    if ("nnls" in json_variables):
        nnls = str(json_variables["nnls"])
        command += " -best"

    """ if ("advanced_input_label" in json_variables):
        advanced_input_label = str(json_variables["advanced_input_label"])
        command += f" -storejava {advanced_input_label}" """

        # PDB stuff
    if ("pdb_flag" in json_variables):
        project_prefix = os.getcwd() + "/"

        #pdb_flag = str(json_variables["pdb_flag"])
        choice = str(json_variables["pdbinput"]) # choice will be c1, c2, or c3

        if (choice == "c1" or choice == "c2"):
            if (os.path.exists("PDB")):
                shutil.rmtree("PDB")
            os.mkdir("PDB")
            if (choice == "c1"):
                files = json_variables["pdblocaldirectory"] # C1 returns a list of files
            if (choice == "c2"):
                files = json_variables["pdblocalfiles"] # C2 returns a list of files

            for filename in files:
                shutil.move(os.path.join(project_prefix, filename), "PDB")

            command += f" -pdb {project_prefix + 'PDB'}"


        elif (choice == "c3"):
            server_path = json_variables["pdbserverpath"][0]

            new_server_path = os.getcwd()
            #printQuit(new_server_path)

            f1 = os.getcwd()
            f2 = server_path.split("users/")[1] # Takes all after users
            new_server_path = joinFilePaths(f1, f2)

            command += f" -pdb {new_server_path}"

        else:
            raise Exception (f"There is no support for the choice '{choice}'")


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

    if ("developer_debug_flag" in json_variables):
        developer_debug_flag = str(json_variables["developer_debug_flag"])
    else:
        developer_debug_flag = ""


    #printQuit(f"Value: {str(developer_debug_flag)} Type: {type(developer_debug_flag)}")



    #command = f"""-K {number_topsolution} -align {align} -best {nnls} -data {datafile} -l0max {l0max}
    #-matrix {matrixfile} -maxsum {maxsum} -out {outputdir} -outalign {outalign} -pdb {pdbdirectory}
    #-precond {precond} -reltol {reltol} -rmsd {rmsd} -top {toperror}

    #java -jar sesgeneral-1.1.jar -out -data results/users/mcasertano/no_project_specified/yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt
    # -matrix results/users/mcasertano/no_project_specified/Amatrix_RDC_prox.txt -out solution -K 100000 -l0max 2147483647
    # -maxsum 99999999999 -precond 0 -reltol 0.0005 -top 0.005


    #printQuit(command)
    if not os.path.exists(runname):
        os.mkdir(runname)
    os.chdir(runname)

    if not os.path.exists("SES"):
        os.mkdir("SES")
    os.chdir("SES")

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    # Runs the command and captures the output
    cmd_output = str(subprocess.run(str(command) + " 2>&1", stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8'))
    output = f"Command: {command}      Output: {cmd_output}"
    lines = cmd_output.split("\n")
    #printQuit(cmd_output)
    try:
        y = [x for x in lines if x.startswith("L-curve")][0]
    except Exception as e:
        outputStr = ""
        outputStr += "Error with the command output. Please see output below: \n"
        outputStr += cmd_output
        printQuit(outputStr)

    cmd_output.split("===Best L0")[0]


    matrix = np.loadtxt(matrixfile)

    data = np.loadtxt(datafile)
    experimental_values = data[:, 0]
    experimental_errors = data[:, 1]

    # AT THIS POINT IN THE CODE, WE HAVE THE COMMAND AND THE OUTPUT
    # Go to the end of this function for an explanation of why we are generating "newOutput"

    # WE ARE GOING TO CREATE THE LO FOLDER
    L0_folder_name = "L0_Solutions"
    if os.path.exists(L0_folder_name):
        shutil.rmtree(L0_folder_name)

    time.sleep(.1) # Sometimes a delay is needed for the folder to completely delete
    os.mkdir(L0_folder_name)
    os.chdir(L0_folder_name)

    L0_solution_lines_start = lines.index("Best solutions by l0-norm (column index starts at 1):")
    L0_solution_lines_end = [lines.index(l) for l in lines if l.startswith("L-curve: ")][0]
    L0_solution_lines = list(range(L0_solution_lines_start + 1, L0_solution_lines_end - 1))

    newOutput = ""

    L0_numbers = []
    Chi2L_numbers = []
    rel_error_numbers = []

    for i, L0_solution_line in enumerate(L0_solution_lines):
        L0_solution_line = lines[L0_solution_line]
        try:
            m = re.match(r'l0-norm=(.*): Relative Error=(.*), Columns=(.*), Weights=(.*)', L0_solution_line)
            l0_norm = m.group(1)
            rel_error_value = m.group(2)
            columns = np.array(stringToList(m.group(3))).astype(int)-1
            weights = np.array(stringToList(m.group(4))).astype(float)
        except:
            printQuit (f"""L0_solution_line: {L0_solution_line}\n\n""")


        #output = f"{str(matrix)} \n {matrix.dtype} \n \n {str(columns)} \n {columns.dtype} \n \n {str(weights)} \n {weights.dtype} \n \n"
        #printQuit(output)

        #printQuit(f"{matrix} \n\n {columns} {matrix[:, columns]}")
        #predicted_values = np.multiply(matrix[:, columns], weights)

        try:
            if (columns.shape[0] == 1):
                predicted_values = np.multiply(matrix[:, columns], weights)
            else:
                predicted_values = np.matmul(matrix[:, columns], weights)
        except Exception as e:
            printQuit(f"""\n\nMatrix: {matrix}
            \n\nColumns: {columns}
            \n\nWeights: {weights}
            \n\nMatrix[:, columns] shape: {matrix[:, columns].shape}
            \n\nWeights shape: {weights.shape}
            \n\nException: {e}""")

        predicted_values = predicted_values.flatten()
        #printQuit(f"{experimental_errors.dtype}{np.square(np.subtract(experimental_values, predicted_values)).dtype}")
        #printQuit(np.square(np.subtract(experimental_values, predicted_values)))

        try:
            Chi2 = np.sum(np.square(np.divide(np.subtract(experimental_values, predicted_values), experimental_errors)))
        except:
            printQuit(f"""Experimental values:{experimental_values}
            \n\nExperimental errors: {experimental_errors}
            \n\nMatrix: {matrix}
            \n\nColumns: {columns}
            \n\nWeights: {weights}
            \n\nPredicted values: {predicted_values}
            \n\nMatrix[:, columns] shape: {matrix[:, columns].shape}
            \n\nWeights shape: {weights.shape}""")
        L = np.shape(data)[0] # Length

        indices = np.array(list(range(1, L+1)))
        rel_error_values = np.divide(np.subtract(predicted_values, experimental_values), experimental_errors)

        newOutput += f"===Best L0={l0_norm} Solution===\nRelative Error: {rel_error_value}\nChi2: {Chi2}\nL: {L}, Chi2/L: {Chi2/L}"

        #output = f"{tableString} \n\n Predicted Values: {str(predicted_values)} \n\n Experimental Values: {str(experimental_values)}"

        columns = ["<index>", "<exp. value>", "<pred. value>", "<relative err.>"]
        data = data[:, 0]

        # AT THIS POINT, ALL THE NEEDED LISTS HAVE BEEN GENERATED, SO WE CAN PRINT SOME OUTPUT IF DEBUGGING
        """
        printQuit(f"Data (shape: {data.shape}): {data}\n\nExperimental Values (shape: {experimental_values.shape}): {experimental_values}\n\nPredicted Values (shape: {predicted_values.shape}): {predicted_values}\n\nRelative Error Values (shape: {rel_error_values.shape}): {rel_error_values}\n\n")
        """

        data = np.column_stack((indices, experimental_values, predicted_values, rel_error_values))
        data = np.round(data, 3)

        template = "{0:10}|{1:15}|{2:15}|{3:15}"
        newOutput += "\n" + template.format(*columns) + "\n"
        for row in data:
            newOutput += template.format(*row) + "\n"

        newOutput += "\n\n"

        # SAVE IN FILE
        df = pd.DataFrame(data=data[:, 1:], index=data[:, 0], columns=columns[1:])
        if (l0_norm == 0):
            df.to_csv("BestPossibleSolution.txt", sep='\t')
        else:
            df.to_csv("L0={}.txt".format(l0_norm), sep='\t')


        #for i in range (L):
            #newOutput += "\n" + "%.5f %.5f %.5f %.5f" % (indices[i], experimental_values[i], predicted_values[i], rel_error_value[i])

        L0_numbers.append(int(l0_norm))
        Chi2L_numbers.append(float(Chi2/L))
        rel_error_numbers.append(float(rel_error_value))

    # The below code adds the best possible solution as an L0-norm=0 line to the extracted L0_solution_lines
    # This is because it can be worked with in the same way
    if ("nnls" in json_variables):
        columns = ["Index", "Exp. Data Info", "Exp. Value" "Pred. Value", "Relative Err."]

        table = []
        start_line = [lines.index(line) for line in lines if line.startswith("===Best x>0 Possible Solution===")][0]
        end_line = [lines.index(line) for line in lines if line.startswith("[")][0]
        relevant_lines = lines[start_line+5:end_line]

        rel_error_line = lines[start_line+1]
        rel_error_numbers.insert(0, float(re.match("Relative Error: (.*)", rel_error_line).group(1)))

        Chi2L_line = lines[start_line+3]
        Chi2L_numbers.insert(0, float(re.match("L: (.*), Chi2/L: (.*)", Chi2L_line).group(2)))

        L0_numbers.insert(0, 0)

        #printQuit(f"{L0_numbers}\n{rel_error_numbers}\n{Chi2L_numbers}")

        for line in relevant_lines:
            line = list(map(lambda x: x.strip(), line.split("\t")))
            table.append(line)

        table = np.array(table)

        new_columns = ["<index>", "<exp. value>", "<pred. value>", "<relative err.>"]

        df = pd.DataFrame(data=table[:, 2:], index=table[:, 0], columns=columns[1:])
        df.to_csv("BestPossibleSolution.txt", sep='\t')

    os.chdir("..")

    #printQuit(newOutput)
    #printQuit(f"{newOutput} \n {outputStr}")

    y = y[10:-1].split(", ")
    Lcurve_y = list(map(float, y))
    """ if (abs(float(Lcurve_y[0] - 1.0)) > 0.001):
        Lcurve_y.insert(0, 1.0) """
    Lcurve_x = list(range(1, len(Lcurve_y) + 2))
    plot = getPlot (Lcurve_x, Lcurve_y)

    """ except:
        pass """
    """ try:
        cmd_output = subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8')
    except Exception as e:
        return ("ERROR => There was an issue in running the below command.\n\t {} \n\nThe error is printed below: \n{}\n\nThis command was based on your input. Please double check that the JSON file is correct.".format(command, e)) """


    #newOutput = cmd_output.split("===Best L0")[0] + newOutput + "Generating the output" + cmd_output.split("Generating the output")[-1]

    # Below code basically replaces the middle part of the Java generated output with the "newOutput"
    # This is because the Java generated output stops listing solutions after a certain arbitrary L0
    # Therefore, newOutput contains the values for all solutions, but the beginning and middle part of the Java generated output
    # are still used because those don't have to change. They basically just re-print the output and other basic information.
    completeOutput = "JSON input to executable:\n" + json.dumps(json_variables, indent=4) + "\n\n"
    completeOutput += "JSON output from executable:\n" + cmd_output.split("===Best L0")[0] + newOutput
    completeOutput += "Generating the output" + cmd_output.split("Generating the output")[-1]
    return completeOutput, Lcurve_y, plot, runname, L0_numbers, Chi2L_numbers, rel_error_numbers


if __name__ == '__main__':
    completeOutput, Lcurve, plot, runname, L0_numbers, Chi2L_numbers, rel_error_numbers = main()

    f = open("{}.txt".format("output"), "w")
    f.write(completeOutput)
    f.close()

    f = open("Lcurve.txt", "w")
    f.write("L-curve: {}".format(str(Lcurve)))
    f.close()

    # Saves L0 solutions in their respective files
    fig = saveAndPlot(L0_numbers, Chi2L_numbers, rel_error_numbers)

    # Wraps string into the text area box and returns it
    output = {}
    output['_textarea'] = completeOutput #output_string
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




"""
java -jar sesgeneral-1.1.jar -data "yRDC_struct1_25_50_100_alt_noise_0p1_prox.txt" -matrix "Amatrix_RDC_prox.txt" -out "output" -K 100000 -best -l0max 2147483647 -maxsum 99999999999 -precond 0 -reltol 0.0005 -top 0.005
"""
