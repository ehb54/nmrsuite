#!/opt/miniconda3/bin/python

import json, sys, time, os, re, shutil
from io import StringIO ## for Python3
from genapp3 import genapp ## python3
from predict_scripts.predictPCS import predictPCS
from predict_scripts.predictPRE import predictPRE
from predict_scripts.predictRDC import predictRDC
from predict_scripts.components2tensor import components2tensor

import numpy as np

# function to check if a checkbox was checked by seeing if its id exists in the JSON
is_checked = lambda json_variables, checkbox_id : checkbox_id in json_variables.keys()

line_csv_to_arr = lambda line_csv : np.array(line_csv.split(",")).astype(float)

def getTensorFromEuler (user_input):
    inputs = user_input.split(",")
    inputs = list(map(float, inputs))
    [Aax, Arh, a, b, g] = inputs
    susceptibility_tensor = components2tensor(Aax, Arh, a, b, g)
    return susceptibility_tensor

def inputToList (input_text):
    text = re.sub(r'\s*:\s*', ':', input_text)
    text = re.findall(r'\s|,|[^,\s]+', text)
    output_list = []
    for element in text:
        if (element.contains(":")):
            [range_start, range_end] = element.split(":")
            output_list.extend(list(range(int(range_start.strip()), int(range_end.strip())+1)))
        else:
            output_list.append(int(element.strip()))

    return output_list # Not numpy array
            

def getResidueNumbers (inclusion_text, exclusion_text):
    inclusion_set = set(inputToList(inclusion_text))
    exclusion_set = set(inputToList(exclusion_text))
    residue_numbers = np.array(list(inclusion_set - exclusion_set))
    
    return residue_numbers

if __name__=='__main__':

        argv_io_string = StringIO(sys.argv[1])
        json_variables = json.load(argv_io_string)

        output_str = ""
        output = {}
        run_directory = json_variables["runname"]

        ### initialize the genapp object
        ga = genapp( json_variables )

        pre_chosen = is_checked(json_variables, "choose_pre")
        pcs_chosen = is_checked(json_variables, "choose_pcs")
        rdc_chosen = is_checked(json_variables, "choose_rdc")

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
                shutil.move(filename, "PDB")#shutil.move(os.path.join(project_prefix, filename), "PDB")

            #command += f" -pdb {project_prefix + 'PDB'}"

        elif (choice == "c3"):
            server_path = json_variables["pdbserverpath"][0]

            new_server_path = os.getcwd()
            #printQuit(new_server_path)

            f1 = os.getcwd()
            f2 = server_path.split("users/")[1] # Takes all after users
            new_server_path = joinFilePaths(f1, f2)

            #command += f" -pdb {new_server_path}"

        else:
            raise Exception (f"There is no support for the choice '{choice}'")

            #command += f" -pdb {project_prefix + 'PDB'}"

        output_str = ""

        if (pre_chosen or pcs_chosen):
            SL_position = line_csv_to_arr(json_variables["paramagnetic_coord"])

        if (pre_chosen or rdc_chosen):
            freq = float(json_variables["frequency"])

        if (pre_chosen):
            T2dia = float(json_variables["T2_dia"])
            TAUc = float(json_variables["tauc"])
            Htime = float(json_variables["exptime_1H"])
            spin = float(json_variables["spin"])
            gammaRatio = 1 # Not sure what to do with this

            output_str += str(predictPRE (T2dia, TAUc, Htime, SL_position, freq, pdb_filename, spin, gammaRatio))

        if (pcs_chosen or rdc_chosen):
            pdb_filename = os.path.join(os.path.join(project_prefix, "PDB"), "mnodesS00002.pdb")#json_variables["pdbfile"][0]

            reslist = []

            if (json_variables["residue_flag"] == "c2"):
                inclusion_file = json_variables["inclusionfile"]
                exclusion_file = json_variables["exclusionfile"]
                inclusion_text = ", ".join(open("inclusion_file", "r").readlines())
                exclusion_text = ", ".join(open("inclusion_file", "r").readlines())
                reslist = getResidueNumbers(inclusion_text, exclusion_text) # Originally named residue_numbers

            if (json_variables["residue_flag"] == "c3"):
                inclusion_text = json_variables["inclusiontext"]
                exclusion_text = json_variables["exclusiontext"]
                reslist = getResidueNumbers(inclusion_text, exclusion_text) # Originally named residue_numbers
                
            #reslist = np.array([-0.2073, 0.2761, 0.1499, -0.1655, 0.3134]) #FIX
            pdb_model = 1 # Not sure what to do with this
            if (json_variables["listboxsuscept"] == "c1"):
                susceptibility_tensor = line_csv_to_arr(json_variables["susceptibility_c1"])
            else:
                susceptibility_tensor = getTensorFromEuler(json_variables["susceptibility_c2"])

        if (pcs_chosen):
            atom_type = json_variables["atomtypes"]
            output_str += str(predictPCS(susceptibility_tensor, reslist, pdb_filename, SL_position, atom_type, pdb_model))
        
        if (rdc_chosen):
            [atom1, atom2] = json_variables["atomtypes"].split(",")
            Sorder = float(json_variables["order_parameter"])
            
            output_str += str(predictRDC(susceptibility_tensor, reslist, pdb_filename, freq, atom1, atom2, Sorder, pdb_model))
        
        
        output['_textarea'] = output_str
        
        print(json.dumps(output))
