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
    print (text)
    for element in text:
        if (element.strip() == "" or element.strip() == ","): # Disregard if blank
            continue
        elif (":" in element):
            [range_start, range_end] = element.split(":")
            output_list.extend(list(range(int(range_start.strip()), int(range_end.strip())+1)))
        else:
            output_list.append(int(element.strip()))

    return output_list # Not numpy array

def printQuit (string): # Prints a desired string
    output = {}
    output["_textarea"] = str(string)
    print (json.dumps(output))
    quit()

def joinFilePaths (f1, f2):
    return f1 + '/' + '/'.join([i for i in f2.split('/') if i not in f1.split('/')])


def getResidueNumbers (inclusion_text, exclusion_text):
    inclusion_set = set(inputToList(inclusion_text))
    exclusion_set = set(inputToList(exclusion_text))
    residue_numbers = np.array(list(inclusion_set - exclusion_set))

    return residue_numbers

def makeBarChart (predicted_data, data_type):
    residue_numbers = predicted_data[:,0].tolist()
    chart = {
            "data": [
                {
                    "x": residue_numbers, # Convert to string
                    "y": predicted_data[:,1].tolist(),
                    "text": residue_numbers,
                    "type": "bar",
                    "marker": {
                        "color": "Red"
                    }
                }
            ],
            "layout": {
                    "title": f"Predicted {data_type} Data for the First Structure",
                    "xaxis": {
                        "title": "Residue Number"
                    },
                    "yaxis": {
                        "title": f"Predicted {data_type} Value",
                        #"range": [np.min(leg[:,0])-0.05*np.min(leg[:,0]), np.max(leg[:,0])+0.1*np.max(leg[:,0])]
                    },
            }
        }
    
    return chart

if __name__=='__main__':

        argv_io_string = StringIO(sys.argv[1])
        json_variables = json.load(argv_io_string)

        output_str = ""
        output = {}
        runname = json_variables["runname"]

        ### initialize the genapp object
        ga = genapp( json_variables )

        pre_chosen = is_checked(json_variables, "choose_pre")
        pcs_chosen = is_checked(json_variables, "choose_pcs")
        rdc_chosen = is_checked(json_variables, "choose_rdc")

        if not os.path.exists(runname):
            os.mkdir(runname)
        os.chdir(runname)

        if not os.path.exists("Predict"):
            os.mkdir("Predict")
        os.chdir("Predict")

        if os.path.exists("PDB"):
            shutil.rmtree("PDB")
        os.mkdir("PDB")

        os.chdir("../")

        choice = str(json_variables["pdbinput"]) # choice will be c1, c2, or c3

        if (choice == "c1" or choice == "c2"):
            if (choice == "c1"):
                files = json_variables["pdblocaldirectory"] # C1 returns a list of files

            if (choice == "c2"):
                files = json_variables["pdblocalfiles"] # C2 returns a list of files

            for filename in files:
                shutil.move(filename, os.path.join("Predict", "PDB"))#shutil.move(os.path.join(project_prefix, filename), "PDB")

        elif (choice == "c3"):
            server_path = json_variables["pdbserverpath"][0]

            new_server_path = os.getcwd()
            #printQuit(new_server_path)

            f1 = os.getcwd()
            f2 = server_path.split("users/")[1] # Takes all after users
            new_server_path = joinFilePaths(f1, f2)

        else:
            raise Exception (f"There is no support for the choice '{choice}'")


        os.chdir("Predict")

        pdb_filenames = []
        for file in files:
            filename = file.rsplit("/", 1)[-1]
            new_filename = "PDB/" + filename
            pdb_filenames.append(new_filename)

        pdb_model = int(json_variables["pdb_model"]) - 1
        chainID = str(json_variables["chainID"]).strip()

        reslist = np.array([])

        if (json_variables["residue_flag"] == "c2"):
            inclusion_file = json_variables["inclusionfile"]
            exclusion_file = json_variables["exclusionfile"]
            inclusion_text = ", ".join(open(inclusion_file[0], "r").readlines())
            exclusion_text = ", ".join(open(exclusion_file[0], "r").readlines())
            reslist = getResidueNumbers(inclusion_text, exclusion_text) # Originally named residue_numbers

        if (json_variables["residue_flag"] == "c3"):
            inclusion_text = json_variables["inclusiontext"]
            exclusion_text = json_variables["exclusiontext"]
            reslist = getResidueNumbers(inclusion_text, exclusion_text) # Originally named residue_numbers


        #printQuit(pdb_filenames)

        '''
        pdb_filename = os.path.join("PDB", "mnodesS00001.pdb")#json_variables["pdbfile"][0]
        pdb_number = re.sub(r'[^0-9]', '', pdb_filename) # Finds number in string (e.g. "mnodesS00001" is "00001")
        pdb_number_count = len(re.sub(r'[^0-9]', '', pdb_filename)) # Counts length of number)

        num_files = 100

        #string = f"Current Working Directory: {os.getcwd()}\nResult of 'ls': {os.listdir()}"
        #printQuit(string)

        for num in range (2, num_files + 1):
            pdb_filename_current = pdb_filename.replace(pdb_number, str(num).zfill(pdb_number_count))
            pdb_filenames.append(pdb_filename_current)
        '''

        if (pre_chosen or pcs_chosen):
            SL_position = line_csv_to_arr(json_variables["paramagnetic_coord"])
            atom_type = json_variables["atomtypes"]

            if ("," in atom_type or " " in atom_type): #If more than one atom type was given, will only pick the first
                atom_type = atom_type.split(",")[0].strip()

        if (pre_chosen or rdc_chosen):
            freq = float(json_variables["frequency"])

        if (pre_chosen):
            T2dia = float(json_variables["T2_dia"])
            TAUc = float(json_variables["tauc"])
            Htime = float(json_variables["exptime_1H"])
            spin = float(json_variables["spin"])
            gammaRatio = 1 # Not sure what to do with this

            predictedPREIsEmpty = True

            predictedPRE = ""
            for pdb_filename in pdb_filenames:
                current_predictedPRE = predictPRE(reslist, T2dia, TAUc, Htime, SL_position, freq, pdb_filename, spin, gammaRatio, pdb_model, chainID, atom_type)
                if (predictedPREIsEmpty):
                    predictedPRE = current_predictedPRE
                    #printQuit(predictedPRE)
                    predictedPREIsEmpty = False
                else:
                    predictedPRE = np.concatenate((predictedPRE, np.reshape(current_predictedPRE[:, -1], (-1, 1))), axis=1)

            np.savetxt("predictedPRE.txt", predictedPRE)
            output_str += "PRE:\n" + str(predictedPRE) + "\n"

        if (pcs_chosen or rdc_chosen):
            #reslist = np.array([-0.2073, 0.2761, 0.1499, -0.1655, 0.3134]) #FIX
            if (json_variables["listboxsuscept"] == "c1"):
                susceptibility_tensor = line_csv_to_arr(json_variables["susceptibility_c1"])
            else:
                susceptibility_tensor = getTensorFromEuler(json_variables["susceptibility_c2"])


        if (pcs_chosen):

            predictedPCS = ""
            predictedPCSIsEmpty = True

            for pdb_filename in pdb_filenames:
                current_predictedPCS = predictPCS(susceptibility_tensor, reslist, pdb_filename, SL_position, atom_type, pdb_model, chainID)
                #printQuit(current_predictedPCS)
                if (predictedPCSIsEmpty):
                    predictedPCS = current_predictedPCS
                    predictedPCSIsEmpty = False
                    #printQuit(predictedPCS)
                else:
                    predictedPCS = np.concatenate((predictedPCS, np.reshape(current_predictedPCS[:, -1], (-1, 1))), axis=1)
            #printQuit(predictedPCS)
            np.savetxt("predictedPCS.txt", predictedPCS)
            output_str += "PCS:\n" + str(predictedPCS) + "\n"

        if (rdc_chosen):
            [atom1, atom2] = json_variables["atomtypes"].split(",")
            atom1 = atom1.strip()
            atom2 = atom2.strip()
            Sorder = float(json_variables["order_parameter"])
            T = float(json_variables["temp"])

            predictedRDC = ""
            predictedRDCIsEmpty = True

            #printQuit(str(atom1) + "\n" + str(atom2))

            for pdb_filename in pdb_filenames:
                current_predictedRDC = predictRDC(susceptibility_tensor, reslist, pdb_filename, freq, atom1, atom2, Sorder, pdb_model, chainID, T)
                if (predictedRDCIsEmpty):
                    predictedRDC = current_predictedRDC
                    predictedRDCIsEmpty = False
                else:
                    #string = f"First: {np.shape(predictedRDC)} \nSecond: {np.shape(np.reshape(current_predictedRDC[:, 1], (-1, 1)))} \nThird: {np.shape(current_predictedRDC[:, 1])}"
                    #printQuit(string)
                    predictedRDC = np.concatenate((predictedRDC, np.reshape(current_predictedRDC[:, 1], (-1, 1))), axis=1)

            np.savetxt("predictedRDC.txt", predictedRDC)
            output_str += "RDC:\n" + str(predictedRDC) + "\n"

        output['_textarea'] = output_str

        if (pre_chosen): output['pre_graph'] = makeBarChart(predictedPRE, "PRE")
        if (pcs_chosen): output['pcs_graph'] = makeBarChart(predictedPCS, "PCS")
        if (rdc_chosen): output['rdc_graph'] = makeBarChart(predictedRDC, "RDC")

        print(json.dumps(output))
