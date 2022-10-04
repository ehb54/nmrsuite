#!/opt/miniconda3/bin/python

from Bio import PDB
import numpy as np
import re, json, sys, os, shutil, itertools


from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

from pathlib import Path
from io import StringIO
from genapp3 import genapp ## python3

import numpy as np

# function to check if a checkbox was checked by seeing if its id exists in the JSON
is_checked = lambda json_variables, checkbox_id : checkbox_id in json_variables.keys()

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

def get_PDB_number (pdb_filename):
    import re
    number = int(re.findall(r'[0-9]+\.', pdb_filename)[0][:-1])
    return number
    
if __name__=='__main__':

        argv_io_string = StringIO(sys.argv[1])
        json_variables = json.load(argv_io_string)


        runname = json_variables["runname"]

        ### initialize the genapp object
        ga = genapp( json_variables )

        pre_chosen = is_checked(json_variables, "choose_pre")
        pcs_chosen = is_checked(json_variables, "choose_pcs")
        rdc_chosen = is_checked(json_variables, "choose_rdc")

        if not os.path.exists(runname):
            os.mkdir(runname)
        os.chdir(runname)
        
        if not os.path.exists("RMSD_Pairwise"):
            os.mkdir("RMSD_Pairwise")
        os.chdir("RMSD_Pairwise")
        
        if os.path.exists("PDB"):
            shutil.rmtree("PDB")
        os.mkdir("PDB")

        os.chdir("../")

        # PDB File Input

        #pdb_flag = str(json_variables["pdb_flag"])
        choice = str(json_variables["pdbinput"]) # choice will be c1, c2, or c3

        if (choice == "c1" or choice == "c2"):
            if (choice == "c1"):
                files = json_variables["pdblocaldirectory"] # C1 returns a list of files
    
            if (choice == "c2"):
                files = json_variables["pdblocalfiles"] # C2 returns a list of files

            #printQuit(str(files) + str(os.getcwd()) + str(os.listdir()))

            for filename in files:
                shutil.move(filename, os.path.join("RMSD_Pairwise", "PDB"))#shutil.move(os.path.join(project_prefix, filename), "PDB")
            
            #printQuit(files)
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
        
        os.chdir("RMSD_Pairwise")

        pdb_filenames = []
        for file in files:
            filename = file.rsplit("/", 1)[-1]
            new_filename = "PDB/" + filename
            pdb_filenames.append(new_filename)

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
        
    
        num_files = len(pdb_filenames)
        rmsd_array = np.zeros((num_files, num_files))

        reslistIsEmpty = (np.size(reslist) == 0)

        for (filename_a, filename_b) in list(itertools.combinations(pdb_filenames, 2)):
            parser = PDB.PDBParser(QUIET=True)
    
            struct_a = parser.get_structure('file', filename_a)
            struct_b = parser.get_structure('file', filename_b)

            model_a = struct_a[0]
            model_b = struct_b[0]

            atoms_a = np.array([])
            atoms_b = np.array([])

            for chain in model_a:
                for residue in chain:
                    if (reslistIsEmpty or np.isin(int(residue.id[1]), reslist)):
                        atoms_a = np.append(atoms_a, residue['CA'])#['CA'])
        
            for chain in model_b:
                for residue in chain:
                    if (reslistIsEmpty or np.isin(int(residue.id[1]), reslist)):
                        atoms_b = np.append(atoms_b, residue['CA'])#['CA'])


            super_imposer = PDB.Superimposer()
            #printQuit(str(np.size(atoms_a)) + "\n" + str(np.size(atoms_b)))
            super_imposer.set_atoms(atoms_a, atoms_b)
            super_imposer.apply(model_b.get_atoms())

            rmsd_array[get_PDB_number(filename_a)-1][get_PDB_number(filename_b)-1] = super_imposer.rms

        v = np.array([]) #v is going to be the vector that contains all the pairwise RMSD
        for i in range(num_files-1):
            v = np.append(v, rmsd_array[i, i+1:num_files])

        tree = linkage(v,'average')

        output = {}
        output_str = ""#str(rmsd_array) + "\n" + str(tree)
        np.savetxt("rmsd.txt", rmsd_array)
        np.savetxt("tree.txt", tree)
        output['_textarea'] = output_str

        print(json.dumps(output))
