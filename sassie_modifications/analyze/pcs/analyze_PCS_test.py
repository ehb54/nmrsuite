from pcs2tensor import pcs2tensor
from generate_vertices import generate_vertices
from write_pymol import write_pymol
import numpy as np
import os

'''
This program is a modified version of analyze_PCS
Used for testing during development in Windows environment
'''

'''
Current file names:
hydrogen_coord_file - ../testset_Tm/coordinates_H.dat
all_coord_file - ../testset_Tm/coordinates_all.dat
PCS_data_file - "../testset_Tm/PCS_data.dat"
'''

def run_PCS (pdb = None, PCS_data = None, exclusion = None):
    # handle exclusion
    try:
        #Read in the data files and coordinate files
        f = np.loadtxt("../../testset_Tm/coordinates_H.dat") # temporary read algorithm 
        coord_H = f[:,1:4]
        coord = np.loadtxt("../../testset_Tm/coordinates_all.dat") # temporary read algorithm
        pcs_raw = np.loadtxt("../../testset_Tm/PCS_data.dat")
        pcs = pcs_raw[:,1]
    except:
        print ("Error in reading the input files. Please confirm that the input has been correctly entered.")

    try:
        #Generate initial starting points for the simplex minimization (optimization)
        vertices = generate_vertices(coord, 1) #Generate guesses using all coordinates
    except:
        print ("Error in generating vertices from the input. Please confirm that the input has been correctly entered.")

    try:
        #Perform simplex minimization to determine tensor with best fit and perform other maths
        backcalc_PCS, position_vect, x_vect, Chi2, cond_num, euler = pcs2tensor(coord_H, pcs, vertices, 33, 1e-5)
    except Exception as e:
        print ("Error in backcalculating to create the tensor. Please confirm that the input has been correctly entered.")
        #print (e)

    #Define number of digits to show in output file
    round_digits = 2

    #File 1 - PCS_results/Plotdata File

    #Design structure of output file
    col_1 = pcs_raw[:, 0] #Residue number
    col_2 = pcs #Experimental PCS values
    col_3 = backcalc_PCS #Calculated (predicted) PCS values
    col_4 = col_2 - col_3 #Delta between experimental and calculated PCS values

    #Create PCS file array
    plotarr = np.vstack([col_1, col_2, col_3, col_4]).transpose()

    #Save PCS array to PCS output file
    np.savetxt("PCS_results.txt", plotarr, header = "Residue Number  Experimental PCS Values  Calculated PCS Values  Delta")


    #File 2 - PCS Summary File
    summary_file = open("PCS_summary.txt", "w")


    """ runname = ""
    app3 = ""
    output_dir = os.path.join(runname, app3)
    pdb_file = os.path.join(output_dir, pdb)
    PCS_file = os.path.join(output_dir, PCS_data)
    exclusion_file = os.path.join(output_dir, exclusion) """
    
    output = ""
    output += "\n"
    output += ("4. XYZ Position Vector: " + str(position_vect))
    output += "\n"
    output += ("5. Tensor Vector: " + str(x_vect))
    output += "\n"
    output += ("6. Chi2: " + str(Chi2))
    output += "\n"
    output += ("7. Condition Number: " + str(cond_num))
    output += "\n"
    output += ("8. Euler Angles: " + euler)

    summary_file.write(output)
    summary_file.close()
    
    '''
    inputs and where they are found
    other options that are checked
    position (XYZ)
    tensor
    chi2
    cond num
    euler angles

    for i in range(len(output)):
        output[i] = str(i) + .
    '''

    #File 3 - PyMol File
    

    pymol_script = open("disp_SL.py", "w")

    output = write_pymol(999, position_vect)

    pymol_script.write(output)
    pymol_script.close()

run_PCS()