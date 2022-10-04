#from sassie.pcs2tensor import pcs2tensor
#from sassie.generate_vertice import generate_vertices
#from sassie.write_pymol import write_pymol

from pcs2tensor import pcs2tensor
from generate_vertices import generate_vertices
from write_pymol import write_pymol
import numpy as np
import os

'''
This program performs PCS data analysis
Input files: 
    1) pdb file (all coordinates and hydrogen coordinates)
    2) Edited PCS data (columns contain: residue number, shifts, error)
    3) List of data points to be excluded
    4) Condition number cutoff
    5) Tolerances (in X and in function)
Calls generate_vertices to determine starting points for simplex
    -Starting points determined using all coordinates
Calls pcs2tensor to calculate tensor and other properties
    1) Calculates tensor using 3 parameters (x,y,z coordinates)
    2) Performs optimization of 8 parameters (coordinates and tensor)
    3) Calculates tensor eigenvalues
    4) Calculates axial and rhombic components
    5) Calculates euler angles
    6) Displays standard outputs of the 5 components above
Outputs a text file with columns containing:
    1) Residue
    2) Experimental PCS values
    3) Calculated PCS values
    4) Delta between experimental and calculated PCS values
'''


def run_PCS (pdb = None, PCS_data = None, exclusion = None, cond_num_cutoff = None, tolerance = None):
    # handle exclusion
    # set defaults from advanced options
    if cond_num_cutoff is None:
        cond_num_cutoff = 33 #default cond_num_cutoff
    if tolerance is None:
        tolerance = 1e-02 #default tolerance
    try:
        #Read in the data files and coordinate files, can be replaced
        f = np.loadtxt(hydrogen_coord_file) # temporary read algorithm 
        coord_H = f[:,1:4]
        coord = np.loadtxt(all_coord_file) # temporary read algorithm
        pcs_raw = np.loadtxt(PCS_data_file)
        pcs = pcs_raw[:,1]
        #Note: must add error column with values of 1 if error column is missing
    except:
        print ("Error in reading the input files. Please confirm that the input has been correctly entered.")
        #print ("Coord_H Dimensions:", coord_H.shape)
        #print ("Coord Dimensions:", coord.shape)
        #print ("PCS Dimensions:", pcs_raw.shape)
    try:
        #Generate initial starting points for the simplex minimization (optimization)
        vertices = generate_vertices(coord, 1) #Generate guesses using all coordinates
    except:
        print ("Error in generating vertices from the input. Please confirm that the input has been correctly entered.")

    try:
        #Perform simplex minimization to determine tensor with best fit and perform other maths
        backcalc_PCS, position_vect, x_vect, Chi2, cond_num, euler = pcs2tensor(coord_H, pcs, vertices, cond_num_cutoff, tolerance)
    except:
        print ("Error in backcalculating to create the tensor. Please confirm that the input has been correctly entered.")

    #Define number of digits to show in output file


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
    runname = ""
    app3 = ""
    output_dir = os.path.join(runname, app3)
    pdb_file = os.path.join(output_dir, pdb)
    PCS_file = os.path.join(output_dir, PCS_data)
    exclusion_file = os.path.join(output_dir, exclusion)
    
    output = ""
    output += ("1. PDB Input File: " + pdb_file)
    output += "\n"
    output += ("2. PCS Input File: " + PCS_file)
    output += "\n"
    output += ("3. Exclusion File: " + exclusion_file)
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
    '''

    '''
    for i in range(len(output)):
        output[i] = str(i) + .
    '''


    #File 3 - PyMol File
    
    pymol_script = open("disp_SL.py", "w")
    runname = ""
    app = ""
    output_dir = os.path.join(runname, app)

    output = write_pymol(999, position_vect)

    pymol_script.write(output)
    pymol_script.close()

run_PCS()