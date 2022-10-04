import sys, string, locale, os, re, json
import time
import numpy as np
import math

import sasmol.sasmol as sasmol

from sassie.analyze.pcs.pcs2tensor import pcs2tensor
from sassie.analyze.pcs.generate_vertices import generate_vertices
import sassie.analyze.pcs.pcs_bokeh_4plot as plot 
#import sassie.analyze.pcs.pcs_bokeh as plot
from sassie.analyze.pcs.write_pymol import write_pymol

#For local testing
#import pcs_bokeh as plot 
#from pcs2tensor import pcs2tensor
#from generate_vertices import generate_vertices
#from write_pymol import write_pymol

#sys.path.append('./')

def pcs_core(self, app, frame, plotQueues):

    mvars = self.mvars
    frame = frame

    cond_num_cutoff = mvars.cond_num_cutoff
    tolerance = mvars.tolerance
    self.log.debug('in pcs_core')
    pgui = self.run_utils.print_gui

    default_pcs_error = 1.0

    read_in_data(self,frame,default_pcs_error)

    # handle exclusion
    try:
        coord_all_raw, coord_h_raw, pcs_raw = read_in_data(self,frame,default_pcs_error)
        coord_all = coord_all_raw[:,0:3].astype(np.float)
        coord_h = coord_h_raw[:,4:7].astype(np.float)
        pcs = pcs_raw[:,4].astype(np.float)
        #hydrogenCoordFile = "../testset_Tm/coordinates_H.dat"
        #allCoordFile = "../testset_Tm/coordinates_all.dat"
        #PCSDataFile = "../testset_Tm/PCS_data.dat"
        #Read in the data files and coordinate files
        #f = np.loadtxt(hydrogenCoordFile) # temporary read algorithm 
        #coord_H = f[:,1:4]
        #coord = np.loadtxt(allCoordFile) # temporary read algorithm
        #pcs_raw = np.loadtxt(PCSDataFile)
        #pcs = pcs_raw[:,1]
    except:
        message = "Error in reading the input files. Please confirm that the input has been correctly entered."
        pgui(message)
        exit()
    try:
        #Generate initial starting points for the simplex minimization (optimization)
        vertices = generate_vertices(coord_all, 1) #Generate guesses using all coordinates
    except:
        
        message = "Error in generating vertices from the input. Please confirm that the input has been correctly entered."
        pgui(message)
        exit()
    try:
        pgui("\n<<<<  Model " + str(frame + 1) + " calculation starts >>>>\n")
        #Perform simplex minimization to determine tensor with best fit and perform other maths
        backcalc_PCS, position_vect, x_vect, Chi2, cond_num, euler = pcs2tensor(self,coord_h, pcs, vertices, cond_num_cutoff,tolerance) # , 0)
    except:
        message = "Error in backcalculating to create the tensor. Please confirm that the input has been correctly entered."
        pgui(message)
        exit()

#    backcalc_PCS, position_vect, x_vect, Chi2, cond_num, euler = pcs2tensor(self,coord_h, pcs, vertices, cond_num_cutoff,tolerance)
    #Design structure of output file
#    col_1 = pcs_raw[:, 2].astype(float) #Residue number
#    col_2 = pcs #Experimental PCS values
#    col_3 = backcalc_PCS #Calculated (predicted) PCS values
#    col_4 = col_2 - col_3 #Delta between experimental and calculated PCS values

    #Plot the results
#    plotarr = np.vstack([col_1, col_2, col_3, col_4]).transpose()
#    plotarr = np.round(plotarr, round_digits)

    #Save PCS array to PCS output file

    sl_vector = np.copy(coord_h)
    for i in range(3):
        sl_vector[:,i] -=  position_vect[i]

    mvars.sl_distance = np.sqrt(np.sum(np.square(sl_vector), axis = 1))
    
    output_dir = os.path.join(mvars.runname, app)

    mvars.file_calc_pcs[frame] = os.path.join(output_dir, mvars.runname+'_' + "calc_pcs_" + str(frame+1).zfill(5)) + '.txt'
    mvars.file_summary[frame] = os.path.join(output_dir, mvars.runname+'_' + "summary_" + str(frame+1).zfill(5)) + '.txt'
    mvars.file_pymol[frame] = os.path.join(output_dir, mvars.runname+'_' + "disp_SL_" + str(frame+1).zfill(5)) + '.py'

# TODO: SL distance files are needed?
#    mvars.file_sl_distance[frame] = os.path.join(output_dir, mvars.runname+'_' + "SL_distance" + str(frame+1).zfill(5)) + '.txt'

#    np.savetxt(mvars.file_calc_pcs, plotarr, header = "Residue Number  Experimental PCS Values  Calculated PCS Values  Delta")

    save_pcs(self, pcs_raw, backcalc_PCS, frame, app)

    #File 2 - PCS Summary File
    summary_file = open(mvars.file_summary[frame], "w")
    pdb_file = mvars.pdbfile.split("/")[-1]
    dcd_file = mvars.dcdfile.split("/")[-1]
    PCS_file = mvars.pcs_input_file.split("/")[-1]
    if mvars.residue_exclusion_file_flag:
        exclusion_file = mvars.residue_exclusion_file.split("/")[-1] 
    else: 
        exclusion_file = 'Not used'

    output = ""
    output += ("1. PDB Input File: " + pdb_file)
    output += "\n"
    output += ("2. DCD Input File: " + dcd_file)
    output += "\n"
    output += ("3. PCS Input File: " + PCS_file)
    output += "\n"
    output += ("4. Input Exclusion File: " + exclusion_file)
    output += "\n"
    output += ("5. Condition number cutoff: " + str(cond_num_cutoff))
    output += "\n"
    output += ("6. Tolerance: " + str(tolerance))
    output += "\n"
    output += ("7. XYZ Position Vector [Angstrom]: " + str(position_vect))
    output += "\n"
    output += ("8. Tensor Vector: " + str(x_vect))
    output += "\n"
    output += ("9. Chi square: " + str(Chi2))
    output += "\n"
    output += ("10. Condition Number: " + str(cond_num))
    output += "\n"
    output += ("11. Euler Angles [degrees]: " + euler)

    summary_file.write(output)
    summary_file.close()

    #File 3 - PyMol File


    pymol_script = open(mvars.file_pymol[frame], "w")

    output = write_pymol(999, position_vect)

    pymol_script.write(output)
    pymol_script.close()

    if (frame == mvars.number_of_frames -1 ):
        generate_plot(self, app, plotQueues)

    return

def get_vector_from_pdb(self,pcs_array,frame_input):

    """ 
    Functionality: Find unit vectors from an user input pdb file 
    type   pcs_array  : list  
                        [pcs_type_id, chain_id, resid, atom_name_1, pcs, pcs_error]
    @rtype vector_coor: numpy array 
    @return           : unit vector of protein in input pdb file    
    """

    frame = frame_input
    pgui = self.run_utils.print_gui

    mvars = self.mvars
    mol = sasmol.SasMol(0)
    mol.read_pdb(mvars.dcdfile)
    natoms = mol.natoms()
    coor = mol.coor()
#    local_nh_vector_coor = numpy.zeros((numres_vnh,4))

    len_pcs = len(pcs_array)

#    coord_h = numpy.zeros((len_pcs,7)) 
    coord_h = [[] for i in range(len_pcs)]
    coord_all = np.zeros((natoms,3))

    # Get coordinates for all atoms
    for i in range(0,natoms):
        atom_name = mol.name()[i]

        for j in range(0,3):
            coord_all[i][j] = round(coor[frame, i, j],3) 

        for count in range(0,len_pcs):
            chain_id, res_id, type1 = pcs_array[count][1:4]
            if (mol.chain()[i] == chain_id and mol.resid()[i] == res_id):
                if (atom_name == type1 ):
                    for j in range(4):
                        #print(coord_h[count][j],  pcs_array[count][j])
                        coord_h[count].append(pcs_array[count][j])
                    # coord_h[count][0:4] = pcs_array[count][0:4]
                    for j in range(0,3):
                        #coord_h[count][j+4] = coor[frame,i,j]
                        coord_h[count].append(round(coor[frame,i,j],3))
                        
    coord_h_np = np.asarray(coord_h)

    return coord_all, coord_h_np
 
def read_in_data(self,frame,default_pcs_error):

    """
    Function:
    Read user input PCS file, input residue exclusion file

    1.Input PCS file format: Experimental PCS data file
      Header should start with # or >
      # chain_id  atom_name_1
        resid   PCS PCS_Error(optional, default = 1Hz)
         :
        resid   PCS PCS_Error(optional, default = 1Hz)
      # chain_id atom_name_1
        resid
         :
        resid
    2. Residue exclusion file: List of residues to exclude from PCS analysis (Optional input)
       # chain_id atom_name_1 
       list (1 * N_list)
       # chain_id atom_name_1
       list (1 * N_list)

    return
    @rtype coord_all = [x,y,z] for all atoms in pdb file 
    @rtype coord_h   = [pcstype, chain, resid, atomname, x, y, z ] for all atoms defined in the input pcs and exclusion files
    @rtype pcs_raw   = [pcstype, chain, resid, atomname, pcs, error ] for all atoms defined in the input pcs and exclusion files
    """

    mvars = self.mvars
    frame = frame

    if mvars.residue_exclusion_file_flag:
        mvars.residue_list = open(mvars.residue_exclusion_file,'rU').readlines()

    pgui = self.run_utils.print_gui

    pcs_type = []
    local_pcs = []
   
    pcs_type_id = -1 
    count_pcs = 0

    pcs_in = open(mvars.pcs_input_file,'rU').readlines()
    count = 0
    for line in pcs_in:
        lin = string.split(line)
        if (lin[0][0] == "#" and len(lin[0]) > 1):
            lin[0] = lin[0].replace("#","")
            if (lin[0] == ""):
                lin[0] = "#"
            else:
                lin.insert(0,"#")            
        if (lin[0] == "#"):
            pcs_type_id += 1 
            chain_id, atom_name_1 = lin[1:3] 
            pcs_type.append([pcs_type_id, chain_id, atom_name_1])
        else:
            if (len(lin) == 2):   # if pcs_error doesn't exist, use default_pcs_error 
                local_pcs.append([pcs_type_id, chain_id, int(locale.atof(lin[0])), atom_name_1, float(locale.atof(lin[1])), default_pcs_error])
            if (len(lin) > 2 ):
                local_pcs.append([pcs_type_id, chain_id, int(locale.atof(lin[0])), atom_name_1, float(locale.atof(lin[1])), float(locale.atof(lin[2])) ])
            count += 1

    # Read optional resiude list file to select resiudes for pcs analysis 
    local_residue_list=[]
#    residue_list_type=[]
    local_exclusion_list=[]

    type_id = -1
    if mvars.residue_exclusion_file_flag:
        local_residue_list = []
        for line in mvars.residue_list:
            lin = string.split(line)
            if (lin[0][0] == "#" and len(lin[0]) > 1):
                lin[0] = lin[0].replace("#","")
                if (lin[0] == ""):
                    lin[0] = "#"
                else:
                    lin.insert(0,"#")

            if lin[0] == "#" :
                type_id += 1
                chain_id, atom_name_1 = lin[1:3]
            else:
                    local_exclusion_list.append([type_id,chain_id,int(locale.atof(lin[0])),atom_name_1])

        if (len(local_exclusion_list) == 0):
            message = "The exclusion list did not contain any list of residues. Edit the exclusion file and try again."
            pgui(message)
            sys.exit()

    message, pcs_selected = check_pcs_input(self,local_pcs,local_exclusion_list)

    if (len(message) > 0):
        pgui(message)
        sys.exit()

#    len_pcs_selected = len(pcs_selected) 

    # Get vectors from 1st frame of input pdb or dcd
    coord_all, coord_h = get_vector_from_pdb(self,pcs_selected,frame)
    pcs_raw = np.asarray(pcs_selected)

    return coord_all, coord_h, pcs_raw

def check_pcs_input(self,pcs_array, exclusion_array):

    """
    Function
    Find residue ids that exist both at input pcs file and residue list file
    ==> Modified to exclude residues given in residue list file 
    @type:  pcs_array
            [pcs_type_id, chain_id, resid, atom_name_1, pcs, pcs_error]
    @type:  exclusion_array
            [pcs_type_id, chain_id, resid]
    @rtype: pcs_selected 
            [pcs_type_id, chain_id, resid, atom_name_1, pcs, pcs_error]
    @rtype: reslist_selected
            [pcs_type_id, chain_id, resid]
    """
    mvars = self.mvars

    message = ""

    # Input variables

    pcs_input = pcs_array
    exclusion_input = exclusion_array

    len_pcs_input = len(pcs_input)
    len_exclusion_input = len(exclusion_input)

    pcs_selected = []
#    reslist_selected = []

    if (len_exclusion_input == 0):
        for i in range(0,len_pcs_input):
            pcs_selected.append(pcs_input[i])
#            reslist_selected.append(pcs_input[i][:5])

    else:
        count = 0
        for i in range(0,len_pcs_input):
            exclusion_exist = 0
            for j in range(0,len_exclusion_input):
                if (pcs_input[i][1:3] == exclusion_input[j][1:3]):
                    if (pcs_input[i][3] == exclusion_input[j][3]): 
                        count += 1
                        exclusion_exist = 1 
            if (exclusion_exist == 0):
                pcs_selected.append(pcs_input[i])
#                reslist_selected.append(pcs_input[i][:5])

        if (count == 0):
            message += "could not find any residue to exclude"

    return message, pcs_selected #, reslist_selected

def generate_plot(self, app, plotQueues):

    '''
    objects for bokeh:

        bokeh_plot_1: mvars.script_first_frame, mvars.div_first_frame : Plot of first frame

    If multi-frame case:

        bokeh_plot_2: mvars.script_last_frame, mvars.div_last_frame : Plot of first frame

        bokeh_plot_3: mvars.script_multi_frame, mvars.div_multi_frame: Plot for values as a function of frame
    '''

    pgui = self.run_utils.print_gui
    mvars = self.mvars

    # for frame = 0
    f = 0
    message = ""

#    try:
#        mvars.script_first_frame, mvars.div_first_frame = plot.single_model_plot(self,f,400,500,app)
#    except:
#        message += "PCS failed to generate first frame plot."
#        pgui(message)
#        sys.exit()

    mvars.script_first_frame, mvars.div_first_frame = plot.single_model_plot(self,f,400,500,app)

    plotQueues['bokeh_plot_1'].put(mvars.script_first_frame)

#    print json.dumps(plotQueues['bokeh_plot_1'])

    if (mvars.number_of_frames > 1):
        f = mvars.number_of_frames - 1
        try:
            mvars.script_last_frame, mvars.div_last_frame = plot.single_model_plot(self,f,400,500,app)
        except:
            message += "PCS failed to generate last frame plot."
            pgui(message)
            sys.exit()
        plotQueues['bokeh_plot_2'].put(mvars.script_last_frame)

        try:
            mvars.script_multi_frame, div_multi_frame = plot.multi_model_plot(self,400,500,app)
        except:
            message += "PCS failed to generate multi_frame plot."
            pgui(message)
            sys.exit()
# TODO: Decide which properties need to be monitored as a function of frame
#        plotQueues['bokeh_plot_3'].put(mvars.script_multi_frame)

    time.sleep(1)

    return

def save_pcs(self, pcs_raw, backcalc_PCS, frame, app):

    message = ""
    pgui = self.run_utils.print_gui
    mvars = self.mvars

    output_dir = os.path.join(mvars.runname, app)

#   Calculate correlation coefficient and r_factor

    pcs_exp_array = pcs_raw[:,4].astype(np.float)

    corr_coef = np.corrcoef(pcs_exp_array,backcalc_PCS)
    r_factor = calc_r(pcs_exp_array,backcalc_PCS)

    fout = open(mvars.file_calc_pcs[frame], 'w')

#    if (mvars.number_of_frames > 1):
    if (frame == 0):
        mvars.file_out_summary = os.path.join(output_dir, mvars.runname+'_'+'results_per_model')+'.txt'
        fout_stack = open(mvars.file_out_summary, 'w')
#        fout_stack.write("%12s %12s %12s %12s %12s %12s %12s %12s %12s \n" \
#            %("Model", "Corr_coeff", "R-factor", "Axx[*10^-3]", "Ayy[*10^-3]", "Azz[*10^-3]", "Alpha", "Beta", "Gamma" ))
        fout_stack.write("%12s %12s %12s \n" \
            %("Model", "Corr_coeff", "R-factor"))

    else:
        fout_stack = open(mvars.file_out_summary, 'a')

    len_pcs = len(pcs_raw) 

# TODO: add experimental error
    fout.write("%12s %12s %12s %12s %12s %12s %18s \n" %("ID", "type", "Chain", "Residue","PCS_exp", "PCS_calc", "delta(Exp-calc)"))

    for i in range(len_pcs):
        nid = str(i + 1)
        ptype = str(pcs_raw[i][3])
        chain = str(pcs_raw[i][1])
        resid = str(pcs_raw[i][2])
        pcs_exp = float(pcs_raw[i][4])
        pcs_calc = backcalc_PCS[i]
        pcs_delta = pcs_exp - pcs_calc
        fout.write("%12s %12s %12s %12s %12.5f %12.5f %12.5f \n" %(nid, ptype, chain, resid, pcs_exp, pcs_calc, pcs_delta))

    fout.close()

#    fout_stack.write("%12s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n" \
#        %(int(frame+1), corr_coef[0][1], r_factor, s_eigenvalue[0]*sfactor, s_eigenvalue[1]*sfactor, s_eigenvalue[2]*sfactor, \
#        e_ang[0],e_ang[1],e_ang[2]))
#    fout_stack.close()

    fout_stack.write("%12s %12.5f %12.5f \n" \
        %(int(frame+1), corr_coef[0][1], r_factor ))
    fout_stack.close()

    return

def calc_r(exp,calc):

    '''
    Calculate correlation coefficeint between experimental and calculated RDCs.
    @type  exp : list
    @param exp : Experimental input RDCs
    @type  calc : list
    @param calc : Calculated RDCs
    @rtype  r_coeff : float
    @return r_coeff : Correlation coeff.
    '''

    dobs_dcalc = exp - calc
    num = np.matmul(dobs_dcalc,dobs_dcalc)
    num = np.mean(num)

    denom = 2.*np.mean(np.matmul(exp,exp))

    r_coeff = math.sqrt(num/denom)

    return r_coeff

