import sys, string, locale, os, re, json
import time
import numpy as np
import math
from StringIO import StringIO

import sasmol.sasmol as sasmol

from sassie.analyze.pre.slfit import slfit
#import sassie.analyze.pre.pre_bokeh as plot
#import sassie.analyze.pre.slfit as slfit
#For local testing

#sys.path.append('./')
#import pre_bokeh as plot
#from slfit import slfit

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

def getResidueNumbers (inclusion_text, exclusion_text):
    inclusion_set = set(inputToList(inclusion_text))
    exclusion_set = set(inputToList(exclusion_text))
    residue_numbers = np.array(list(inclusion_set - exclusion_set))

    return residue_numbers


def pre_core(self, app, frame, plotQueues):

    mvars = self.mvars
    frame = frame

    ratiofile = mvars.ratiofile
    pdbfile = mvars.pdbfile
    T2dia = mvars.T2dia
    Htime = mvars.Htime
    freq  = mvars.freq
    TAUc  = mvars.TAUc
    argv_io_string = StringIO(sys.argv[1])
    json_variables = json.load(argv_io_string)
    scaling = float(json_variables["scaling"])
    model = float(json_variables["model"])

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


    self.log.debug('in pre_core')
    pgui = self.run_utils.print_gui


    default_pre_error = 1.0
#    read_in_data(self,frame,default_pcs_error)

    ratio = np.loadtxt(ratiofile)

    coord_all_raw, coord_h_raw, pre_raw = read_in_data(self,frame,default_pre_error)
    chainID = pre_raw[1][1]

    output_dir = os.path.join(mvars.runname, app)

    mvars.file_calc_pre[frame] = os.path.join(output_dir, mvars.runname+'_' + "calc_pre_" + str(frame+1).zfill(5)) + '.txt'
    mvars.file_pymol[frame] = os.path.join(output_dir, mvars.runname+'_' + "disp_SL_" + str(frame+1).zfill(5)) + '.py'
    mvars.file_summary[frame] = os.path.join(output_dir, mvars.runname+'_' + "summary_" + str(frame+1).zfill(5)) + '.txt'

    script,div,position,ratio_dist,Chi2,qR_fact,corrcoef = slfit(self, app, ratio, T2dia, Htime, freq, TAUc, pdbfile, scaling = scaling, model = model, reslist = reslist, chainID = chainID)

    #PRE Summary File
    summary_file = open(mvars.file_summary[frame], "w")
    pdb_file = mvars.pdbfile.split("/")[-1]
    dcd_file = mvars.dcdfile.split("/")[-1]
    PRE_file = mvars.ratiofile.split("/")[-1]

    exclusion_file = 'Not used'

    output = ""
    output += ("1. PDB Input File: " + pdb_file)
    output += "\n"
    output += ("2. DCD Input File: " + dcd_file)
    output += "\n"
    output += ("3. PRE Input File: " + PRE_file)
    output += "\n"
    output += ("4. T2_dia time [ms]: " + str(T2dia))
    output += "\n"
    output += ("5. Time during experiment [ms]: " + str(Htime))
    output += "\n"
    output += ("6. Spectrometer Frequency [MHz]: " + str(freq))
    output += "\n"        
    output += ("7. XYZ Position Vector [Angstrom]: [" + str.format('{0:.3f}',position[0]) + "," + str.format('{0:.3f}',position[1]) + "," + str.format('{0:.3f}',position[2]) +"]") 
    output += "\n"
    output += ("8. Chi square: " + str.format('{0:.3f}',Chi2))
    output += "\n"
    output += ("9. qR factor: " + str.format('{0:.3f}',qR_fact))
    output += "\n"
    output += ("10. Corr coeff: " + str.format('{0:.3f}',corrcoef))
    output += "\n"

    summary_file.write(output)
    summary_file.close()

    #pymol file
    pymol_script = open(mvars.file_pymol[frame], "w")

    output = write_pymol(999, position)

    pymol_script.write(output)
    pymol_script.close()

    save_pre(self, ratio_dist, frame, app, pre_raw)

    plotQueues['bokeh_plot_1'].put(script)


# handle exclusion
#    try:
#        slfit (ratio, 50e-3, 5e-3, 600.13, 4.5, "1D3Z_f.pdb")
#        slfit (ratio, T2dia, Htime, freq, TAUc, pdbfile)
#    except:
#        message = "Error in PRE module."
#        pgui(message)
#        exit()

#    generate_plot(self, app, plotQueues)

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

def read_in_data(self,frame,default_pre_error):

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

    pgui = self.run_utils.print_gui

    pre_type = []
    local_pre = []

    pre_type_id = -1
    count_pre = 0

    pre_in = open(mvars.ratiofile,'rU').readlines()
    count = 0
    for line in pre_in:
        lin = string.split(line)
        if (lin[0][0] == "#" and len(lin[0]) > 1):
            lin[0] = lin[0].replace("#","")
            if (lin[0] == ""):
                lin[0] = "#"
            else:
                lin.insert(0,"#")
        if (lin[0] == "#"):
            pre_type_id += 1
            chain_id, atom_name_1 = lin[1:3]
            pre_type.append([pre_type_id, chain_id, atom_name_1])
        else:
            if (len(lin) == 2):   # if pcs_error doesn't exist, use default_pcs_error
                local_pre.append([pre_type_id, chain_id, int(locale.atof(lin[0])), atom_name_1, float(locale.atof(lin[1])), default_pre_error])
            if (len(lin) > 2 ):
                local_pre.append([pre_type_id, chain_id, int(locale.atof(lin[0])), atom_name_1, float(locale.atof(lin[1])), float(locale.atof(lin[2])) ])
            count += 1


    local_exclusion_list=[]
#    len_pcs_selected = len(pcs_selected)
    message, pre_selected = check_pre_input(self,local_pre,local_exclusion_list)

    # Get vectors from 1st frame of input pdb or dcd
    coord_all, coord_h = get_vector_from_pdb(self,pre_selected,frame)
    pre_raw = np.asarray(pre_selected)

    return coord_all, coord_h, pre_raw

def check_pre_input(self,pre_array, exclusion_array):

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

    pre_input = pre_array
    exclusion_input = exclusion_array

    len_pre_input = len(pre_input)
    len_exclusion_input = len(exclusion_input)

    pre_selected = []
#    reslist_selected = []

    if (len_exclusion_input == 0):
        for i in range(0,len_pre_input):
            pre_selected.append(pre_input[i])
#            reslist_selected.append(pcs_input[i][:5])

    else:
        count = 0
        for i in range(0,len_pre_input):
            exclusion_exist = 0
            for j in range(0,len_exclusion_input):
                if (pre_input[i][1:3] == exclusion_input[j][1:3]):
                    if (pre_input[i][3] == exclusion_input[j][3]):
                        count += 1
                        exclusion_exist = 1
            if (exclusion_exist == 0):
                pre_selected.append(pre_input[i])
#                reslist_selected.append(pre_input[i][:5])

        if (count == 0):
            message += "could not find any residue to exclude"

    return message, pre_selected #, reslist_selected

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

def save_pre(self, ratio_dist, frame, app, pre_raw):

    message = ""
    pgui = self.run_utils.print_gui
    mvars = self.mvars

    output_dir = os.path.join(mvars.runname, app)

#   Calculate correlation coefficient and r_factor

    exp_pre = ratio_dist[:,3]
    backcalc_PRE = ratio_dist[:,1]
    res = ratio_dist[:,0]

    pre_exp_array = exp_pre.astype(np.float)

    corr_coef = np.corrcoef(pre_exp_array,backcalc_PRE)
    r_factor = calc_r(pre_exp_array,backcalc_PRE)

    fout = open(mvars.file_calc_pre[frame], 'w')

    len_pre = len(exp_pre)

# TODO: add experimental error
    fout.write("%12s %12s %12s %12s %12s %12s %12s \n" %("ID", "type", "Chain", "Residue", "PRE_exp", "PRE_calc", "delta(Exp-calc)"))

    for i in range(len_pre):
        nid = str(i + 1)
        ptype = str(pre_raw[1][3])
        chain = str(pre_raw[1][1])
        resid = str(np.int(res[i]))
        pre_exp = exp_pre[i]
        pre_calc = backcalc_PRE[i]
        pre_delta = pre_exp - pre_calc
        fout.write("%12s %12s %12s %12s %12.5f %12.5f %12.5f \n" %(nid, ptype, chain, resid, pre_exp, pre_calc, pre_delta))


    fout.close()


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

import numpy as np
import sys

def read_vertices(self):
    #--------------------------------------------------------------
    #   Given protein coordinates generate vertices of a cube around the protein
    #   INPUT: user_guess_num: integer
    #          user_guess_x, user_guess_y, user_guess_z: float array
    #--------------------------------------------------------------
    mvars = self.mvars
    npoint = mvars.user_guess_num
    x_array = mvars.user_guess_x
    y_array = mvars.user_guess_y
    z_array = mvars.user_guess_z

    vertices = np.zeros([npoint,3])
    for i in range(npoint):
        vertices[i][0] = x_array[i]
        vertices[i][1] = y_array[i]
        vertices[i][2] = z_array[i]


    return vertices

#=====================================================

def generate_vertices(self, coord=None, k_face=None):
    #--------------------------------------------------------------
    #   Given protein coordinates generate vertices of a cube around the protein
    #   INPUT:
    #       1) coord = [x,y,z]  nx3 array, n = number of residues
    #       2) k_face = 0 (default) only generate 8 corners
    #                   1 add face-points (shifted 0.7 away from the face)
    #--------------------------------------------------------------
    if k_face is None:
        k_face = 0
    #default: only corners

    if (coord.shape[1] >= 4):
        coord = np.delete(coord, 0, axis=1) #If coordinate input contains more than 3 columns, assume the first column is residue number

    coord = np.array(coord) # converts standard array to numpy array

    mx = np.max(coord, axis=0)
    mn = np.min(coord, axis=0)
    mc = np.mean(coord, axis=0)

    vertices = np.array([mn, [mn[0], mn[1], mx[2]], [mn[0], mx[1], mn[2]], [mx[0], mn[1], mn[2]], [mn[0], mx[1], mx[2]], [mx[0], mn[1], mx[2]], [mx[0], mx[1], mn[2]], mx])

    if k_face == 1:
        faces = np.array([[mn[0] - (mc[0] - mn[0]) * 0.7, mc[1], mc[2]], [mx[0] + (mx[0] - mc[0]) * 0.7, mc[1], mc[2]], [mc[0], mn[1] - (mc[1] - mn[1]) * 0.7, mc[2]], [mc[0], mx[1] + (mx[1] - mc[1]) * 0.7, mc[2]], [mc[0], mc[1], mn[2] - (mc[2] - mn[2]) * 0.7], [mc[0], mc[1], mx[2] + (mx[2] - mc[2]) * 0.7]])
        vertices = np.concatenate([vertices, faces])


    return vertices

#=====================================================

def write_pymol(resnum = None, atcoord = None, atsize = 10, atcolor = 'orange', atlabel = '', chainID = 'A', param_type = 0):
    #TODO, Cheol: ChainID should come from pdb file
    if param_type == 0:
        atelement = 'None'

    output = ""
    output += "#Python script to display SL in PyMol"
    output += "\n"
    output += "from pymol import cmd"
    output += "\n"
    output += "#paramagnetic center parameters"
    output += "\n"
    output += ("resnum = " + str(resnum))
    output += "\n"
    #output += ("coor = [" + str(atcoord[0]) + "," + str(atcoord[1]) + "," + str(atcoord[2]) +"]")
    output += ("coor = [" + str.format('{0:.3f}',atcoord[0]) + "," + str.format('{0:.3f}',atcoord[1]) + "," + str.format('{0:.3f}',atcoord[2]) +"]")
    output += "\n"
    output += ("atsize = " + str(atsize))
    output += "\n"
    output += ("chainID = \"" + str(chainID) + "\"")
    output += "\n"
    output += ("atlabel = \"" + str(atlabel) + "\"")
    output += "\n"
    output += ("atcolor = \"" + atcolor + "\"")
    output += "\n"
    output += ("atelement = \"" + atelement + "\"")
    output += "\n"
    output += "#create pseudoatom and show"
    output += "\n"
    output += "SL_atom = \"SL\""
    output += "\n"
    output += "cmd.pseudoatom(SL_atom, resi = resnum, chain = chainID, b = atsize, color = atcolor, pos = coor, elem = atelement, label = atlabel)"
    output += "\n"
    output += "cmd.show(\"spheres\",SL_atom)"

    return output
