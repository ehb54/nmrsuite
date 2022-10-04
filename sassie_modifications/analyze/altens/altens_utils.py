import sys, string, locale, os, re, json
import time
import numpy
import math

import sasmol.sasmol as sasmol
#import sassie.analyze.altens.altens_bokeh as plot 
import altens_bokeh as plot

sys.path.append('./')

def get_nh_vector(self,rdc_array,frame_input):

    """ 
    Functionality: Find unit vectors from an user input pdb file 
    TODO: Need to make this work for atom pairs defined in header section input RDC file in future 
    == > Done
    Nov. 30. 2018: Resiude list file is for exclusion of residues for ALTENS analysis 
    type   rdc_array: list  
    @rtype vector_coor: numpy array 
    @return : unit vector of protein in input pdb file    
    """

    frame = frame_input
    pgui = self.run_utils.print_gui

    mvars = self.mvars
    mol = sasmol.SasMol(0)
    mol.read_pdb(mvars.dcdfile)
    natoms = mol.natoms()
    coor = mol.coor()
#    local_nh_vector_coor = numpy.zeros((numres_vnh,4))

    len_rdc = len(rdc_array)

    vector_list = numpy.zeros((len_rdc,4)) 
    vector_coor = numpy.zeros((len_rdc,4)) 

    # Check atom naming convention in rdc and pdb files

    for i in range(0,natoms):
        atom_name = mol.name()[i]
        for count in range(0,len_rdc):
            chain_id, res_id, type1, type2 = rdc_array[count][1:5]
            if (mol.chain()[i] == chain_id and mol.resid()[i] == res_id):
                if (atom_name == type1 ):
                    vector_list[count][0] = i
                    vector_list[count][2] = 1    # if atom_name in rdc input exists in pdb file
                if (atom_name == type2 ):
                    vector_list[count][1] = i
                    vector_list[count][3] = 1
                # if (mol.resname()[i] == "GLY" and atom_name[0] == "CA" ):
                # Case for CH vector in GLY 
                if ( mol.resname()[i] == "GLY" and (type1 == "CA" or type2 == "CA")):     
                    if ( atom_name == "CA" ): 
                        vector_list[count][0] = i
                        vector_list[count][2] = 1
                    if (atom_name == "1HA" or atom_name == "HA1"):
                        vector_list[count][1] = i
                        vector_list[count][3] = 1

    for i in range(0,len_rdc):
        if (vector_list[i][2]*vector_list[i][3] == 0):
            message = "Check input rdc file if atom types in a residue " + str(rdc_array[i][1:5]) + " are correct or remove this residue from input rdc file"
            #print(rdc_array[i][1:5])
            pgui(message)
            sys.exit()
                          
    for i in range(0,len_rdc):
        n1 = int(vector_list[i][0])
        n2 = int(vector_list[i][1])
        atom_name_i = mol.name()[n1]
        atom_name_j = mol.name()[n2]
        xx = coor[frame,n2,0] - coor[frame,n1,0]
        yy = coor[frame,n2,1] - coor[frame,n1,1]
        zz = coor[frame,n2,2] - coor[frame,n1,2]
        rr = numpy.sqrt(xx*xx + yy*yy + zz*zz)
        xx = xx/rr
        yy = yy/rr
        zz = zz/rr
#        print (n1,n2, atom_name_i, atom_name_j,coor[frame,n1, 0], coor[frame,n2,0])       
        vector_coor[i][0] = i
#        vector_coor[i][0] = mol.resid()[n1] 
        vector_coor[i][1] = xx
        vector_coor[i][2] = yy
        vector_coor[i][3] = zz

    return vector_coor
 
def read_in_data(self,frame,default_rdc_error):

    """
    Function:
    Read user input RDC file, input residue list file, and
    Set up random seed if use_monte_carlo_flag = True
    Note rdc data is normalized by dmax defined in get_dmax

    1.Input RDC file format: Experimental RDC data file
      Header should start with # or >
      # chain_id  atom_name_1 atom_name_2
        resid   RDC RDC_Error(optional, default = 1Hz)
         :
        resid   RDC RDC_Error(optional, default = 1Hz)
      # chain_id atom_name_1 atom_name_2
        resid
         :
        resid
    2. Residue list file: List of residues to exclude from RDC analysis (Optional input)
       # chain_id atom_name_1 atom_name_2
       list (1 * N_list)
       # chain_id atom_name_1 atom_name_2
       list (1 * N_list)
    """

    mvars = self.mvars
    frame = frame

    mvars.rdc = open(mvars.rdc_input_file,'rU').readlines()    
    mvars.dmax_array = []
    mvars.di_error = []

    if mvars.residue_list_file_flag:
        mvars.residue_list = open(mvars.residue_list_file,'rU').readlines()

    pgui = self.run_utils.print_gui

    rdc_type = []
    local_rdc = []
   
    rdc_type_id = -1 
    count_rdc = 0

    count = 0
    for line in mvars.rdc:
        lin = string.split(line)
        if (lin[0][0] == "#" and len(lin[0]) > 1):
            lin[0] = lin[0].replace("#","")
            if (lin[0] == ""):
                lin[0] = "#"
            else:
                lin.insert(0,"#")            
        if (lin[0] == "#"):
            rdc_type_id += 1 
            chain_id, atom_name_1, atom_name_2 = lin[1:4] 
            rdc_type.append([rdc_type_id, chain_id, atom_name_1, atom_name_2])
        else:
            if (len(lin) == 2):   # if rdc_error doesn't exist, use default_rdc_error 
                local_rdc.append([rdc_type_id, chain_id, int(locale.atof(lin[0])), atom_name_1, atom_name_2, locale.atof(lin[1]), default_rdc_error])
            if (len(lin) > 2 ):
                local_rdc.append([rdc_type_id, chain_id, int(locale.atof(lin[0])), atom_name_1, atom_name_2, locale.atof(lin[1]), locale.atof(lin[2])])
            count += 1

    # Read optional resiude list file to select resiudes for rdc analysis 
    local_residue_list=[]
#    residue_list_type=[]
    local_exclusion_list=[]

    type_id = -1
    if mvars.residue_list_file_flag:
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
                chain_id, atom_name_1, atom_name_2 = lin[1:4]
            else:
                    local_exclusion_list.append([type_id,chain_id,int(locale.atof(lin[0])),atom_name_1,atom_name_2])

        if (len(local_exclusion_list) == 0):
            message = "The exclusion list did not contain any list of residues. Edit the input file and try again."
            pgui(message)
            sys.exit()

    message, rdc_selected, residue_list_selected = check_rdc_input(self,local_rdc,local_exclusion_list)
    if (len(message) > 0):
        pgui(message)
        sys.exit()

    len_rdc_selected = len(rdc_selected) 

#    print (rdc_selected)

    # Get vectors from 1st frame of input pdb or dcd
    mvars.vcoor = get_nh_vector(self,rdc_selected,frame)

    # Converting list to numpy array
    di = numpy.zeros((len_rdc_selected,3))
    for i in range(0,len_rdc_selected):
        di[i][0] = i
#        dmax, err, message = get_dmax(local_rdc[i][3], local_rdc[i][4])
        dmax, err, message = get_dmax(rdc_selected[i][3], rdc_selected[i][4])
        if (err): 
            pgui(message)
            sys.exit() 
        else:
            di[i][1] = float(rdc_selected[i][5])/dmax
            mvars.dmax_array.append(dmax)
            mvars.di_error.append(float(rdc_selected[i][6]))
      
    mvars.di = di 
    mvars.rlist = numpy.asarray(residue_list_selected)

#    print (mvars.di_error)
#    print (mvars.dmax_array)
    # set up random seed

    if mvars.use_monte_carlo_flag:
        if mvars.seed[0] == 1:
            from numpy.random import RandomState
            mvars.seed_object = RandomState(mvars.seed[1])
        else:
            mvars.seed_object = -1

    return

def get_dmax(atom1,atom2):

    """
    Function
    Get atom types and return Dmax value if types exist.
    type atom1, atom2: string
    @rtype dmax: float
    @rtype err: boolean 
    @rtype message: string  
    """
    dmax_table = {
        "HN": 21647.48153, "NH": 21647.48153, 
        "HC": 47949.08021, "CH": 47949.08021,
        "CN":  2786.77704, "NC":  2786.77704,
        "CC":  4499.897848
        }

    err = True
    message = ""
    dmax = 0

    name_of_vector = atom1[0] + atom2[0]

    for key in dmax_table:
        if (key == name_of_vector):
            dmax =dmax_table[key] 
            err = False
            break

    return dmax,err,message
 
def check_rdc_input(self,rdc_array, exclusion_array):

    """
    Function
    Find residue ids that exist both at input rdc file and residue list file
    ==> Modified to exclude residues given in residue list file 
    @rtype: rdc_selected 
    @rtype: reslist_selected
    """
    mvars = self.mvars

    message = ""

    # Input variables

    rdc_input = rdc_array
    exclusion_input = exclusion_array

    len_rdc_input = len(rdc_input)
    len_exclusion_input = len(exclusion_input)

    rdc_selected = []
    reslist_selected = []

    if (len_exclusion_input == 0):
        for i in range(0,len_rdc_input):
            rdc_selected.append(rdc_input[i])
            reslist_selected.append(rdc_input[i][:6])

    else:
        count = 0
        for i in range(0,len_rdc_input):
            exclusion_exist = 0
            for j in range(0,len_exclusion_input):
                if (rdc_input[i][1:3] == exclusion_input[j][1:3]):
                    if (rdc_input[i][3] == exclusion_input[j][3] and rdc_input[i][4] == exclusion_input[j][4]):
                        count += 1
                        exclusion_exist = 1 
                    if (rdc_input[i][3] == exclusion_input[j][4] and rdc_input[i][4] == exclusion_input[j][3]):
                        count += 1
                        exclusion_exist = 1
            if (exclusion_exist == 0):
                rdc_selected.append(rdc_input[i])
                reslist_selected.append(rdc_input[i][:6])

        if (count == 0):
            message += "could not find any residue to exclude"

    return message, rdc_selected, reslist_selected

def dir_cos(self):

    """
    Calculate direction cosines from N-H unit vectors
    """

    mvars = self.mvars

    local_input = mvars.vcoor
    local_input_len = len(local_input)
    local_output = numpy.zeros((len(local_input),3))
 
    for ii in range(0,local_input_len):
         xx = local_input[ii][1]
         yy = local_input[ii][2]
         zz = local_input[ii][3]
         local_output[ii][0] = math.acos(xx)
         local_output[ii][1] = math.acos(yy)
         local_output[ii][2] = math.acos(zz)

    mvars.dircosangle = local_output
  
    return 

def extract_par(mat_s):

    """
    Extract 5 independent elements from Alignment tensor S
    @type  mat_s : list
    @param mat_s : Alignment tensor S
    @rtype extract_m : list
    @rturn : extract tensor components to recalculate RDCs 
    """  
    a = mat_s[1][1]
    b = mat_s[2][2]
    c = mat_s[0][1]
    d = mat_s[0][2]
    e = mat_s[1][2]

    extract_m = numpy.array([a,b,c,d,e])
    return extract_m

def gen_a(self):

    """ 
    Calculates the A matrix from the input co-ordinate
    set and obtains the singular value decomposition
    A = U.T.transpose(V). The elements of T are the
    singular values of A. The inverse of A is obtained
    as inverse(A)=V.diagonal(1/T).transpose(U).
    """
    mvars = self.mvars

    tol = 0.01
   
    local_dircosangle = mvars.dircosangle 
    local_dircosangle_len = len(local_dircosangle)
    local_mat_a = numpy.zeros((len(local_dircosangle),5))

#    local_output = numpy.ones((local_input_len,local_input_len))
    for ii in range(0,local_dircosangle_len):
        phix = math.cos(local_dircosangle[ii][0])
        phiy = math.cos(local_dircosangle[ii][1])
        phiz = math.cos(local_dircosangle[ii][2])

        local_mat_a[ii][0] = phiy*phiy - phix*phix
        local_mat_a[ii][1] = phiz*phiz - phix*phix
        local_mat_a[ii][2] = 2.*phix*phiy
        local_mat_a[ii][3] = 2.*phix*phiz
        local_mat_a[ii][4] = 2.*phiy*phiz

    mvars.mat_a = local_mat_a

    mvars.con_num = numpy.linalg.cond(local_mat_a)

    # single value decomposition 
    # numpy svd: u * t * v = a  where t is diagonal element of T matrix
    # matlab svd: u,t,v = svd(a) while u*t*v_transpose = a 

    local_u,local_t, local_v = numpy.linalg.svd(local_mat_a,full_matrices=True) 
    local_t_len = len(local_t)


    local_mat_t = numpy.zeros((local_dircosangle_len,local_t_len))

    for ii in range(0,local_t_len):
        local_mat_t[ii][ii] = local_t[ii]

    v_trans = numpy.transpose(local_v)
    u_trans = numpy.transpose(local_u)


    mvars.mat_t = local_mat_t
    mvars.mat_v_trans = v_trans

    local_s = numpy.zeros((local_t_len,local_dircosangle_len))

    local_max_t = numpy.max(numpy.max(local_t))
    tol_t = tol*local_max_t

    for ii in range(0,5):
        if local_mat_t[ii][ii] > tol_t:
           local_s[ii][ii] = 1./local_mat_t[ii][ii] 
        else:
           local_s[ii][ii] = 0.

    mat_a_inv = numpy.matmul(v_trans,local_s)

    mvars.mat_a_inv = numpy.matmul(mat_a_inv,u_trans)

    return


def x_tensor3(a_inv,d_eff):

    '''
    Calculate alignment tensor S = Inv(A)*RDC and give its eigenvalues eigenvectors.
    @type  a_inv : list
    @param a_inv : Inverse of matrix A
    @type  d_eff : list
    @param d_eff : Normalized RDCs 
    @rtype  local_s : list 
    @return local_s : Alignment tensor S
    @rtype  local_v : list
    @return local_v : Eigenvectors of alignment tensor S
    @rtype  local_d : list
    @return local_d : Eigenvalues of alignment tensor S
    '''
    local_a_inv = a_inv
    local_di = d_eff 

    shape_a_inv = local_a_inv.shape
    len_a_inv = shape_a_inv[1]

    local_q_v = numpy.zeros((5))

    local_s = numpy.zeros((3,3))
    for ii in range(0,5):
        for jj in range(0,len_a_inv):
            local_q_v[ii] = local_q_v[ii] + local_a_inv[ii][jj]*local_di[jj][1]

    local_s[0][0] = -local_q_v[0] - local_q_v[1]
    local_s[1][1] = local_q_v[0]
    local_s[2][2] = local_q_v[1]
    local_s[0][1] = local_q_v[2]
    local_s[0][2] = local_q_v[3]
    local_s[1][2] = local_q_v[4]
    local_s[2][0] = local_s[0][2]
    local_s[2][1] = local_s[1][2]
    local_s[1][0] = local_s[0][1]

### Eigen value, eigen vectors

    local_d, local_v = numpy.linalg.eig(local_s)

    local_v = numpy.real(local_v)
    local_d = numpy.real(local_d)
     
    index = numpy.argsort(numpy.abs(local_d))

    local_v = local_v[:,index]

    local_d_copy = local_d

    local_d_copy[:] = local_d[index[:]]

#    mvars.s_tensor = local_s
#    mvars.s_val = local_d_copy
#    mvars.rot = local_v
   
    return local_v, local_d_copy,local_s


def euler_angle(vt):

    '''
    Calculate Euler angles from input eigenvectors.
    @type  vt : list
    @param vt : Input eigenvectors
    @rtype  eul_angle : list
    @return eul_angle : Euler angles 
    '''
    flag = 0

    # d : rotation matrix composed of eigenvectors of alignment tensor S

    d = numpy.array([vt,vt,vt,vt,vt,vt,vt,vt])

    d[1][2][:] = -d[0][2][:]
    d[2][1][:] = -d[0][1][:]
    d[3][0][:] = -d[0][0][:]
    d[4][1][:] = -d[0][1][:]
    d[4][2][:] = -d[0][2][:]
    d[5][0][:] = -d[0][0][:]
    d[5][2][:] = -d[0][2][:]
    d[6][0][:] = -d[0][0][:]
    d[6][1][:] = -d[0][1][:]
    d[7][:][:] = -d[0][:][:]

    rin = vt

    for ii in range(0,8):
        if flag > 0:
           break
        rin = d[ii,:] 
        eul_angle,flag = euler_parm(rin)

    return eul_angle 


def euler2all(a,b,c):

    '''
    For a given set of euler angles (in degrees),
    produce all other possible angles
    based on the symmetry properties of
    the diffusion tensor
    {a,b,g}={a,b,g+180o}={a+180o,180o-b,180o-g}={a+180o,180o-b,360o-g}
    and select the set that has all three angles in the [0 180]-cube
    @type  a,b,c : float
    @param a,b,c : input euler angle (alpha, beta, gamma) in degree
    @rtype     z : list
    @return    z : Possible sets of Euler angles
    @rtype  zsel : list
    @return zsel : Selected euler angles from 0 to 180 degree
    '''

    # Assume input vector

    z = numpy.ones((4,3))

    zsel = [a,b,c]
    z[0][:] = [a, b, c]
    z[1][:] = [a, b, c+180.]
    z[2][:] = [a+180., 180.-b, 180.-c]
    z[3][:] = [a+180., 180.-b, 360.-c]

    # Convert angles to normal range

    for ii in range(0,4):
        for jj in range(0,3):
            if z[ii][jj] >= 270.:
               z[ii][jj] = z[ii][jj] -360.
            if z[ii][jj] <= -270.:
               z[ii][jj] = z[ii][jj] +360.

    # Select angles 

    for ii in range(0,4):
        if (z[ii][0] >= 0.) & (z[ii][0] <= 180.) & (z[ii][1] >= 0.) & (z[ii][1] <= 180.) & (z[ii][2] >= 0.) & (z[ii][2] <= 180.): 
            zsel[:] = z[ii][:]                  
      
    return z, zsel

def euler_parm(rin):

    '''
    Finds the Euler angles corresponding to the input rotationi matrix
    REMEMBER that any row multiplied by -1 is also a valid solution
    @type  rin : list
    @param rin : input rotation matrix made by eigenvectors of alignment tensor S 
    @rtype  angle : list
    @return angle : Corrected euler angle
    @rtype  zflag : integer
    @return zflag : Flag to call euler_parm in the module, euler_angle, if zflag > 0  
    '''

    r = rin

    tol = 1.0e-12
    tol_sin = 1.0e-12

    angle = numpy.zeros((3))
    zflag = 0

    # TODO: Check negative_flag is necessary. 

    negative_flag = 1

    if negative_flag > 0:
       if numpy.linalg.det(r) < 0:
          r[2][:] = -rin[2][:]

    t = math.acos(r[2][2])

    sin_t = math.sin(t)

    if abs(sin_t) > tol_sin:
       f1 = math.asin(r[2][1]/sin_t)
       f2 = math.pi - f1
       k1 = math.asin(r[1][2]/sin_t)
       k2 = math.pi - k1

       flag = 0

       flag = choose_euler(r,flag,f1,t,k1,0)
       flag = choose_euler(r,flag,f1,t,k2,1)
       flag = choose_euler(r,flag,f2,t,k1,2)
       flag = choose_euler(r,flag,f2,t,k2,3)

       result = numpy.array([[0,f1,t,k1],[1,f1,t,k2],[2,f2,t,k1],[3,f2,t,k2]])

# This may be redundant

#       for ii in range (0,4):
#           index = numpy.where(flag == result[ii][0])
#

#       print "In Euler_parm index = ", index

#       if index[0].size != 0:
#          idd = index[0][0]

       for jj in range (0,3):
           angle[jj] = 180.*(result[flag][jj+1])/math.pi

       for ii in range(0,3):
          if angle[ii] < 0.:
             angle[ii] = 360. - numpy.abs(angle[ii])
    else:

        if math.cos(t)+1 < tol:
           angle[1] = 180.
           angle[0] = (math.acos(r[0][0]))*180./math.pi
           angle[2] = 0.
        if math.cos(t)-1 < tol:
           angle[1] = 0.
           angle[0] = (math.acos(r[0][0]))*180./math.pi
           angle[2] = 0.
    if (angle[0] != 0.)&(angle[1] != 0.)& (angle[2] != 0.):
       zflag = 1

    return angle,zflag 

def choose_euler(r,flag,a,b,c,num):

    '''
    Pick an euler angle to manipulate in module, euler_param.
    @type  r : list
    @param r : Input rotation matrix 
    @type  a,b,c : float
    @param a,b,c : Input euler angles (alpha, beta, gamma) in radian
    @type  num : integer
    @param num : Index for 4 possible euler angles for a prinicpal axis
    @type  flag : integer
    @param flag : Flag to choose an index of the euler angle to be considered as a solution
    @rtype  flag_out : integer
    @return flag_out : Index of appropriate euler angle among 4 possible angles 
    '''
    tol = 1.0e-12

    flag_out = flag

    if flag == 0:
       if (math.cos(a)*math.cos(b)*math.cos(c)-math.sin(a)*math.sin(c)) - r[0][0] < tol:
          if (math.sin(a)*math.cos(b)*math.cos(c)+math.cos(a)*math.sin(c)) - r[0][1] < tol:
             if (-math.cos(a)*math.cos(b)*math.sin(c)-math.sin(a)*math.cos(c)) - r[1][0] < tol:
                if (-math.sin(a)*math.cos(b)*math.sin(c)+math.cos(a)*math.cos(c)) - r[1][1] < tol:
                   if (math.cos(a)*math.sin(b)) - r[2][0] < tol:
                      flag_out = num  
    
    return flag_out


def recalc_d(s,a):

    '''
    Calculate matrix d = mat(A)*mat(S) for test the quality of S.
    @type  s : list
    @param s : Alignment tensor S
    @type  a : list
    @type  a : Matrix A
    @rtype  new_d : numpy array 
    @return new_d : Back-calulated RDC  
    '''

    new_d = numpy.matmul(a,s)

    return new_d

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
    num = numpy.matmul(dobs_dcalc,dobs_dcalc)
    num = numpy.mean(num)

    denom = 2.*numpy.mean(numpy.matmul(exp,exp))

    r_coeff = math.sqrt(num/denom)

    return r_coeff


def calc_ani(s):

    '''
    Calculate anisotropy parameter for alignment tensor
    @type  s : list
    @param s : Alignment tensor S
    @rtype  aniso : float
    @return aniso : Anisotropy parameter
    '''
    aniso = (s[0]-s[1])/s[2]

    return aniso


def mc_table(self,d):

    mvars = self.mvars
    num = len(d[:,1])
    d_mc = numpy.zeros((num,2)) 

#   Generate random errors for RDCs 
    ran_num = mvars.seed_object.randn(num,2)
#    d_noise[:,1] = d[:,1] + ran_num[:,1]
    d_mc[:,1] = d[:,1] + ran_num[:,1]*mvars.di_error[:]/mvars.dmax_array[:]   # Correct the inclusion of dmax

    return d_mc

def initialize_average(self,new_d, s_tensor, s_eigenvalue, e_ang, corr_coef, r_factor, an_param):
    
    mvars = self.mvars
    mvars.ave_new_d = numpy.copy(new_d)
    mvars.ave_s_tensor = numpy.copy(s_tensor)
    mvars.ave_s_eigenvalue = numpy.copy(s_eigenvalue)
    mvars.ave_e_ang = numpy.copy(e_ang)
    mvars.ave_corr_coef = numpy.copy(corr_coef)
    mvars.ave_r_factor = numpy.copy(r_factor)
    mvars.ave_an_param = numpy.copy(an_param)

    mvars.std_new_d = numpy.multiply(new_d,new_d)
    mvars.std_s_tensor = numpy.multiply(s_tensor,s_tensor)
    mvars.std_s_eigenvalue = numpy.multiply(s_eigenvalue,s_eigenvalue)
    mvars.std_e_ang = numpy.multiply(e_ang,e_ang)
    mvars.std_an_param = numpy.multiply(an_param, an_param)

    return

def calc_average(self, new_d, s_tensor, s_eigenvalue, e_ang, corr_coef, r_factor, an_param, frame):

    mvars = self.mvars
    number_of_frames = mvars.number_of_frames

    if (frame != 0): 
        mvars.ave_new_d = mvars.ave_new_d + new_d
        mvars.ave_s_tensor = mvars.ave_s_tensor + s_tensor
        mvars.ave_s_eigenvalue = mvars.ave_s_eigenvalue + s_eigenvalue
        mvars.ave_e_ang = mvars.ave_e_ang + e_ang
        mvars.ave_corr_coef = mvars.ave_corr_coef + corr_coef
        mvars.ave_r_factor = mvars.ave_r_factor + r_factor
        mvars.ave_an_param = mvars.ave_an_param + an_param

        mvars.std_new_d = mvars.std_new_d + numpy.multiply(new_d,new_d)
        mvars.std_s_tensor = mvars.std_s_tensor + numpy.multiply(s_tensor,s_tensor)
        mvars.std_s_eigenvalue = mvars.std_s_eigenvalue + numpy.multiply(s_eigenvalue,s_eigenvalue)
        mvars.std_e_ang = mvars.std_e_ang + numpy.multiply(e_ang,e_ang)
        mvars.std_an_param = mvars.std_an_param + numpy.multiply(an_param, an_param)

    if (frame == number_of_frames - 1 ):
        mvars.ave_new_d = mvars.ave_new_d / float(number_of_frames)
        mvars.ave_s_tensor = mvars.ave_s_tensor / float(number_of_frames)
        mvars.ave_s_eigenvalue = mvars.ave_s_eigenvalue / float(number_of_frames)
        mvars.ave_e_ang = mvars.ave_e_ang / float(number_of_frames)
        mvars.ave_corr_coef = numpy.corrcoef(mvars.di[:,1]*mvars.dmax_array[:],mvars.ave_new_d)
        mvars.ave_r_factor = calc_r(mvars.di[:,1]*mvars.dmax_array[:],mvars.ave_new_d)       
#        mvars.ave_corr_coef = mvars.ave_corr_coef / float(number_of_frames)
#        mvars.ave_r_factor = mvars.ave_r_factor / float(number_of_frames)
        mvars.ave_an_param = mvars.ave_an_param / float(number_of_frames)
        
        mvars.std_new_d = mvars.std_new_d/float(number_of_frames) - numpy.multiply(mvars.ave_new_d,mvars.ave_new_d)
        mvars.std_s_tensor = mvars.std_s_tensor/float(number_of_frames) - numpy.multiply(mvars.ave_s_tensor,mvars.ave_s_tensor)
        mvars.std_s_eigenvalue = mvars.std_s_eigenvalue/float(number_of_frames) - numpy.multiply(mvars.ave_s_eigenvalue, mvars.ave_s_eigenvalue)
        mvars.std_e_ang = mvars.std_e_ang/float(number_of_frames) - numpy.multiply(mvars.ave_e_ang, mvars.ave_e_ang)
        mvars.std_an_param = mvars.std_an_param/float(number_of_frames) - numpy.multiply(mvars.ave_an_param, mvars.ave_an_param)

        mvars.std_new_d = numpy.sqrt(mvars.std_new_d)
        mvars.std_s_tensor = numpy.sqrt(mvars.std_s_tensor)
        mvars.std_s_eigenvalue = numpy.sqrt(mvars.std_s_eigenvalue)
        mvars.std_e_ang = numpy.sqrt(mvars.std_e_ang)
        mvars.std_an_param = numpy.sqrt(mvars.std_an_param)

#   Check Average Correlation Coeff
#        test_corr_coef = numpy.corrcoef(mvars.di[:,1]*mvars.dmax_array[:],mvars.ave_new_d)
#        test_r_factor = calc_r(mvars.di[:,1]*mvars.dmax_array[:],mvars.ave_new_d)
#        print (test_corr_coef[0][1], test_r_factor )
#        print (mvars.ave_corr_coef[0][1], mvars.ave_r_factor)
    return

def generate_plot(self, app, plotQueues):

    '''
    objects for bokeh: 
    
    mvars.script_first_frame, mvars.div_first_frame : Plot of first frame
    mvars.script_first_frame_mc, mvars.div_first_frame_mc : Optional Plot of first frame

    If multi-frame case:

    mvars.script_last_frame, mvars.div_last_frame : Plot of first frame
    mvars.script_last_frame_mc, mvars.div_last_frame_mc : Optional Plot of first frame ( may not be necessary )

    mvars.script_average, mvars.div_average: Plot for averages over frame
    mvars.script_multi_frame, mvars.div_multi_frame: Plot for values as a function of frame
    '''

    pgui = self.run_utils.print_gui
    mvars = self.mvars

    # for frame = 0
    f = 0
    message = ""

    try: 
        mvars.script_first_frame, mvars.div_first_frame = plot.single_model_plot(self,f,400,500,app)
    except:
        message += "ALTENS failed to generate first frame plot."
        pgui(message)
        sys.exit()

    plotQueues['bokeh_plot_1'].put(mvars.script_first_frame)

#    print json.dumps(plotQueues['bokeh_plot_1'])

    if mvars.use_monte_carlo_flag:
#        try:
        mvars.script_first_frame_mc, mvars.div_first_frame_mc = plot.generate_mc_histogram(self,f,400,500,app)
#        except:
#            message += "ALTENS failed to generate MC plot."
#            pgui(message)
#            sys.exit()

        plotQueues['bokeh_plot_2'].put(mvars.script_first_frame_mc)

    if (mvars.number_of_frames > 1):
        f = mvars.number_of_frames - 1
        try:
            mvars.script_last_frame, mvars.div_last_frame = plot.single_model_plot(self,f,400,500,app)
        except:
            message += "ALTENS failed to generate last frame plot."
            pgui(message)
            sys.exit()
        plotQueues['bokeh_plot_3'].put(mvars.script_last_frame)

        try:        
            mvars.script_average, mvars.div_average = plot.single_model_plot(self,mvars.number_of_frames,400,500,app)
        except:
            message += "ALTENS failed to generate avereage plot."
            pgui(message)
            sys.exit()

        plotQueues['bokeh_plot_4'].put(mvars.script_average)
        #if mvars.use_monte_carlo_flag:
        #    mvars.script_last_frame_mc, mvars.div_last_frame_mc = plot.generate_mc_histogram(self,f,400,500,app)
        try:
            mvars.script_multi_frame, div_multi_frame = plot.multi_model_plot(self,400,500,app)
        except:
            message += "ALTENS failed to generate multi_frame plot."
            pgui(message)
            sys.exit()

        plotQueues['bokeh_plot_5'].put(mvars.script_multi_frame)

    time.sleep(1)

    return

def save_results(self,new_d, corr_coef, r_factor, s_tensor, s_eigenvalue, e_ang, frame, app):
 
    message = ""
    pgui = self.run_utils.print_gui
    mvars = self.mvars

    sfactor = 1000   

#   Write data
    output_dir = os.path.join(mvars.runname, app)
    output_rdc = mvars.runname+'_calc_rdc_'+str(frame+1).zfill(5)+'.txt'

    mvars.file_out[frame] = os.path.join(output_dir, mvars.runname+'_' + "results_" + str(frame+1).zfill(5)) + '.txt'
    mvars.file_calcrdc[frame] = os.path.join(output_dir, mvars.runname+'_calc_rdc_'+str(frame+1).zfill(5)) + '.txt'

    fout = open(mvars.file_out[frame],'w')
    frdc = open(mvars.file_calcrdc[frame],'w')  

#    if (mvars.number_of_frames > 1):
    if (frame == 0):
        mvars.file_out_summary = os.path.join(output_dir, mvars.runname+'_'+'results_per_model')+'.txt'
        fout_stack = open(mvars.file_out_summary, 'w')
        fout_stack.write("%12s %12s %12s %12s %12s %12s %12s %12s %12s \n" \
            %("Model", "Corr_coeff", "R-factor", "Axx[*10^-3]", "Ayy[*10^-3]", "Azz[*10^-3]", "Alpha", "Beta", "Gamma" )) 
    else:
        fout_stack = open(mvars.file_out_summary, 'a')

    fout.write("====================================================================\n")
    fout.write("======================== Altens Information ========================\n")
    fout.write("====================================================================\n")
    fout.write(" 1. PDB input file = %s\n"%mvars.pdbfile.split("/")[-1])
    fout.write(" 2. DCD input file = %s\n"%mvars.dcdfile.split("/")[-1])
    fout.write(" 3. RDC input file = %s\n"%mvars.rdc_input_file.split("/")[-1])
    if mvars.residue_list_file_flag:
        fout.write(" 4. Exclusion file = %s\n"%mvars.residue_list_file.split("/")[-1])
    else:
        fout.write(" 4. Exclusion file = None \n")
    fout.write(" 5. Monte Carlo analysis = %s\n"%mvars.use_monte_carlo_flag)
    fout.write(" 6. Back-calculated RDC file = %s\n"%output_rdc)
    fout.write(" 7. Correlation coefficient between RDC(calc) and RDC(exp) = %.5f Hz\n"%corr_coef[0][1])
    fout.write(" 8. R-factor for agreement between RDC(calc) and RDC(exp) = %.5f\n"%r_factor)
    fout.write(" 9. Order matrix [*10^-3] \n")
    fout.write("%12.5f %12.5f  %12.5f\n"%(s_tensor[0][0]*sfactor, s_tensor[0][1]*sfactor, s_tensor[0][2]*sfactor))
    fout.write("%12.5f %12.5f  %12.5f\n"%(s_tensor[1][0]*sfactor, s_tensor[1][1]*sfactor, s_tensor[1][2]*sfactor))
    fout.write("%12.5f %12.5f  %12.5f\n"%(s_tensor[2][0]*sfactor, s_tensor[2][1]*sfactor, s_tensor[2][2]*sfactor))

#   calculate alignment tensor (atensor)
    atensor_scale = 2.0/3.0
    atensor_xx = atensor_scale * s_eigenvalue[0]
    atensor_yy = atensor_scale * s_eigenvalue[1]
    atensor_zz = atensor_scale * s_eigenvalue[2]

    atensor_axial = 3.0/2.0*atensor_zz
    atensor_rhomb = atensor_xx - atensor_yy
    try:
        rhombicity = atensor_rhomb/atensor_axial
    except:
        message += "Poorly defined axial component of alignment tensor close to 0 was observed. Check the input RDC file or pdb file."  
        pgui(message)
        sys.exit()

#    fout.write("10. Alignment Tensor Eigenvalues [*10^-3]: Axx, Ayy, Azz = %.5f, %.5f, %.5f\n"%(s_eigenvalue[0]*sfactor,s_eigenvalue[1]*sfactor,s_eigenvalue[2]*sfactor))
    fout.write("10. Alignment tensor\n") 
    fout.write("    Eigenvalues [*10^-3]: Axx, Ayy, Azz = %.5f, %.5f, %.5f\n"%(atensor_xx*sfactor, atensor_yy*sfactor, atensor_zz*sfactor))
    fout.write("    Axial component (Aa) [*10^-3] = %.5f\n"%(atensor_axial*sfactor))
    fout.write("    Rhombic component (Ar) [*10^-3] = %.5f\n"%(atensor_rhomb*sfactor))
    fout.write("    Rhombicity (Ar/Aa)  = %.5f\n"%(rhombicity))
    fout.write("11. Euler angles [degrees] \n")
    fout.write("    Alpha = %.5f \n"%(e_ang[0]))
    fout.write("    Beta  = %.5f \n"%(e_ang[1]))
    fout.write("    Gamma = %.5f "%(e_ang[2]))
    fout.close()

#   Write output rdcs
    len_out_rdc = len(mvars.di)
    frdc.write("%12s %12s %12s %12s %12s %12s %12s %12s \n" %("ID", "Vectors", "type", "Chain", "Residue","RDC_exp", "RDC_calc", "delta(Exp-calc)"))
    for ii in xrange(len_out_rdc):
        atom_detail = str(mvars.rlist[ii][1])+str(mvars.rlist[ii][2])+str(mvars.rlist[ii][3])+str(mvars.rlist[ii][4])
        chain_id = str(mvars.rlist[ii][1])
        residue_id = str(mvars.rlist[ii][2])
        v_type = str(mvars.rlist[ii][3])+str(mvars.rlist[ii][4])
        frdc.write("%12s %12s %12s %12s %12s %12.5f %12.5f %12.5f \n"%(int(mvars.di[ii][0])+1, str(atom_detail), \
             str(v_type), str(chain_id), str(residue_id), mvars.di[ii][1]*mvars.dmax_array[ii], new_d[ii], \
             mvars.di[ii][1]*mvars.dmax_array[ii] - new_d[ii]))
    frdc.close()
#   Write summary 
    fout_stack.write("%12s %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f \n" \
        %(int(frame+1), corr_coef[0][1], r_factor, s_eigenvalue[0]*sfactor, s_eigenvalue[1]*sfactor, s_eigenvalue[2]*sfactor, \
        e_ang[0],e_ang[1],e_ang[2]))
    fout_stack.close()

    return

def save_average_results(self, app):

    message = ""
    mvars = self.mvars

    sfactor = 1000

    ''' write summary '''
#  frame = last frame    
#  frame + 1 ==> for average

    output_dir = os.path.join(mvars.runname, app)
    output_rdc = mvars.runname+'_calc_rdc_average'+'.txt'

    aveframe = mvars.number_of_frames
    mvars.file_out[aveframe] = os.path.join(output_dir, mvars.runname)+'_results_average.txt'
    mvars.file_calcrdc[aveframe] = os.path.join(output_dir, mvars.runname)+'_calc_rdc_average.txt'

    fout = open(mvars.file_out[aveframe], 'w')
    frdc = open(mvars.file_calcrdc[aveframe], 'w')

    fout.write("====================================================================\n")
    fout.write("=== Altens Average Information (except Monte Carlo Error Report) ===\n")
    fout.write("====================================================================\n")
    fout.write(" 1. PDB input file = %s\n"%mvars.pdbfile.split("/")[-1])
    fout.write(" 2. DCD input file = %s\n"%mvars.dcdfile.split("/")[-1])
    fout.write(" 3. RDC input file = %s\n"%mvars.rdc_input_file.split("/")[-1])
    if mvars.residue_list_file_flag:
        fout.write(" 4. Exclusion file = %s\n"%mvars.residue_list_file.split("/")[-1])
    else:
        fout.write(" 4. Exclusion file = None \n")
    fout.write(" 5. Monte Carlo analysis = %s\n"%mvars.use_monte_carlo_flag)
    fout.write(" 6. Calculated RDC file = %s\n"%output_rdc)
    fout.write(" 7. Correlation coefficient between RDC(exp) and averaged RDC(calc) = %.5f Hz\n"%mvars.ave_corr_coef[0][1])
    fout.write(" 8. R-factor for agreement between RDC(exp) and averaged RDC(calc) = %.5f\n"%mvars.ave_r_factor)
    fout.write(" 9. Order matrix [*10^-3] \n")
    fout.write("%12.5f(+-%.5f) %12.5f(+-%.5f) %12.5f(+-%.5f)\n"%(mvars.ave_s_tensor[0][0]*sfactor, mvars.std_s_tensor[0][0]*sfactor, mvars.ave_s_tensor[0][1]*sfactor, mvars.std_s_tensor[0][1]*sfactor, mvars.ave_s_tensor[0][2]*sfactor, mvars.std_s_tensor[0][2]*sfactor))
    fout.write("%12.5f(+-%.5f) %12.5f(+-%.5f) %12.5f(+-%.5f)\n"%(mvars.ave_s_tensor[1][0]*sfactor, mvars.std_s_tensor[1][0]*sfactor, mvars.ave_s_tensor[1][1]*sfactor, mvars.std_s_tensor[1][1]*sfactor, mvars.ave_s_tensor[1][2]*sfactor, mvars.std_s_tensor[1][2]*sfactor))
    fout.write("%12.5f(+-%.5f) %12.5f(+-%.5f) %12.5f(+-%.5f)\n"%(mvars.ave_s_tensor[2][0]*sfactor, mvars.std_s_tensor[2][0]*sfactor, mvars.ave_s_tensor[2][1]*sfactor, mvars.std_s_tensor[2][1]*sfactor, mvars.ave_s_tensor[2][2]*sfactor, mvars.std_s_tensor[2][2]*sfactor))

#   calculate alignment tensor (atensor)
    atensor_scale = 2.0/3.0
    atensor_xx = mvars.ave_s_eigenvalue[0] * atensor_scale
    atensor_yy = mvars.ave_s_eigenvalue[1] * atensor_scale
    atensor_zz = mvars.ave_s_eigenvalue[2] * atensor_scale

    atensor_axial = 3.0/2.0*atensor_zz
    atensor_rhomb = atensor_xx - atensor_yy

    std_atensor_xx = mvars.std_s_eigenvalue[0] * atensor_scale
    std_atensor_yy = mvars.std_s_eigenvalue[1] * atensor_scale
    std_atensor_zz = mvars.std_s_eigenvalue[2] * atensor_scale

    std_atensor_axial = 3.0/2.0*std_atensor_zz
    std_atensor_rhomb = math.sqrt(std_atensor_xx * std_atensor_xx + std_atensor_yy * std_atensor_yy )
    
    try:
        rhombicity = atensor_rhomb/atensor_axial
        std_rhombicity = abs(rhombicity) * math.sqrt((std_atensor_rhomb/atensor_rhomb)**2 + (std_atensor_axial/atensor_axial)**2 )
    except:
        message += "Poorly defined axial component of alignment tensor close to 0 was observed. Check the input RDC file or pdb file."
        pgui(message)
        sys.exit()

#    fout.write("10. Alignment Tensor Eigenvalues [*10^-3]: Axx, Ayy, Azz = %.5f, %.5f, %.5f\n"%(s_eigenvalue[0]*sfactor,s_eigenvalue[1]*sfactor,s_eigenvalue[2]*sfactor))
    fout.write("10. Alignment tensor\n")
    fout.write("    Eigenvalues [*10^-3]: Axx, Ayy, Azz = %.5f(+-%.5f), %.5f(+-%.5f), %.5f(+-%.5f)\n"%(atensor_xx*sfactor, std_atensor_xx*sfactor, atensor_yy*sfactor, std_atensor_yy*sfactor, atensor_zz*sfactor, std_atensor_zz*sfactor))
    fout.write("    Axial component (Aa) [*10^-3] = %.5f(+-%.5f)\n"%(atensor_axial*sfactor, std_atensor_axial*sfactor))
    fout.write("    Rhombic component (Ar) [*10^-3] = %.5f(+-%.5f)\n"%(atensor_rhomb*sfactor, std_atensor_rhomb*sfactor))
    fout.write("    Rhombicity (Ar/Aa)  = %.5f(+-%.5f)\n"%(rhombicity, std_rhombicity))
    fout.write("11. Euler angles [degrees]\n") 
    fout.write("    Alpha = %.5f (+-%.5f)\n"%(mvars.ave_e_ang[0], mvars.std_e_ang[0]))
    fout.write("    Beta  = %.5f (+-%.5f)\n"%(mvars.ave_e_ang[1], mvars.std_e_ang[1]))
    fout.write("    Gamma = %.5f (+-%.5f)"%(mvars.ave_e_ang[2], mvars.std_e_ang[2]))


    fout.close()

    ''' save output rdc files '''
    len_out_rdc = len(mvars.di)
    frdc.write("%12s %12s %12s %12s %12s %12s %12s %18s %18s \n" %("ID", "Vectors", "type", "Chain", "Residue","RDC_exp", "RDC_calc", "delta(Exp-calc)", "std_RDC_calc"))
    for ii in xrange(len_out_rdc):
        atom_detail = str(mvars.rlist[ii][1])+str(mvars.rlist[ii][2])+str(mvars.rlist[ii][3])+str(mvars.rlist[ii][4])
        chain_id = str(mvars.rlist[ii][1])
        residue_id = str(mvars.rlist[ii][2])
        v_type = str(mvars.rlist[ii][3])+str(mvars.rlist[ii][4])
        frdc.write("%12s %12s %12s %12s %12s %12.5f %12.5f %18.5f %18.5f \n"%(int(mvars.di[ii][0])+1, str(atom_detail), str(v_type), str(chain_id), str(residue_id), mvars.di[ii][1]*mvars.dmax_array[ii], mvars.ave_new_d[ii],mvars.di[ii][1]*mvars.dmax_array[ii]-mvars.ave_new_d[ii], mvars.std_new_d[ii]))
    frdc.close()

#  Write a residue list file for comparision
    sample_res_list = os.path.join(output_dir, mvars.runname)+'_residue_list.txt'
    sample_excl_list = os.path.join(output_dir, mvars.runname)+'_exclusion_list.txt'
    fres = open(sample_res_list,'w')
    fexcl = open(sample_excl_list,'w')

    count_type = 0
    prev_type = ""
    for i in xrange(len(mvars.rlist)):
        if (prev_type != mvars.rlist[i][0]):
            fres.write("# " + str(mvars.rlist[i][1]) + " " + str(mvars.rlist[i][3]) + " " + str(mvars.rlist[i][4]) + "\n" )
            fexcl.write("# " + str(mvars.rlist[i][1]) + " " + str(mvars.rlist[i][3]) + " " + str(mvars.rlist[i][4]) +"\n" )              
            prev_type = mvars.rlist[i][0]
            fres.write(str(int(mvars.rlist[i][2])) + "\n") 
        else:
            fres.write(str(int(mvars.rlist[i][2])) + "\n")
    fres.close()
    fexcl.close()

    return

def altens_core(self,app, frame, plotQueues):

    mvars = self.mvars
    frame = frame

    self.log.debug('in altens_core')
    pgui = self.run_utils.print_gui

    default_rdc_error = 1.0
    read_in_data(self,frame,default_rdc_error)
 
    dir_cos(self)

#    pgui("Calculate direction cosine Done\n")

    gen_a(self)

#    pgui("Generate matrix A Done\n")

    eigenvect, s_eigenvalue, s_tensor = x_tensor3(mvars.mat_a_inv,mvars.di)

    eigenvect_transpose = numpy.transpose(eigenvect)
    if eigenvect_transpose[0][0] > 0. :
        eigenvect_transpose = -1.0 * eigenvect_transpose

    e_ang = euler_angle(eigenvect_transpose)

    all_eulers, select_eulers = euler2all(e_ang[0],e_ang[1],e_ang[2])

    e_ang = select_eulers
    prin_angle = e_ang
    prin_param = s_eigenvalue

#    pgui("Single value decomposition Done\n")

    ''' extract Sij and recalculate Dij '''

    s_elements = extract_par(s_tensor)

    new_d = recalc_d(s_elements,mvars.mat_a)

#    ds = numpy.array([mvars.rlist,mvars.di,new_d,mvars.di-new_d])

    ''' calculate R factor '''

    corr_coef = numpy.corrcoef(mvars.di[:,1],new_d)

    r_factor = calc_r(mvars.di[:,1],new_d)

    an_param = calc_ani(prin_param)

#   Converting calculated di = new_d * dmax
    new_d = numpy.multiply(new_d,mvars.dmax_array)

#    if mvars.use_monte_carlo_flag:
#        mc_test(self,app,prin_param,prin_angle,frame)

    save_results(self,new_d, corr_coef, r_factor, s_tensor, s_eigenvalue, e_ang, frame, app)

    if mvars.use_monte_carlo_flag:
        mc_test(self,app,prin_param,prin_angle,frame)

#   Calculate averages 
    if (frame == 0 ):
        initialize_average(self, new_d, s_tensor, s_eigenvalue, e_ang, corr_coef, r_factor, an_param)    

#   Update sum 
    calc_average(self, new_d, s_tensor, s_eigenvalue, e_ang, corr_coef, r_factor, an_param, frame)

#   Save average values
    if (frame == mvars.number_of_frames -1 ):   
        save_average_results(self,app) 
        generate_plot(self, app, plotQueues)
     
#    pgui(60*"="+"\n")           
    return


def mc_test(self,app,prin_param,prin_angle,frame):

    # Monte Carlo error test
    mvars = self.mvars
    nsteps = mvars.number_of_monte_carlo_steps 
    e_set = numpy.ones((nsteps,3))
    s_set = numpy.ones((nsteps,3))
    n_gen = 0
    k = 0

    sfactor = 1000  # scaling factor for principal values (for printing)

    ''' Consider multi-frame work'''

    ''' setup trajectory file'''
    output_dir = os.path.join(mvars.runname, app)

#    mvars.file_mc_out[frame] = os.path.join(output_dir, mvars.runname+'_'+'mc_analysis_'+str(frame+1).zfill(5))+'.txt'
    mvars.file_mc_trajectory[frame] = os.path.join(output_dir, mvars.runname+'_'+'mc_'+str(frame+1).zfill(5))+'.txt'

#    if (mvars.number_of_frames > 1):
    if (frame == 0):
        mvars.file_mc_out_summary = os.path.join(output_dir, mvars.runname+'_'+'mc_per_model')+'.txt'
        fout_stack = open(mvars.file_mc_out_summary, 'w')
        fout_stack.write("%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s \n" \
            %("Model", "Axx[*10^-3]", "Std_Axx[*10^-3]", "Ayy[*10^-3]", "Std_Ayy[*10^-3]", "Azz[*10^-3]", "Std_Azz[*10^-3]", "Alpha", "Std_Alpha", \
            "Beta", "Std_Beta", "Gamma", "Std_Gamma", "Asym", "Std_Asym"))

    else:
        fout_stack = open(mvars.file_mc_out_summary, 'a')

    fout = open(mvars.file_out[frame],'a')
#    fout = open(mvars.file_mc_out[frame], 'w')
    ftr  = open(mvars.file_mc_trajectory[frame], 'w')

#    ftr = open(os.path.join(output_dir, mvars.runname+'_'+'mc_trajectory_'+str(frame+1).zfill(5))+'.txt', 'w')
#    fout = open(os.path.join(output_dir, mvars.runname+'_'+'mc_analysis_'+str(frame+1).zfill(5))+'.txt', 'w')

    ftr.write("%18s  %18s %18s  %18s  %18s  %18s\n" \
             %("Axx[*10^-3]", "Ayy[*10^-3]", "Azz[*10^-3]", "Alpha[degrees]", "Beta[degrees]", "Gamma[degrees]"))

    atensor_scale = 2.0/3.0
    for ii in range(0,nsteps):
        n_gen = n_gen + 1
        d_old = mvars.di
        d_new_mc = mc_table(self,d_old)
        eigenvect_mc, s_mc, s_tensor_mc = x_tensor3(mvars.mat_a_inv,d_new_mc)

        eigenvect_mc_transpose = numpy.transpose(eigenvect_mc)
        if eigenvect_mc_transpose[0][0] > 0:
            eigenvect_mc_transpose = -1.0 * eigenvect_mc_transpose

        e_ang = euler_angle(eigenvect_mc_transpose)
        all_eulers, select_eulers = euler2all(e_ang[0],e_ang[1],e_ang[2])
        e_ang = select_eulers

        if ii == 0:
           e_ang_full = numpy.array([e_ang])
           tab_param = numpy.array([[s_mc[0],s_mc[1],s_mc[2]]])
           tab_ang  = numpy.array([[e_ang[0],e_ang[1],e_ang[2]]])
        else:
           e_ang_full  = numpy.append(e_ang_full,numpy.array([e_ang]),axis = 0)

#    ''' check if S principal parameters are meaningful '''

        if ((s_mc[0]<=(prin_param[0]+5.0))&(s_mc[0]>=(prin_param[0]-5.0))&(s_mc[1]<=(prin_param[1]+5.0))&(s_mc[1]>=(prin_param[1]-5.0))& \
           (s_mc[2]<=(prin_param[2]+5.0))&(s_mc[2]>=(prin_param[2]-5.0))&(e_ang[0]<=(prin_angle[0]+60.0))&(e_ang[0]>=(prin_angle[0]-60.0)) \
           &(e_ang[1]<=(prin_angle[1]+60.0))&(e_ang[1]>=(prin_angle[1]-60.0))&(e_ang[2]<=(prin_angle[2]+60.0))&(e_ang[2]>=(prin_angle[2]-60.0))):

           k=k+1

#           if ii == 0:
#              tab_param = numpy.array([[s_mc[0],s_mc[1],s_mc[2]]]) 
#              tab_ang  = numpy.array([[e_ang[0],e_ang[1],e_ang[2]]])
           if (ii != 0):
              tab_param = numpy.append(tab_param,numpy.array([[s_mc[0],s_mc[1],s_mc[2]]]),axis=0)
              tab_ang = numpy.append(tab_ang,numpy.array([[e_ang[0],e_ang[1],e_ang[2]]]),axis=0)

        atensor_xx = atensor_scale * s_mc[0]
        atensor_yy = atensor_scale * s_mc[1]
        atensor_zz = atensor_scale * s_mc[2]

#        ftr.write("%12.5f  %12.5f %12.5f  %12.5f  %12.5f  %12.5f\n" \
#             %(s_mc[0]*sfactor,s_mc[1]*sfactor,s_mc[2]*sfactor,e_ang[0],e_ang[1],e_ang[2])) 
        ftr.write("%18.5f  %18.5f %18.5f  %18.5f  %18.5f  %18.5f\n" \
             %(atensor_xx*sfactor,atensor_yy*sfactor,atensor_zz*sfactor,e_ang[0],e_ang[1],e_ang[2]))

    ftr.close()
 
    sxx = numpy.mean(tab_param,axis=0)[0] * atensor_scale
    std_sxx = numpy.std(tab_param,axis=0)[0] * atensor_scale

    syy = numpy.mean(tab_param,axis=0)[1] * atensor_scale
    std_syy = numpy.std(tab_param,axis=0)[1] * atensor_scale

    szz = numpy.mean(tab_param,axis=0)[2] * atensor_scale
    std_szz = numpy.std(tab_param,axis=0)[2] * atensor_scale

    res_alpha = numpy.mean(tab_ang,axis=0)[0] 
    sd_res_alpha = numpy.std(tab_ang,axis=0)[0] 

    res_beta = numpy.mean(tab_ang,axis=0)[1]
    sd_res_beta = numpy.std(tab_ang,axis=0)[1] 

    res_gamma = numpy.mean(tab_ang,axis=0)[2] 
    sd_res_gamma = numpy.std(tab_ang,axis=0)[2] 

    param_mc = [sxx,syy,szz]
    an_param_mc = calc_ani(param_mc)
    error_ani_mc = an_param_mc*math.sqrt(((std_sxx*std_sxx + std_syy*std_syy)/(sxx-syy)/(sxx-syy)) + \
                       (std_szz*std_szz/szz/szz))

    fout.write("\n====================================================================\n")
    fout.write("==================== Monte Carlo Error Report ======================\n")
    fout.write("====================================================================\n")

    fout.write(" 1. Number of generated parameters = %s\n"%n_gen)
    fout.write(" 2. Axx[*10^-3] = %.5f +- %.5f\n"%(sxx*sfactor,std_sxx*sfactor))
    fout.write(" 3. Ayy[*10^-3] = %.5f +- %.5f\n"%(syy*sfactor,std_syy*sfactor))
    fout.write(" 4. Azz[*10^-3] = %.5f +- %.5f\n"%(szz*sfactor,std_szz*sfactor))

    aa_mc = szz*3.0/2.0
    std_aa_mc = std_szz*3.0/2.0

    ar_mc = sxx - syy
    std_ar_mc = math.sqrt(std_sxx**2 + std_szz**2) 

    rhomb_mc = ar_mc/aa_mc
    std_rhomb_mc = rhomb_mc * (math.sqrt((std_aa_mc/aa_mc)**2+(std_ar_mc/ar_mc)**2))

    fout.write(" 5. Aa[*10^-3] = %.5f +- %.5f\n"%(aa_mc*sfactor,std_aa_mc*sfactor))
    fout.write(" 6. Ar[*10^-3] = %.5f +- %.5f\n"%(ar_mc*sfactor,std_ar_mc*sfactor))
    fout.write(" 7. Rhombicity (Ar/Aa) = %.5f +- %.5f\n"%(rhomb_mc,std_rhomb_mc))

    fout.write(" 8. Alpha [degrees] = %.5f +- %.5f\n"%(res_alpha,sd_res_alpha))
    fout.write(" 9. Beta  [degrees] = %.5f +- %.5f\n"%(res_beta,sd_res_beta))
    fout.write("10. Gamma [degrees] = %.5f +- %.5f\n"%(res_gamma,sd_res_gamma))

#    fout.write("8. Asymmetry parameter = %.5f +- %.5f\n"%(an_param_mc,error_ani_mc))
    
    fout.close()

    fout_stack.write("%15s %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f\n" \
        %(int(frame+1), sxx*sfactor, std_sxx*sfactor, syy*sfactor, std_syy*sfactor, szz*sfactor, std_szz*sfactor, res_alpha, sd_res_alpha, \
        res_beta, sd_res_beta, res_gamma, sd_res_gamma, an_param_mc, error_ani_mc))
    fout_stack.close()
    
    return
