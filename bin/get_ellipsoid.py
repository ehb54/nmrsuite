import re
import os, fnmatch
import math

def calc_mass(pdb_file):
    '''
    input: pdb file
    output: center of mass, 3 coordinates, unit: Angstrom
    '''
    with open(pdb_file) as in_f:
        tmp = in_f.readlines()
    ori_x = 0
    ori_y = 0
    ori_z = 0
    count = 0
    for j in range(len(tmp)):
        if tmp[j].startswith('ATOM'):
            count += 1
            if count == 1:
                 start_ind = j
                 start_tmp = tmp[j]
                 start_str = ''.join(start_tmp)
                 start_str = " ".join(start_str.split())
                 start_ls = start_str.split(" ")
                 start_atm = [float(start_ls[6]), float(start_ls[7]), float(start_ls[8])]
            tmp_str = ''.join(tmp[j])
            tmp_str = " ".join(tmp_str.split())
            tmp_ls = tmp_str.split(" ")
            ori_x += float(tmp_ls[6])
            ori_y += float(tmp_ls[7])
            ori_z += float(tmp_ls[8])
    #end_tmp = tmp[start_ind + count-1]
    #end_str = ''.join(end_tmp)
    #end_str = " ".join(end_str.split())
    #end_ls = end_str.split(" ")
    #end_atm = [float(end_ls[6]), float(end_ls[7]), float(end_ls[8])]
    #calculate the protein diameter to get the scale factor
    #diameter = 0
    #for i in range(3):
    #    diameter += (start_atm[i] - end_atm[i])**2
    #math.sqrt(diameter)
    ori_x = ori_x / count
    ori_y = ori_y / count
    ori_z = ori_z / count
    
    return ori_x, ori_y, ori_z

def angle_axe(elm_prediction):
    '''
    input: ELM prediction file
    output: 1) 3 Euler angles (unit: degree)
    '''
    with open(elm_prediction) as in_f:
        tmp = in_f.readlines()
    for i in range(len(tmp)):
        #if tmp[i].startswith('Dx'):
        #    dx = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i])][0]
        #if tmp[i].startswith('Dy'):
        #    dy = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i])][0]
        #if tmp[i].startswith('Dz'):
        #    dz = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i])][0]
        if tmp[i].startswith('alpha'):
            alpha = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i])][0]
        if tmp[i].startswith('beta'):
            beta = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i])][0]
        if tmp[i].startswith('gamma'):
            gamma = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i])][0]
    return alpha, beta, gamma

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def find_n(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def get_abc():
    with open("str_abc.txt") as in_f:
        strabc = in_f.readline()
    fabc = [float(s) for s in strabc.split(' ')]
    return fabc[0], fabc[1], fabc[2]        

def run_elli(base_dir):
    # Input:
    pdb_file = find_n("*.pdb",base_dir)[0]
    elm_prediction = find("ELM_prediction",base_dir) 
    # 1. Center of Mass of the molecule
    cmx, cmy, cmz = calc_mass(pdb_file)
    #2. Ellipsoid Semiaxes (in A)
    #Rotation Input: three euler angles: alpha, beta, gamma (in degrees)
    deg1,deg2,deg3 = angle_axe(elm_prediction)
    rotationInput = [deg1, deg2, deg3]
    color = [0.85, 0.85, 1.00]
    a1, a2, a3 = get_abc()
    to_append = """ 
#1. Center of Mass of the molecule
cmx, cmy, cmz = %0.3f, %0.3f, %0.3f
#2. Ellipsoid semiaxes length (in Angstrom)
a1,a2,a3 = %0.3f, %0.3f, %0.3f
#3. Color: see https://pymolwiki.org/index.php/Color_Values for more options
color = [0.85, 0.85, 1.00]
#4. Rotation Input: three Euler angles: alpha, beta, gamma (in degrees)
rotationInput = [%0.3f, %0.3f, %0.3f]
tmp = drawEllipsoid(color, cmx, cmy, cmz, a1, a2, a3, *rotationMatrix(rotationInput))
cmd.load_cgo(tmp, 'ellipsoid-cgo')
cmd.set('cgo_transparency', 0.5, 'ellipsoid-cgo')"""%(cmx, cmy, cmz, a1, a2, a3, deg1, deg2, deg3)
    with open("ELM_ellipsoid.py","a") as out_f:
        out_f.write(to_append)
