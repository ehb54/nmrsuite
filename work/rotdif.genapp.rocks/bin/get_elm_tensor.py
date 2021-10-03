import re
import operator
import numpy as np
from os.path import join
def elm_tensor(pdb):
    with open('ELM_prediction') as in_f:
        tmp = in_f.readlines()
    for i in range(len(tmp)):
        if tmp[i].startswith('====Diffusion Tensor Sorted Eigendecomposition===='):
            row_dx = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i+1])[:3]]
            row_dy = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i+2])[:3]]
            row_dz = [float(it) for it in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+",tmp[i+3])[:3]]
    dx = [row_dx[0],row_dy[0],row_dz[0]]
    dy = [row_dx[1],row_dy[1],row_dz[1]]
    dz = [row_dx[2],row_dy[2],row_dz[2]]
    #normalize the length of dx, dy, dz by ellipsoid semi axes
    with open("str_abc.txt") as in_f:
        tmp_sa = in_f.readline()
    semiaxe = [float(s) for s in tmp_sa.split(' ')]
    dx = dx / np.linalg.norm(dx) * semiaxe[0]
    dy = dy / np.linalg.norm(dy) * semiaxe[1]
    dz = dz / np.linalg.norm(dz) * semiaxe[2]
    
    with open(pdb) as in_f2:
        tmp2 = in_f2.readlines()
    ori_x = 0
    ori_y = 0
    ori_z = 0
    count = 0
    for j in range(len(tmp2)):
        if tmp2[j].startswith('ATOM'):
            count += 1
            tmp_str = ''.join(tmp2[j])
            tmp_str = " ".join(tmp_str.split())
            tmp_ls = tmp_str.split(" ")
            ori_x += float(tmp_ls[6])
            ori_y += float(tmp_ls[7])
            ori_z += float(tmp_ls[8])
    
    ori_x = ori_x / count
    ori_y = ori_y / count
    ori_z = ori_z / count
    ori_coor = [ori_x, ori_y, ori_z]
    dx_1 = list(map(operator.sub, ori_coor , dx ))
    dx_2 = list(map(operator.add, ori_coor , dx ))
    dy_1 = list(map(operator.sub, ori_coor , dy ))
    dy_2 = list(map(operator.add, ori_coor , dy ))
    dz_1 = list(map(operator.sub, ori_coor , dz ))
    dz_2 = list(map(operator.add, ori_coor , dz ))
    all_x = [str(i) for i in dx_1+dx_2]
    all_y = [str(i) for i in dy_1+dy_2]
    all_z = [str(i) for i in dz_1+dz_2]
    obj = "[ CYLINDER, %s, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,\
         CYLINDER, %s, 0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,\
        CYLINDER, %s,0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,]" % (','.join(all_x), ','.join(all_y), ','.join(all_z))
    
    with open('ELM_tensor_axes.py','w') as out_f:
        out_f.write("""# axes plot script
from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
#create the axes object, draw axes with cylinders colored red (Dx), green (Dy), blue (Dz) \n""" +
"obj = " + obj +
"""\n # then we load it into PyMOL \n
cmd.load_cgo(obj,'diffusion_tensor')""")   
