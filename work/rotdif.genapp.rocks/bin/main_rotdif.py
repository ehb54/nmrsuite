#!/usr/bin/python
from os.path import join
import os
import glob
import re
import shutil
import sys
import subprocess
import json
import cStringIO
import shlex 
import socket
import time
import pty
import select
import atexit
import signal
import functools
import itertools
import math
import pandas as pd
import numpy as np
from itertools import groupby
from matplotlib import cm, colors
from StringIO import StringIO
from subprocess import Popen, PIPE, STDOUT

#import our functions
from run_rotdif import rotdif
from plot_2d_3d import plot_2d_3d
from plot_exp import exp_plot 
from plot_vec import vec_plot
from plot_iso import iso_plot
from plot_axi import axi_plot
from plot_ani import ani_plot
from plot_diso import diso_plot
from plot_daxi import daxi_plot
from plot_dani import dani_plot
from plot_chi2 import chi2_plot
from save_elm import save_elm
from save_dock import save_dock
from indicator_chain import dock
from StringIO import StringIO
import shutil

if __name__ == "__main__":
    color_list=[ '#1f77b4','#ff7f0e','#2ca02c','#d62728', '#9467bd',\
'#8c564b', '#e377c2', '#7f7f7f', '#bcbd22','#17becf']
    json_variables = " " 
    InitialDirectoryStr = os.path.abspath(os.path.dirname(sys.argv[0]))
    argv_io_string = StringIO(sys.argv[1])
    json_variables = json.load(argv_io_string)
    relax_location = json_variables['relax_location']
    old_relax = os.path.basename(relax_location[0])
    to_dock = dock(old_relax)
    pdb_location = json_variables['pdb_location']
    old_pdb = os.path.basename(pdb_location[0])

    run_name = json_variables['run_name']
    sub_dir = run_name + "_ROTDIF"
    if os.path.isdir(sub_dir):
        pass
    else:
        os.mkdir(sub_dir)

    new_relax = join(sub_dir,old_relax)
    new_pdb = join(sub_dir,old_pdb)
    shutil.move(old_relax,new_relax)
    shutil.move(old_pdb, new_pdb)
    os.chdir(sub_dir)
    output_res, base_dir, dyna_flag, elm_flag, elmdock_flag = rotdif(to_dock,old_pdb,old_relax,sub_dir)
    exp_graph,exp_keys,exp_resi,dyna_resi = exp_plot(color_list, dyna_flag, elm_flag, elmdock_flag,old_relax)
    outputmd_arr_files = []
    log_out = [join(str(base_dir),join(sub_dir,'rotdif_log.out'))]
    ani_out = [join(str(base_dir),join(sub_dir, 'out_ani.py'))]
    axi_out = [join(str(base_dir),join(sub_dir, 'out_axial.py'))]

    if dyna_flag == True:
        output_res[ 'diso_plot' ] = diso_plot(color_list, exp_keys, to_dock, dyna_flag, elm_flag, elmdock_flag)
        output_res[ 'daxi_plot' ] = daxi_plot(color_list, exp_keys, to_dock, dyna_flag, elm_flag, elmdock_flag)
        output_res[ 'dani_plot' ] = dani_plot(color_list,exp_keys, to_dock, dyna_flag, elm_flag, elmdock_flag)
 
    if elm_flag == True:
        save_elm(dyna_flag, elm_flag, elmdock_flag)
        elm_out = [join(str(base_dir),join(sub_dir,'ELM_prediction'))]
        output_res[ 'elm_out' ] = elm_out    

    if elmdock_flag == True:
        save_dock(dyna_flag, elm_flag, elmdock_flag)
        dock_out = [join(str(base_dir),join(sub_dir, 'ELMDOCK'))]     
        output_res[ 'elmdock' ] = dock_out
        dock_pdb = [join(str(base_dir),join(sub_dir, 'out_dock.pdb'))]
        output_res[ 'pdb'] = dock_pdb
        output_res['outputpdb'] = dock_pdb

    #get residule numbers for the vector plot
    output_res[ 'outputrotdif' ] = log_out
    output_res[ 'ani_out' ] = ani_out
    output_res[ 'axi_out' ] = axi_out
    output_res[ 'exp_plot' ] = exp_graph
    output_res[ 'vec_plot' ] = vec_plot(exp_resi)
    output_res[ 'chi2_plot' ] = chi2_plot()
    output_res[ 'plot_2d' ] = plot_2d_3d(color_list, exp_keys)[0]
    output_res[ 'plot_3d' ]  = plot_2d_3d(color_list, exp_keys)[1]
    output_res[ 'iso_plot' ] = iso_plot(exp_resi, color_list,exp_keys, dyna_flag, elm_flag, elmdock_flag)
    output_res[ 'axi_plot' ] = axi_plot(exp_resi, color_list,exp_keys, dyna_flag, elm_flag, elmdock_flag)
    output_res[ 'ani_plot' ] = ani_plot(exp_resi, color_list,exp_keys, dyna_flag, elm_flag, elmdock_flag)
    print (json.dumps(output_res))
