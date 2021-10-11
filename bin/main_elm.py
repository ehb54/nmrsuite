#!/opt/miniconda2/bin/python
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
from math import exp, expm1, log10, log, log1p
from subprocess import Popen, PIPE, STDOUT
from os.path import join
import shutil
from get_elm_tensor import elm_tensor
from get_ellipsoid import run_elli

def message_box(text,icon):
    
    _message = {}
    _message['icon'] = icon
    _message['text'] = text
    
    UDP_IP = json_variables['_udphost']
    UDP_PORT = int(json_variables['_udpport'])
    sock = socket.socket(socket.AF_INET, # Internet
                         socket.SOCK_DGRAM) # UDP
    
    socket_dict={}
    socket_dict['_uuid'] = json_variables['_uuid']
    socket_dict['_message'] = _message

    doc_string = json.dumps(socket_dict)
    sock.sendto(doc_string,(UDP_IP,UDP_PORT))
    
    return

def rotdif(Pdbfilename,sub_dir):
    output_res = {}
    json_variables = " "
    InitialDirectoryStr = os.path.abspath(os.path.dirname(sys.argv[0]))
    argv_io_string = StringIO(sys.argv[1])
    json_variables = json.load(argv_io_string)
    #pdb_location   = json_variables['pdb_location']
    #PdbStr = pdb_location[0]
    #Pdbfilename = os.path.basename(PdbStr)
## UDP messaging ##################################################
    UDP_IP = json_variables['_udphost']
    UDP_PORT = int(json_variables['_udpport'])
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    socket_dict={}
    socket_dict['_uuid'] = json_variables['_uuid']
    base_dir = json_variables[ '_base_directory' ]
    #optimization_method_list = json_variables[ 'optimization_method' ]
    ExecutableFileLocationRotdif = '/opt/genapp/rotdif/bin/rotdif-1.1.jar'
    ProcessToCallRotdif = []
    ProcessToCallRotdif.append('/usr/bin/java')
    ProcessToCallRotdif.append('-jar')
    ProcessToCallRotdif.append(ExecutableFileLocationRotdif)
    ProcessToCallRotdif.append('-nogui')
    ProcessToCallRotdif.append('-elm')
    ProcessToCallRotdif.append('-temp')
    ProcessToCallRotdif.append(json_variables['temperature'])
    ProcessToCallRotdif.append('-out')
    ProcessToCallRotdif.append('out')
    #ProcessToCallRotdif.append('-axes')
    #ProcessToCallRotdif.append('-axesl')
    #ProcessToCallRotdif.append(json_variables['axeslength'])
    ProcessToCallRotdif.append('-pdb')
    ProcessToCallRotdif.append(Pdbfilename)
    ProcessToCallRotdif.append('-model')
    ProcessToCallRotdif.append(json_variables['model'])
    #if advanced button is checked:
    if 'adv_option' in json_variables:
        ProcessToCallRotdif.append('-sr')
        ProcessToCallRotdif.append(json_variables['hydro'])
        ProcessToCallRotdif.append('-wr')
        ProcessToCallRotdif.append(json_variables['water'])
    
## RUN ROTDIF MD  ##########################################################
    socket_dict["_textarea"] = 'Starting ELM...\n\n'
    if socket_dict:
        doc_string = json.dumps(socket_dict)
        sock.sendto(doc_string,(UDP_IP,UDP_PORT))
    master_rotdif, slave_rotdif = pty.openpty()
    ProcessRotdif = subprocess.Popen( ProcessToCallRotdif,stdout=slave_rotdif,stdin=PIPE,bufsize=0,close_fds=True)    
    stdout_rotdif = os.fdopen(master_rotdif, 'r', 0)
    rotdif_log = open('all_prediction','w')
    path_to_live_log = join(str(base_dir),join(str(sub_dir),'all_prediction'))
    error_string_md = ''
    timeout = 4 # seconds

    while True:
        ready, _, _ = select.select([master_rotdif], [], [], timeout)
        if ready:
            output = stdout_rotdif.readline()
        #output = stdout_rotdif.read(10)
            if not output:
                break
            if output:
                socket_dict["_textarea"] = output#.strip()
                print >> rotdif_log, output.rstrip()
                rotdif_log.flush()
                if "Percent" in output:
                    output_strip = output.strip()
                    OutArr = re.split(r'[\s]*', output_strip)
                    percent = OutArr[1][:-1]
                    socket_dict['progress_output'] = float(percent)/float(100)
                    socket_dict['progress_text'] = 'ELM calculation progress: ' + str(int (float( percent ))) + '%'
                if socket_dict:
                    doc_string = json.dumps(socket_dict)
                    sock.sendto(doc_string,(UDP_IP,UDP_PORT))
        elif ProcessRotdif.poll() is not None:  
            break
              
    ProcessRotdif.wait()
    os.close(slave_rotdif)
    os.close(master_rotdif)
    output_res[ 'progress_output' ] = str(1.0)
    output_res[ 'progress_text' ] =  'ELM calculation progress: ' + '100%'                         
    socket_dict["_textarea"] = '\nELM Calculations Completed...\n\n'

    if socket_dict:
        doc_string = json.dumps(socket_dict)
        sock.sendto(doc_string,(UDP_IP,UDP_PORT))
                                                    
    output_res[ 'progress_output' ] = str(1.0)
    output_res[ 'progress_text' ] =  'ELM calculation progress: ' + '100%' 
    return output_res, base_dir

def clean_up():
    for CleanUp in glob.glob('./2d_*'):                            
        os.remove(CleanUp)

    for CleanUp in glob.glob('./3d_*'):                            
        os.remove(CleanUp)

    for CleanUp in glob.glob('./*.txt'): 
        if not CleanUp.endswith(Relaxfilename):
            os.remove(CleanUp)
    for CleanUp in glob.glob('./*.out'): 
        os.remove(CleanUp)

    for CleanUp in glob.glob('./*.pdb'):
        if not CleanUp.endswith(Pdbfilename):
            os.remove(CleanUp)

if __name__ == "__main__":
    output_res = {}
    json_variables = " " 
    InitialDirectoryStr = os.path.abspath(os.path.dirname(sys.argv[0]))
    argv_io_string = StringIO(sys.argv[1])
    json_variables = json.load(argv_io_string)
    pdb_location   = json_variables['pdb_location']
    PdbStr = pdb_location[0]
    Pdbfilename = os.path.basename(PdbStr)
    run_name = json_variables['run_name']
    sub_dir = run_name + "_ELM"
    #output_res,base_dir = rotdif(sub_dir)
    if os.path.isdir(sub_dir):
        sys.stderr.write("Run name already existed!")
    else:    
        os.mkdir(sub_dir)
    new_pdb = join(sub_dir,Pdbfilename)
    shutil.move(Pdbfilename,new_pdb)
    os.chdir(sub_dir)
    output_res,base_dir = rotdif(Pdbfilename, sub_dir)
    #save elm
    #elm_out = [join(str(base_dir),join(sub_dir,'ELM_prediction'))]       
    #save params
    #copy template
    shutil.copy("/opt/genapp/rotdif/bin/ELM_ellipsoid.py",join(str(base_dir),join(sub_dir,"ELM_ellipsoid.py")))
    with open(join(str(base_dir),join(sub_dir,"all_prediction"))) as in_f:
        tmp = in_f.readlines()
    tmp[0]="Parameters\n"
    with open(join(str(base_dir),join(sub_dir,"ELM_params")),'w') as out_f:
        for items in tmp[:17]:
            out_f.write(items)
    with open(join(str(base_dir),join(sub_dir,"ELM_prediction")),'w') as out_f2:
        for item in tmp[17:]:
            out_f2.write(item)
    elm_out = [join(str(base_dir),join(sub_dir,'ELM_prediction'))]
    elm_params = [join(str(base_dir),join(sub_dir,"ELM_params"))]
    output_res[ 'elm_out' ] = elm_out
    output_res['elm_params'] = elm_params
    pdb = join(str(base_dir),join(sub_dir,Pdbfilename))
    elm_tensor(pdb)
    output_res['elm_tensor'] = [join(str(base_dir), join(sub_dir,'ELM_tensor_axes.py'))]
    run_elli(join(str(base_dir), sub_dir))
    output_res['elm_ellipsoid'] = [join(str(base_dir), join(sub_dir,'ELM_ellipsoid.py'))]
    print(json.dumps(output_res))
