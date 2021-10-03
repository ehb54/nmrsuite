#!/usr/bin/python
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
#from math import exp, expm1, log10, log, log1p
from subprocess import Popen, PIPE, STDOUT

def FileNameArr():
    ExpArr       = []
    FitArr_upper = []
    FitArr_lower = []
    LabelArr     = []
    Surface3d    = []
    Exp3d        = []
    rotout       = []
    for root, dirs, files in os.walk('./'):
        for name in files:
            filename = os.path.join(root, name)
            
            if '2d_fit_upper' in filename:
                FitArr_upper.append(filename)

            if '2d_fit_lower' in filename:
                FitArr_lower.append(filename)
                
            if '2d_exp_' in filename:
                ExpArr.append(filename)
                
            if '2d_label_' in filename:
                LabelArr.append(filename) 
                
            if '3d_surface_' in filename:
                Surface3d.append(filename)
            
            if '3d_exp_points_' in filename:
                Exp3d.append(filename)
            if 'rotdif' in filename:
                rotout.append(filename)
    ExpArr.sort()
    FitArr_upper.sort()
    FitArr_lower.sort()
    LabelArr.sort()
    Surface3d.sort()
    Exp3d.sort()

    return ExpArr, FitArr_upper, FitArr_lower, LabelArr, Surface3d, Exp3d

def GetData(filename):
    file = open(filename, 'rU')
    data = file.readlines()
    file.close()
        
    xData = []
    yData = []
    eData = []
    tData = []
    
    for line in data:
        xData.append(re.split(r'[\s]*', line)[0])
        yData.append(re.split(r'[\s]*', line)[1])
        if '2d_exp_' in filename:
            eData.append(re.split(r'[\s]*', line)[2])
            
            tooltipArr = re.split(r'[\s]*', line)[-3:]
            tooltipStr = ""
            for ch in range(len(tooltipArr)):
                tooltipStr += tooltipArr[ch] + " "  
            tooltipStr.rstrip() 
            tData.append(tooltipStr)

    return xData, yData, eData, tData

def GetData_3d(filename):
    file = open(filename, 'rU')
    data = file.readlines()
    file.close()
        
    xData = []
    yData = []
    zData = []
       
    for line in data:
        xData.append(re.split(r'[\s]*', line)[0])
        yData.append(re.split(r'[\s]*', line)[1])
        zData.append(re.split(r'[\s]*', line)[2])

    return xData, yData, zData

################################################################################
def GetLabels(filename):
    file = open(filename, 'rU')
    data = file.readlines()
    file.close()
        
    Labels = []
    
    for line in data:
        Labels.append(line)
       
    return Labels

