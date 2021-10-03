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
from math import exp, expm1, log10, log, log1p
from subprocess import Popen, PIPE, STDOUT

##import our functions
from preprocess_2d_3d import GetLabels, GetData_3d, GetData, FileNameArr

def plot_2d_3d(color_list, exp_keys):
    output_plot_2d_fit_upper = []
    output_plot_2d_fit_lower = []
    output_plot_2d_exp       = []
    output_plot_2d_rexp      = []
    labels                   = []
    output_plot_3d_surface   = []
    output_plot_3d_expp      = []
    min_data_x = []
    min_data_y = []
    max_data_x = []
    max_data_y = []
    num_exp       = len(FileNameArr()[0])
    num_fit_upper = len(FileNameArr()[1])
    num_fit_lower = len(FileNameArr()[2])
    # resizing
    for i in range(num_exp):
        output_plot_2d_exp.append([])
        output_plot_2d_rexp.append([])
        output_plot_3d_surface.append([])
        output_plot_3d_expp.append([])
    for i in range(num_fit_upper):
        output_plot_2d_fit_upper.append([])
    for i in range(num_fit_lower):
        output_plot_2d_fit_lower.append([])

#Exp - providing numbers of Exp, Fit and Labels sets ARE the same ...
    for i in range(num_exp):
        xExpArr = GetData(FileNameArr()[0][i])[0]
        yExpArr = GetData(FileNameArr()[0][i])[1]
        eExpArr = GetData(FileNameArr()[0][i])[2]   
        rExpArr = GetData(FileNameArr()[0][i])[3]   
        xFupArr = GetData(FileNameArr()[1][i])[0]
        yFupArr = GetData(FileNameArr()[1][i])[1]
        eFupArr = GetData(FileNameArr()[1][i])[2]   

        xFloArr = GetData(FileNameArr()[2][i])[0]
        yFloArr = GetData(FileNameArr()[2][i])[1]
        eFloArr = GetData(FileNameArr()[2][i])[2] 

        labels.append( GetLabels(FileNameArr()[3][i]) ) 

        xExpArr_3d = GetData_3d(FileNameArr()[4][i])[0]
        yExpArr_3d = GetData_3d(FileNameArr()[4][i])[1]
        zExpArr_3d = GetData_3d(FileNameArr()[4][i])[2]   

        xExpArr_3d_expp = GetData_3d(FileNameArr()[5][i])[0]
        yExpArr_3d_expp = GetData_3d(FileNameArr()[5][i])[1]
        zExpArr_3d_expp = GetData_3d(FileNameArr()[5][i])[2]   
        #xFupArrsorted = xFupArr.sort()
        #yFupArrsorted = yFupArr.sort()    
        min_data_x.append( min(ii for ii in xFupArr ) )
        min_data_y.append( min(ii for ii in yExpArr ) )
        max_data_x.append( max(ii for ii in xFupArr ) )
        max_data_y.append( max(ii for ii in yExpArr ) )
        
        numpoints_exp = len(xExpArr)
        numpoints_fup = len(xFupArr)
        numpoints_flo = len(xFloArr)
            
        for j in range(numpoints_exp):
            output_plot_2d_exp[i].append([xExpArr[j],yExpArr[j],eExpArr[j]])  
            output_plot_2d_rexp[i].append([rExpArr[j]])  
        for j in range(numpoints_fup):
            output_plot_2d_fit_upper[i].append([xFupArr[j],yFupArr[j]])     
        for j in range(numpoints_fup):
            output_plot_2d_fit_lower[i].append([xFloArr[j],yFloArr[j]])  

        output_plot_3d_surface[i].append(xExpArr_3d)
        output_plot_3d_surface[i].append(yExpArr_3d)
        output_plot_3d_surface[i].append(zExpArr_3d)  
        output_plot_3d_expp[i].append(xExpArr_3d_expp)
        output_plot_3d_expp[i].append(yExpArr_3d_expp)
        output_plot_3d_expp[i].append(zExpArr_3d_expp)  
   
    min_x = min(min_data_x)
    min_y = min(min_data_y)
    max_x = max(max_data_x)
    max_y = max(max_data_y)
    min_y1 = float(min_y) - 0.2
    max_y1 = float(max_y) + 0.2

    xscale = ''
    yscale = ''
    data_2d = {
        "options": 
    {
        "title"  : "Experimental &#961;",
        "ymin"   : min_y1,
        "ymax"   : max_y1,
        "xlabel" : "&#952;",
        "ylabel" : "&#961;(Exp.)",
        "legend" :
        {            
            "position"            : "ne",
            "margin"              : [-140, -1],
            "backgroundColor"     : "null",
            "labelBoxBorderColor" : "#000000",
         },
        "grid" : 
        {
            "backgroundColor" : "#ffffff",
            "margin"    : 
                              { 
                                  "top"    : 0, 
                                  "left"   : 0,
                                  "bottom" : 0,
                                  "right"  : 140
                              }
        },
    },
    "data": [ ]
}


## 3D with PLOTLY.JS ############################################################
    data_3d = {
    "layout" : {
        "scene": {
            "xaxis" : { 
                "title"           : "Theta (Deg)",
                "backgroundcolor" : "rgb(200, 200, 230)",
                "gridcolor"       : "rgb(128, 128, 128)",
                "showbackground"  : "true",
                "zerolinecolor"   : "rgb(128, 128, 128)"
            },
            "yaxis" : { 
                "title" : "Phi (Deg)",
                "backgroundcolor" : "rgb(230, 200,230)",
                "gridcolor"       : "rgb(128, 128, 128)",
                "showbackground"  : "true",
                "zerolinecolor"   : "rgb(128, 128, 128)",
                "range"           : [-179, 180]
            },
            "zaxis" : { 
                "title" : "Rho",
                "backgroundcolor" : "rgb(230, 230,200)",
                "gridcolor"       : "rgb(128, 128, 128)",
                "showbackground"  : "true",
                "zerolinecolor"   : "rgb(128, 128, 128)"
            },
            "camera": {
                "center" : { "x": 0, "y": 0, "z": 0 }, 
                "eye"    : { "x": 2, "y": 2, "z": 1 }, 
                "up"     : { "x": 0, "y": 0, "z": 1 }
            }
        },
        "showlegend" : "false",
        "autosize" : "false",
	"width" : 550,
	"height": 550,
	"margin": {
	 "l": 0,
	 "r": 0,
	 "b": 0,
	 "t": 0,
	 "pad": 4
	},
    },
    "data"   : []
}

    for i in range(num_exp):
        current_set_surf = {
        "type"    :  "mesh3d",
        "opacity" :  0.6,
        "x"       :  output_plot_3d_surface[i][0],
        "y"       :  output_plot_3d_surface[i][1],
        "z"       :  output_plot_3d_surface[i][2],
        "color"   :  color_list[i],
        "showlegend" : False
    }
        data_3d['data'].append(current_set_surf)

    for i in range(num_exp):
   # print output_plot_3d_expp[i][0]
   # exit();
        current_set_3d_expp = {
        "type"    :  "scatter3d",
        "mode"    :  "markers",
        "showlegend" : False,
        "marker"  : {
            "color"   : color_list[i],
            "size"    : 2,
            "symbol"  : 'square',
            "line"    :  {
                "color"  : "rgb(0,0,0)",
		"width"  : 1.5 
            },
            "opacity" : 1.0
        },
        "opacity" :  1.0,
        "x"       :  output_plot_3d_expp[i][0],
        "y"       :  output_plot_3d_expp[i][1],
        "z"       :  output_plot_3d_expp[i][2],
        "color"   :  color_list[i],
        "name"    :  labels[i]
    }
        data_3d['data'].append(current_set_3d_expp)
## 3D with PLOTLY.JS ############################################################
    data_3d = {
    "layout" : {
        "scene": {
            "xaxis" : { 
                "title"           : "Theta (Deg)",
                "backgroundcolor" : "rgb(200, 200, 230)",
                "gridcolor"       : "rgb(128, 128, 128)",
                "showbackground"  : "true",
                "zerolinecolor"   : "rgb(128, 128, 128)"
            },
            "yaxis" : { 
                "title" : "Phi (Deg)",
                "backgroundcolor" : "rgb(230, 200,230)",
                "gridcolor"       : "rgb(128, 128, 128)",
                "showbackground"  : "true",
                "zerolinecolor"   : "rgb(128, 128, 128)",
                "range"           : [-179, 180]
            },
            "zaxis" : { 
                "title" : "Rho",
                "backgroundcolor" : "rgb(230, 230,200)",
                "gridcolor"       : "rgb(128, 128, 128)",
                "showbackground"  : "true",
                "zerolinecolor"   : "rgb(128, 128, 128)"
            },
            "camera": {
                "center" : { "x": 0, "y": 0, "z": 0 }, 
                "eye"    : { "x": 2, "y": 2, "z": 1 }, 
                "up"     : { "x": 0, "y": 0, "z": 1 }
            }
        },
        "showlegend" : "false",
        "autosize" : "false",
	"width" : 550,
	"height": 550,
	"margin": {
	 "l": 0,
	 "r": 0,
	 "b": 0,
	 "t": 0,
	 "pad": 4
	},
    },
    "data"   : []
}

    for i in range(num_exp):
        current_set_surf = {
        "type"    :  "mesh3d",
        "opacity" :  0.6,
        "x"       :  output_plot_3d_surface[i][0],
        "y"       :  output_plot_3d_surface[i][1],
        "z"       :  output_plot_3d_surface[i][2],
        "color"   :  color_list[i],
        "showlegend" : False
    }
        data_3d['data'].append(current_set_surf)

    for i in range(num_exp):
        current_set_3d_expp = {
        "type"    :  "scatter3d",
        "mode"    :  "markers",
        "showlegend" : False,
        "marker"  : {
            "color"   : color_list[i],
            "size"    : 2,
            "symbol"  : 'square',
            "line"    :  {
                "color"  : "rgb(0,0,0)",
		"width"  : 1.5 
            },
            "opacity" : 1.0
        },
        "opacity" :  1.0,
        "x"       :  output_plot_3d_expp[i][0],
        "y"       :  output_plot_3d_expp[i][1],
        "z"       :  output_plot_3d_expp[i][2],
        "color"   :  color_list[i],
        "name"    :  exp_keys[i]
    }
        data_3d['data'].append(current_set_3d_expp)
### 2D ###########################################################

    for i in range(num_exp):
    #data_label = "Data" + str(i+1)
        data_label = labels[i]
        current_set_data = {
        "points"  : { "show"      : "true", 
                      "radius"    : 2,
                      "errorbars" : "y", 
                      "yerr"      : { "lowerCap" : "-",
                                      "upperCap" : "-", 
                                      "show"     : "true", 
                                      #"color"    : "red",
                                      "radius"   : 3 
                                  }
                  },
        "label"  : data_label,
        "data"   : output_plot_2d_exp[i],
        "tooltips"  : output_plot_2d_rexp[i],
        "color"  : color_list[i]
    }
        data_2d['data'].append(current_set_data)


    for i in range(num_exp):
        fit_label = "Fit" + str(i+1)
        current_set_fit = {
        "lines"  : { "show" : "true"  },
        "data"   : output_plot_2d_fit_upper[i],
        "color"  : color_list[i]
    }
        data_2d['data'].append(current_set_fit)
  

    for i in range(num_exp):
        fit_label = "Fit" + str(i+1)
        current_set_fit = {
        "lines"  : { "show" : "true" },
        "data"   : output_plot_2d_fit_lower[i],
         "color"  : color_list[i]
    }
        data_2d['data'].append(current_set_fit)
    return data_2d, data_3d
