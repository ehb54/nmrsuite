#!/opt/miniconda2/bin/python
import pandas as pd
import numpy as np
from preprocess_dyna_rho import preprocess_out


def split_data(dyna_flag, elm_flag, elmdock_flag):
    with open('rotdif_results.out') as in_f:
        all_txt = in_f.read().splitlines()
        #get the lines
    dock_start = 0
    dock_end = 0
    diso_start = 0
    diso_end = 0
    daxi_start = 0
    daxi_end = 0
    dani_start = 0
    dani_end = 0
    elm_start = 0
    elm_end = 0
    for i in range(len(all_txt)):
        if 'Full Rotdif Results: Isotropic Solution' in all_txt[i]:
            iso_start = i+2
        elif 'Full Rotdif Results: Axially-Symmetric Solution' in all_txt[i]:
            axi_start = i+2
            iso_end = i-1
        elif 'Full Rotdif Results: Anisotropic Solution' in all_txt[i]:
            ani_start = i+2
            axi_end = i-1
            ani_end = ani_start + (axi_end - axi_start)
        elif 'Full Dynamics Results: Isotropic Solution' in all_txt[i]:
            diso_start = i+2
        elif 'Full Dynamics Results: Axially-Symmetric Solution' in all_txt[i]:
            daxi_start = i+2
            diso_end = i-1
        elif 'Full Dynamics Results: Anisotropic Solution' in all_txt[i]:
            dani_start = i+2
            daxi_end = i-1
            #last two rows are comments
            dani_end = dani_start + (daxi_end - daxi_start)
        elif 'ELM Tensor Prediction' in all_txt[i]:
            elm_start = i
        elif 'Eigendecomposition' in all_txt[i]:
            elm_end = i + 4
        elif 'Docking Results' in all_txt[i]:
            dock_start = i
            dock_end = len(all_txt) - 2
        #starter = [iso_start,axi_start,ani_start,diso_start,daxi_start,dani_start, elm_start]
        #ender = [iso_end,axi_end,ani_end,diso_end,daxi_end,dani_end, elm_end]
    iso_output = preprocess_out(iso_start, iso_end, all_txt,"exp")
    axi_output = preprocess_out(axi_start, axi_end, all_txt,"exp")
    ani_output = preprocess_out(ani_start, ani_end, all_txt,"exp")
    diso_output = []
    daxi_output = []
    dani_output = []
    elm_output = []
    dock_output = []
    if dyna_flag == True:
        diso_output = preprocess_out(diso_start, diso_end, all_txt,"dyna")
        daxi_output = preprocess_out(daxi_start, daxi_end, all_txt,"dyna")
        dani_output = preprocess_out(dani_start, dani_end, all_txt,"dyna")
    if elm_flag == True:
        elm_output = all_txt[elm_start:elm_end]
    if elmdock_flag == True:
        dock_output = all_txt[dock_start:dock_end]
    return iso_output,axi_output,ani_output,diso_output,daxi_output,dani_output, elm_output, dock_output
