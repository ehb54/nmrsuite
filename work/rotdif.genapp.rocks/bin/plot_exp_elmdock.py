#!usr/bin/python
import re
import pandas as pd
import numpy
import os
import sys
#import our functions
from elmdock_split import split_data

def exp_plot(color_list, dyna_flag, elm_flag, elmdock_flag, relax_loc):
    iso_output = split_data(dyna_flag, elm_flag, elmdock_flag)[0]
    #for root, dirs, files in os.walk('./'):
    #    for name in files:
    #        filename = os.path.join(root, name)
    #        if '.txt' in filename:
    #            exp_file = filename
    exp_file = relax_loc
    ##preprocess exp_file in multi-frequencies
    with open(exp_file) as in_f:
        tmp_lists = in_f.read().splitlines()
        exp_list = []
        for items in tmp_lists:
            if len(items) > 0:
                exp_list.append(items.split()) 
    exp_pd = pd.DataFrame(data = exp_list[1:],\
columns = ["residue","chain","atom 1","atom 2","magnet","T1","T1 error","T2","T2 error","NOE","NOE error"])
    resi = list(exp_pd["residue"])
    digi_resi = []
    #eliminate * in residue number
    for list_elem in resi:
        digi_resi.extend(re.findall("\d+", list_elem)) 
    digi_resi = [int(elem) for elem in digi_resi]
    exp_pd["residue"] = digi_resi 
    sort_exp = exp_pd.sort_values(['magnet', 'residue'], ascending=[True, True])
    freq_set = set(list(exp_pd["magnet"]))
    chain_set = set(list(exp_pd["chain"]))
    exp_dict = {} 
    #for multi-field, it should be frequencies; for multi-chains, it should be chains.
    exp_keys = []
    #split data based on frequency difference
    if len(chain_set) == 0:
        sys.stderr.write("No chain information in the data!")
    elif len(chain_set) == 1:    
        for items in freq_set:
            #tmp = sort_exp[sort_exp["magnet"]==items]
            tmp = exp_pd[exp_pd["magnet"]==items]
            exp_dict[items] = tmp
            freq_tag = str(items) + "MHz"
            exp_keys.append(freq_tag)
#from previous block, iso, axi, ani all have the same exp data
        rho_exp = iso_output[0]
        rho_err = iso_output[2]
        rho_block = []
#rho data has the same order as rho exp
        for i in range(len(exp_keys)):
            rho_data ={
         "x" : list(range(len(rho_exp[i]))),
         "y" : rho_exp[i],
         "error_y": {
            "array" : rho_err[i],
            "visible" : True 
         }, 
         "xaxis": "x4",
         "yaxis": "y4",  
         "mode" : "lines+markers",
         "line": {"color": color_list[i]},
         "marker": {"color": color_list[i]},
         "name" : exp_keys[i],
         "legendgroup": str(i),
         "showlegend" : True
       }    
            rho_block.append(rho_data)
        T1_block = []
        T2_block = []
        NOE_block = []
        j = 0
        for items in freq_set:
            T1_list = [float(ele) for ele in list(exp_dict[items]["T1"])]
            T1_err_list = [float(ele) for ele in list(exp_dict[items]["T1 error"])]
            T2_list = [float(ele) for ele in list(exp_dict[items]["T2"])]
            T2_err_list = [float(ele) for ele in list(exp_dict[items]["T2 error"])]
            NOE_list = [float(ele) for ele in list(exp_dict[items]["NOE"])]
            NOE_err_list = [float(ele) for ele in list(exp_dict[items]["NOE error"])]
        #for legend
            T1_data ={
         "x" : list(exp_dict[items]["residue"]),
         "y" : T1_list,
         "error_y": {
            "array" : T1_err_list,
            "visible" : True
         },
         "xaxis": "x",
         "yaxis": "y",
         "mode" : "lines+markers",
         "line": {"color": color_list[j]},
         "marker": {"color": color_list[j]},
         "name" : items,
         "legendgroup" : str(j),
         "showlegend" : False
       }
            T1_block.append(T1_data)

            T2_data ={
         "x" : list(exp_dict[items]["residue"]),
         "y" : T2_list,
         "error_y": {
            "array" : T2_err_list,
            "visible" : True 
         },
         "xaxis": "x3",
         "yaxis": "y3",   
         "mode" : "lines+markers",
         "line": {"color": color_list[j]},
         "marker": {"color": color_list[j]},
         "name" : items,
         "legendgroup": str(j),
         "showlegend" : False
       }
            T2_block.append(T2_data)

            NOE_data ={
         "x" : list(exp_dict[items]["residue"]),
         "y" : NOE_list,
         "error_y": {
            "array" : NOE_err_list,
            "visible" : True
         },
         "xaxis": "x2",
         "yaxis": "y2",
         "mode" : "lines+markers",
         "line": {"color": color_list[j]},
         "marker": {"color": color_list[j]},
         "name" : items,
         "legendgroup" : str(j),
         "showlegend" : False
       }
            NOE_block.append(NOE_data)
            j += 1

#flat list of lists
        list_of_list = [rho_block,T1_block,T2_block,NOE_block]
        list_of_exp = [item for sublist in list_of_list for item in sublist]
        exp_plotly = {
   "data" : list_of_exp,
   "layout" : {
       "title" : "Experimental Data",
       "xaxis" : {
         "title" : "residue number",
         "domain" : [0,0.45]
       },
       "yaxis" : {
         "title" : "R1_exp",
         "domain" : [0.55,1]
       },
       "xaxis2" : {
         "title" : "residue number",
         "domain": [0.55,1],
         "anchor": "y2"
       },
       "yaxis2" : {
         "title" : "NOE_exp",
         "domain": [0,0.35],
         "anchor": "x2"
       },
       "xaxis3" : {
         "title": "residue number",
         "domain": [0,0.45],
         "anchor": "y3"
       },
       "yaxis3" : {
         "title": "R2_exp",
         "domain": [0,0.35],
         "anchor": "x3"
       },
       "xaxis4" : {
         "title": "residue number",
         "domain": [0.55,1],
         "anchor": "y4"
       },
       "yaxis4" : {
         "title": "Rho_exp",
         "domain": [0.65,1],
         "anchor": "x4"
       }
   }
}
    elif len(chain_set) == 2:
        for items in chain_set:
            chain_tag = "chain " + str(items)
            exp_keys.append(chain_tag)
            tmp = sort_exp[sort_exp["chain"]==items]
            exp_dict[items] = tmp
#from previous block, iso, axi, ani all have the same exp data
        rho_exp = iso_output[0]
        rho_err = iso_output[2]
        rho_block = []
#rho data has the same order as rho exp
        for i in range(len(exp_keys)):
            rho_data ={
         "x" : list(range(len(rho_exp[i]))),
         "y" : rho_exp[i],
         "error_y": {
            "array" : rho_err[i],
            "visible" : True 
         },   
         "xaxis":"x4",
         "yaxis":"y4",
         "mode" : "lines+markers",
         "line": {"color": color_list[i]},
         "marker": {"color": color_list[i]},
         "name" : exp_keys[i],
         "legendgroup": str(i),
         "showlegend" : True
       }    
            rho_block.append(rho_data)
        T1_block = []
        T2_block = []
        NOE_block = []
        j = 0 
        for items in chain_set:
            T1_list = [float(ele) for ele in list(exp_dict[items]["T1"])]
            T1_err_list = [float(ele) for ele in list(exp_dict[items]["T1 error"])]
            T2_list = [float(ele) for ele in list(exp_dict[items]["T2"])]
            T2_err_list = [float(ele) for ele in list(exp_dict[items]["T2 error"])]
            NOE_list = [float(ele) for ele in list(exp_dict[items]["NOE"])]
            NOE_err_list = [float(ele) for ele in list(exp_dict[items]["NOE error"])]
        #for legend
            T1_data ={
         "x" : list(exp_dict[items]["residue"]),
         "y" : T1_list,
         "error_y": {
            "array" : T1_err_list,
            "visible" : True
         },
         "xaxis": "x",
         "yaxis": "y",
         "mode" : "lines+markers",
         "line": {"color": color_list[j]},
         "marker": {"color": color_list[j]},
         "name" : items,
         "legendgroup" : str(j),
         "showlegend" : False
       }
            T1_block.append(T1_data)

            T2_data ={
         "x" : list(exp_dict[items]["residue"]),
         "y" : T2_list,
         "error_y": {
            "array" : T2_err_list,
            "visible" : True
         },
         "xaxis": "x3",
         "yaxis": "y3",
         "mode" : "lines+markers",
         "line": {"color": color_list[j]},
         "marker": {"color": color_list[j]},
         "name" : items,
         "legendgroup": str(j),
         "showlegend" : False
       }
            T2_block.append(T2_data)

            NOE_data ={
         "x" : list(exp_dict[items]["residue"]),
         "y" : NOE_list,
         "error_y": {
            "array" : NOE_err_list,
            "visible" : True
         },
         "xaxis": "x2",
         "yaxis": "y2",
         "mode" : "lines+markers",
         "line": {"color": color_list[j]},
         "marker": {"color": color_list[j]},
         "name" : items,
         "legendgroup" : str(j),
         "showlegend" : False
       }
            NOE_block.append(NOE_data)
            j += 1

#flat list of lists
        list_of_list = [rho_block,T1_block,T2_block,NOE_block]
        list_of_exp = [item for sublist in list_of_list for item in sublist]
        exp_plotly = {
   "data" : list_of_exp,
   "layout" : {
       "title" : "Experiment Data Plot",
       "xaxis" : {
         "title" : "residue number",
         "domain" : [0,0.45]
       },
       "yaxis" : {
         "title" : "R1_exp(1/s)",
         "domain" : [0.55,1]
       },
       "xaxis2" : {
         "title" : "residue number",
         "domain": [0.55,1],
         "anchor": "y2"
       },
       "yaxis2" : {
         "title" : "NOE_exp",
         "domain": [0,0.35],
         "anchor": "x2"
       },
       "xaxis3" : {
         "title": "residue number",
         "domain": [0,0.45],
         "anchor": "y3"
       },
       "yaxis3" : {
         "title": "R2_exp(1/s)",
         "domain": [0,0.35],
         "anchor": "x3"
       },
       "xaxis4" : {
         "title": "residue number",
         "domain": [0.55,1],
         "anchor": "y4"
       },
       "yaxis4" : {
         "title": "Rho_exp",
         "domain": [0.65,1],
         "anchor": "x4"
       }
   }
}
             
    else:
        sys.stderr.write("More than 2 chains are not supported")

    return exp_plotly,exp_keys
