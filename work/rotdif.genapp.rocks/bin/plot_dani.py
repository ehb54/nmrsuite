#!/usr/bin/python
##order in dynamics output is s2, tauloc,rex,s2_fast
from preprocess_split import split_data

def dani_plot(color_list, exp_keys, in_dock, dyna_flag, elm_flag, elmdock_flag):
    dani_output = split_data(dyna_flag, elm_flag, elmdock_flag)[5]
    dani_data = []
    if in_dock == True:
        for i in range(len(exp_keys)):
            dani_tmp = [
       {
         "x" : dani_output[0][i],
         "y" : dani_output[1][i],
         "mode" : "lines+markers",
         "marker" : {
           "color" : color_list[i],
           "size": 8
         },
         "line" : {
         "color" : color_list[i],
         "width" : 1
         },
         "name": str(exp_keys[i])
       },
       {
         "x" : dani_output[0][i],
         "y" : dani_output[2][i],
         "type" : "bar",
         "marker": {"color": color_list[i]},
         "showlegend": False,
         "xaxis": "x2",
         "yaxis": "y2"
       },
       {
         "x" : dani_output[0][i],
         "y" : dani_output[3][i],
         "type" : "bar",
         "marker": {"color": color_list[i]},
         "showlegend": False,
         "xaxis" : "x3",
         "yaxis" : "y3",
         "type"  : "bar"
       },
       {
         "x": dani_output[0][i],
         "y": dani_output[4][i],
         "xaxis": "x4",
         "yaxis": "y4",
         "mode": "lines+markers",
         "showlegend": False,
         "marker": {
           "color" : color_list[i],
           "size" : 8
         },
         "line": {
           "color": color_list[i],
           "width":1
         }
       }    
       ]
            dani_data.extend(dani_tmp)

    else:
        dani_data = [
       {
         "x" : dani_output[0],
         "y" : dani_output[1],
         "mode" : "lines+markers",
         "marker" : {
           "color" : color_list[5],
           "size": 8
         },  
         "line" : { 
         "color" : color_list[5],
         "width" : 1 
         },  
         "name": "All Frequencies"
       },  
       {   
         "x" : dani_output[0],
         "y" : dani_output[2],
         "type" : "bar",
         "marker": {"color": color_list[5]},
         "showlegend": False,
         "xaxis": "x2",
         "yaxis": "y2"
       },  
       {   
         "x" : dani_output[0],
         "y" : dani_output[3],
         "type" : "bar",
         "marker": {"color": color_list[5]},
         "showlegend": False,
         "xaxis" : "x3",
         "yaxis" : "y3",
         "type"  : "bar"
       },  
       {   
         "x": dani_output[0],
         "y": dani_output[4],
         "xaxis": "x4",
         "yaxis": "y4",
         "mode": "lines+markers",
         "showlegend": False,
         "marker": {
           "color" : color_list[5],
           "size" : 8
         }
        }]    

    dani_plotly ={
   "data" : dani_data,
   "layout" : {
       "title" : "Dynamics: Fully Anisotropic Model",
       "xaxis" : {
         "title" : "residue number",
         "domain" : [0,0.45]
       },
       "yaxis" : {
         "title" : "S^2",
         "domain" : [0.65,1]
       },
       "xaxis2" : { 
         "title" : "residue number",
         "domain": [0.55,1],
         "anchor": "y2"
       },  
       "yaxis2" : { 
         "title" : "tau_loc(ns)",
         "domain": [0,0.35],
         "anchor": "x2"
       },
       "xaxis3" : {
         "title": "residue number",
         "domain": [0,0.45],
         "anchor": "y3"   
       },
       "yaxis3" : {
         "title": "Rex(1/s)",
         "domain": [0,0.35],
         "anchor": "x3"   
       },
       "xaxis4" : {
         "title": "residue number",
         "domain": [0.55,1],
         "anchor": "y4"
       },
       "yaxis4" : {
         "title": "S^2_fast",
         "domain": [0.65,1],
         "anchor": "x4"
       }
   }
}
    return dani_plotly

