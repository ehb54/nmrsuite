#!/opt/miniconda2/bin/python
##order in dynamics output is s2, tauloc,rex,s2_fast
from preprocess_split import split_data

def daxi_plot(color_list, exp_keys, in_dock, dyna_flag, elm_flag, elmdock_flag):
    daxi_output = split_data(dyna_flag, elm_flag, elmdock_flag)[4]
    daxi_data = []
    if in_dock == True:
        for i in range(len(exp_keys)):
            daxi_tmp = [
       {
         "x" : daxi_output[0][i],
         "y" : daxi_output[1][i],
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
         "x" : daxi_output[0][i],
         "y" : daxi_output[2][i],
         "type" : "bar",
         "marker": {"color": color_list[i]},
         "showlegend": False,
         "xaxis": "x2",
         "yaxis": "y2"
       },
       {
         "x" : daxi_output[0][i],
         "y" : daxi_output[3][i],
         "type" : "bar",
         "marker": {"color": color_list[i]},
         "showlegend": False,
         "xaxis" : "x3",
         "yaxis" : "y3",
         "type"  : "bar"
       },
       {
         "x": daxi_output[0][i],
         "y": daxi_output[4][i],
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
            daxi_data.extend(daxi_tmp)

    else:
        daxi_data = [
       {
         "x" : daxi_output[0],
         "y" : daxi_output[1],
         "mode" : "lines+markers",
         "marker" : {
           "color" : color_list[4],
           "size": 8
         },  
         "line" : { 
         "color" : color_list[4],
         "width" : 1 
         },  
         "name": "All Frequencies"
       },  
       {   
         "x" : daxi_output[0],
         "y" : daxi_output[2],
         "type" : "bar",
         "marker": {"color": color_list[4]},
         "showlegend": False,
         "xaxis": "x2",
         "yaxis": "y2"
       },  
       {   
         "x" : daxi_output[0],
         "y" : daxi_output[3],
         "type" : "bar",
         "marker": {"color": color_list[4]},
         "showlegend": False,
         "xaxis" : "x3",
         "yaxis" : "y3",
         "type"  : "bar"
       },  
       {   
         "x": daxi_output[0],
         "y": daxi_output[4],
         "xaxis": "x4",
         "yaxis": "y4",
         "mode": "lines+markers",
         "showlegend": False,
         "marker": {
           "color" : color_list[4],
           "size" : 8
         }
        }]    

    daxi_plotly ={
   "data" : daxi_data,
   "layout" : {
       "title" : "Dynamics: Axially Symmetric Model",
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
    return daxi_plotly

