#!/usr/bin/python
import operator
from preprocess_split import split_data
import statistics as stat
def axi_plot(exp_resi,color_list, exp_keys, dyna_flag, elm_flag, elmdock_flag):
    axi_output = split_data(dyna_flag, elm_flag, elmdock_flag)[1]
    axi_block = []

    ind_a = []
    to_index = axi_output[0][0]
    for ii in range(len(exp_resi)):
        ind_a.append(to_index.index(int(exp_resi[ii])))
    
    to_loop = len(axi_output[0])

    for i in range(to_loop):
        axi_output[0][i] = [axi_output[0][i][jj] for jj in ind_a]
        axi_output[1][i] = [axi_output[1][i][jj] for jj in ind_a]
        axi_output[2][i] = [axi_output[2][i][jj] for jj in ind_a]
        axi_output[3][i] = [axi_output[3][i][jj] for jj in ind_a]    

        err_dev = stat.stdev(axi_output[3][i])
        dif = list(map(operator.sub,axi_output[2][i],axi_output[1][i]))
        tmp_data = [
        {
          "x": axi_output[1][i],
          "y": axi_output[2][i],
          "text": ["residue "+ str(int(items)) for items in axi_output[0][i]],
          "error_x": {
            "array" : axi_output[3][i], 
            "visible" : False
},
          "mode": "markers",
          "marker": {
            "color": color_list[i],
            "size": 12
          },
          "name" : exp_keys[i],
          "legendgroup" : str(i),
          "showlegend" : True
        },
        {
          "x" : axi_output[1][i],
          "y" : axi_output[1][i],
          "mode" : "lines",
          "line" : {
            "color" : color_list[i],
            "width": 3
          },
          "name" : exp_keys[i],
          "legendgroup" : str(i),
          "showlegend" : False
        },
        {
          "x": axi_output[0][i],
          "y": [x / err_dev for x in dif],
          "type": "bar",
          "marker": {"color": color_list[i]},
          "xaxis": "x2",
          "yaxis": "y2",
          "name" : exp_keys[i],  
          "legendgroup" : str(i),
          "showlegend" : False
        }
        ]
        axi_block.extend(tmp_data)

    axi_plotly = { 
    "data" : axi_block, 
    "layout" : { 
        "title" : "Axially Symmetric Model Fit",
        "xaxis" : {
        "domain": [0,0.4],
        "title" : "Rho_exp"
        },
        "yaxis" : {
        "title" : "Rho_pred"
        },
        "yaxis2" : {
        "anchor" : "x2",
        "title" : "Residuals/sigma"
        },
        "xaxis2" : {
        "domain" : [0.6,1],
        "title" : "residue number"
        }     
    }
    }
    return axi_plotly
   
