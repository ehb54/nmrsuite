#!/usr/bin/python
import operator
from preprocess_split import split_data
import statistics as stat

def iso_plot(exp_resi,color_list, exp_keys, dyna_flag, elm_flag, elmdock_flag):
    iso_output = split_data(dyna_flag, elm_flag, elmdock_flag)[0]
    iso_block = []
    ind_a = []
    to_index = iso_output[0][0]
    for ii in range(len(exp_resi)):
        ind_a.append(to_index.index(int(exp_resi[ii])))
    
    to_loop = len(iso_output[0])

    for i in range(to_loop):
        iso_output[0][i] = [iso_output[0][i][jj] for jj in ind_a]
        iso_output[1][i] = [iso_output[1][i][jj] for jj in ind_a]
        iso_output[2][i] = [iso_output[2][i][jj] for jj in ind_a]
        iso_output[3][i] = [iso_output[3][i][jj] for jj in ind_a]        

        err_dev = stat.stdev(iso_output[3][i])
        dif = list(map(operator.sub,iso_output[2][i],iso_output[1][i]))
        tmp_data = [
        {
          "x": iso_output[1][i],
          "y": iso_output[2][i],
          "text": ["residue "+ str(int(items)) for items in iso_output[0][i]],
          "error_x": {
            "array" : iso_output[3][i], 
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
          "x" : iso_output[1][i],
          "y" : iso_output[1][i],
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
          "x": iso_output[0][i],
          "y": [x / err_dev for x in dif],
          "type": "bar",
          "xaxis": "x2",
          "yaxis": "y2",
          "name" : exp_keys[i],
          "marker": {"color": color_list[i]},
          "legendgroup" : str(i),
          "showlegend" : False
        }
        ]
        iso_block.extend(tmp_data)

    iso_plotly = { 
    "data" : iso_block, 
    "layout" : { 
        "title" : "Isotropic Model Fit",
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
    return iso_plotly
   
