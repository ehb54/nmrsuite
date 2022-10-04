import os
import re
import shutil
import plotly
import json
import time
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


L0_folder_name = "LO_Solutions_New"

def printQuit (string):
    output = {}
    output["_textarea"] = str(string)
    print (json.dumps(output))
    quit()

def L0Plot (L0_number):
    os.chdir(L0_folder_name)

    if (L0_number == 0):
        df = pd.read_csv("BestPossibleSolution.txt", delimiter="\t")
    else:
        df = pd.read_csv("L0={}.txt".format(L0_number), delimiter="\t")

    array = df.values

    index = array[:, 0].astype(int)
    exp_value = array[:, 1].astype(float)
    pred_value = array[:, 2].astype(float)

    #output = f"Index: {index}; Type: {type(index)}.        Exp Value Type :{type(exp_value)}"
    #printQuit(output)
    
    # Finds the min values and max values for the graph
    min_x = np.min(exp_value)
    min_y = np.min(pred_value)
    if min_x < min_y:
        min_abs = min_x
        max_abs = np.max(exp_value)
    else:
        min_abs = min_y
        max_abs = np.max(pred_value)

    # Determines correlation coefficient, prints error if calculation fails for whatever reson
    try:
        corr_coeff = np.corrcoef(exp_value, pred_value)[1, 0]
    except:
        output = f"Exp Value: {exp_value}       Pred Value: {pred_value}"
        printQuit(output)

    # Calculates relative error
    a = sum(np.subtract(exp_value, pred_value) ** 2)
    b = exp_value ** 2
    c = a / float(sum(b))
    rel_error = np.sqrt(c)
    #rel_error = np.linalg.norm(np.divide(np.subtract(exp_value, pred_value), pred_value))

    os.chdir("..")

    # Prints coor coeff and rel error, the two things we need    
    return exp_value.tolist(), pred_value.tolist(), index.tolist(), min_abs, max_abs, corr_coeff, rel_error#, corr_coeff

#saveL0Solutions(open("output.txt", "r").readlines())
#L0Plot(1)

def saveAndPlot(L0_numbers, Chi2L_numbers, rel_error_numbers):
    #L0_numbers, Chi2L_numbers, rel_error_numbers = saveL0Solutions(output)
    num_L0s = len(L0_numbers)

    L0_to_Chi2L =  dict(zip(L0_numbers, Chi2L_numbers))
    L0_to_rel_error =  dict(zip(L0_numbers, rel_error_numbers))

    height = num_L0s * 500
    
    subplot_titles = []
    for L0_number in L0_numbers:
        if (L0_number == 0):
            subplot_titles.append("Best possible x>0 solution")
        else:
            subplot_titles.append("L0-norm={}".format(L0_number))


    fig = make_subplots(rows=num_L0s, cols=1, subplot_titles=tuple(subplot_titles)) # Sets up graph

    shapes = []

    #rel_error_min = min(rel_error_arr)

    additions_to_titles = {} # L0_number: [R, Q]

    for i, L0_number in enumerate(L0_numbers):
        exp_value, pred_value, index, min_abs, max_abs, corr_coeff, rel_error = L0Plot(L0_number)
        additions_to_titles[L0_number] = [corr_coeff, rel_error] # Adds correlation coefficient and relative error to list so those two variables can be shown for each graph

        fig.add_trace(go.Scatter(x = exp_value, y = pred_value, mode = 'markers', text = index), row = i+1, col = 1)

        fig.update_xaxes(title_text="Experimental Data", row=i+1, col=1)
        fig.update_yaxes(title_text="Predicted Data", row=i+1, col=1)

        xref = "x{}".format(i+1)
        yref = "y{}".format(i+1)

        shapes.append({'type': 'line', 'x0': min_abs, 'y0': min_abs, 'x1': max_abs, 'y1': max_abs, 'xref': xref, 'yref': yref})

    #fig['layout'].update(shapes=shapes)

    fig.update_layout(height=height, width=1000, title_text="Agreement between Experimental and Predicted Data", shapes=shapes, showlegend = False)
    fig.update_traces(marker=dict(color='Green'))


    for annotation in fig["layout"]["annotations"]:
        text = annotation["text"]
        if (text == "Best possible x>0 solution"):
            L0_number = 0
        else:
            L0_number = int(re.sub("[^0-9]", "", text))
        
        [R, Q] = additions_to_titles[L0_number]
        Chi2L = L0_to_Chi2L[L0_number]
        rel_error = L0_to_rel_error[L0_number]

        annotation["text"] += f" (Corr.Coeff={round(R, 5)}; Q-factor={round(Q, 5)}; Rel.Error={round(rel_error, 5)}; Chi2/L={round(Chi2L, 5)})"

    fig_data = []
    for i in range(num_L0s, 0, -1):
        fig_data.append(fig.data[-i])

    """ if (0 in L0_numbers):
        fig_data.append(fig_data.pop(0))  """

    fig.data = fig_data
    #fig["data"] = fig["data"].append(fig["data"].pop(0))

    return fig.to_dict()

    # Look up "alias" |-><-| check running a .sh file (e.g. mygit.sh) .bat instead of .sh  also command lit
