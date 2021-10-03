#!/usr/bin/python
from preprocess_split import split_data

def save_elm(dyna_flag, elm_flag, elmdock_flag):
    elm_pred = split_data(dyna_flag, elm_flag, elmdock_flag)[6]
    with open('ELM_prediction','w') as out_f:
        for items in elm_pred:
            out_f.write(items+'\n')

