#!/opt/miniconda2/bin/python
from preprocess_split import split_data

def save_dock(dyna_flag, elm_flag, elmdock_flag):
    elm_pred = split_data(dyna_flag, elm_flag, elmdock_flag)[7]
    with open('elmdock_transformations.out','w') as out_f:
        for items in elm_pred:
            out_f.write(items+'\n')

