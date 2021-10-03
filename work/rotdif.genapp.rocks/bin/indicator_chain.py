#!/usr/bin/python
import os,sys
import pandas as pd
import re

def dock(relax_loc):
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
    if len(chain_set) == 0:
        to_dock = False
        sys.stderr.write("must be at least 1 chain")
    elif len(chain_set) == 1:
        to_dock = False
    elif len(chain_set) == 2:
        to_dock = True
    elif len(chain_set) > 2:
        sys.stderr.write("More than 2 chains are not supported")
    return to_dock
    
