#!/opt/miniconda2/bin/python
import pandas as pd
import numpy as np
import sys
import re
def preprocess_out(start_ind,end_ind,all_txt, dat_type):
    part_data = all_txt[start_ind:end_ind]
    data_list = []
    for items in part_data:
        if len(items) > 0:
            data_list.append(items.split())
    if dat_type == "exp":
        col_name = ["residue","chain","atom1","atom2","freq","rho_exp","rho_pred","rho_error"]
    elif dat_type == "dyna":
        col_name = ["residue","chain","atom1","atom2","freq","csa","s2","s2_err","tau_loc","Rex","s2_fast","tau_fast","chi2"]
    else:
        print("Unsupported data type %s." % dat_type)
    output_pd = pd.DataFrame(data = data_list, columns = col_name)
    resi = list(output_pd["residue"])
    dyna_resi = []
    exp_resi = []
    ori_resi = []
    #deal with comment symbols in residue number
    for list_elem in resi:
	if list_elem.startswith("%"):
            pass
        elif list_elem.startswith("*"):
            dyna_resi.append(list_elem[1:])
            exp_resi.append(list_elem)
        else:
            dyna_resi.append(list_elem)
            exp_resi.append(list_elem) 
    ori_resi = dyna_resi

    if dat_type == "exp":
        output_pd = output_pd[[x in exp_resi for x in ori_resi]]
    else:
        output_pd = output_pd[[x in dyna_resi for x in ori_resi]]
    
    resi_pd = output_pd['residue']
    int_resi = [int(items) for items in resi_pd]
    output_pd['residue'] = int_resi
    sort_output = output_pd.sort_values(['freq','chain','residue'], ascending=[True,True,True]) 
    freq_set = set(list(sort_output["freq"]))
    chain_set = set(list(sort_output["chain"]))
    output_dict = {}
    output_keys = []
#no docking involved
    if len(chain_set) == 1:
        for items in freq_set:
            tmp = sort_output[sort_output["freq"]==items]    
            freq_tag = str(items) + "MHz"
            output_keys.append(freq_tag)
            output_dict[freq_tag] = tmp 
        all_rho_exp = []
        all_rho_pred = []
        all_rho_err = []
        all_rho_resi = []
        if dat_type =="exp":
            for items in output_keys:
                rho_resi_list = [int(ele) for ele in list(output_dict[items]["residue"])]    
                rho_exp_list = [float(ele) for ele in list(output_dict[items]["rho_exp"])]
                rho_pred_list = [float(ele) for ele in list(output_dict[items]["rho_pred"])]
                rho_err_list = [float(ele) for ele in list(output_dict[items]["rho_error"])]
                all_rho_resi.append(rho_resi_list)
                all_rho_exp.append(rho_exp_list)
                all_rho_pred.append(rho_pred_list)
                all_rho_err.append(rho_err_list)
            all_output = [all_rho_resi, all_rho_exp,all_rho_pred,all_rho_err]
        elif dat_type =="dyna":
        #all frequencies have the same dyanmics data
            resi_list = [int(ele) for ele in list(output_dict[output_keys[0]]["residue"])]
            s2_list = [float(ele) for ele in list(output_dict[output_keys[0]]["s2"])]
            tau_loc_list = [float(ele) for ele in list(output_dict[output_keys[0]]["tau_loc"])]
            rex_list = [float(ele) for ele in list(output_dict[output_keys[0]]["Rex"])]
            s2_fast_list = [float(ele) for ele in list(output_dict[output_keys[0]]["s2_fast"])]
            all_output =[resi_list, s2_list, tau_loc_list,rex_list,s2_fast_list]
    elif len(chain_set) == 2:
        for items in chain_set:
            tmp = sort_output[sort_output["chain"]==items]    
            chain_tag = str(items) + "MHz"
            output_keys.append(chain_tag)
            output_dict[chain_tag] = tmp 
        all_s2 = []
        all_tau_loc = []
        all_rex = []
        all_s2_fast = []
        all_rho_exp = []
        all_rho_pred = []
        all_rho_err = []
        all_rho_resi = []
        if dat_type =="exp":
            for items in output_keys:
                rho_resi_list = [int(ele) for ele in list(output_dict[output_keys[0]]["residue"])]    
                rho_exp_list = [float(ele) for ele in list(output_dict[items]["rho_exp"])]
                rho_pred_list = [float(ele) for ele in list(output_dict[items]["rho_pred"])]
                rho_err_list = [float(ele) for ele in list(output_dict[items]["rho_error"])]
                all_rho_resi.append(rho_resi_list)
                all_rho_exp.append(rho_exp_list)
                all_rho_pred.append(rho_pred_list)
                all_rho_err.append(rho_err_list)
            all_output = [all_rho_resi, all_rho_exp,all_rho_pred,all_rho_err]
        elif dat_type =="dyna":
        #different chains have different dyanmics data
            for items in output_keys:
                resi_list = [int(ele) for ele in list(output_dict[output_keys[0]]["residue"])]
                s2_list = [float(ele) for ele in list(output_dict[items]["s2"])]
                tau_loc_list = [float(ele) for ele in list(output_dict[items]["tau_loc"])]
                rex_list = [float(ele) for ele in list(output_dict[items]["Rex"])]
                s2_fast_list = [float(ele) for ele in list(output_dict[items]["s2_fast"])]
                all_rho_resi.append(resi_list)
                all_s2.append(s2_list)
                all_tau_loc.append(tau_loc_list)
                all_rex.append(rex_list)
                all_s2_fast.append(s2_fast_list)
            all_output =[all_rho_resi, all_s2,all_tau_loc,all_rex,all_s2_fast]
        
    else:
        sys.stderr.write("More than 2 chains or fewer than 1 chain are not supported\n")
    return all_output
