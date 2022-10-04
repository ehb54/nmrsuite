import numpy as np
import math as m

def select_best_solution(res=None, cond_num_cutoff=None, Chi2_colmn=None, cond_num_colmn=None):

    nc = res.shape[1]
    if cond_num_colmn is None: cond_num_colmn = nc - 1 #default: last entry in results
    if Chi2_colmn is None: Chi2_colmn = nc - 2 #default: penultimate entry in results
    if cond_num_cutoff is None: cond_num_cutoff = 6 #default

    solution = np.array([])

    ind_cond = np.argwhere(res[:, cond_num_colmn] < cond_num_cutoff)
    if ind_cond.size > 0:
        Chi2min = min(res[ind_cond, Chi2_colmn])
        indmin = np.where(res[:, Chi2_colmn] == Chi2min)
        solution = np.copy(res[indmin[0], :]).flatten()
    return solution

                    #======================================================
