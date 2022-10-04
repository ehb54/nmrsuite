import numpy as np

#import generate_vertices as genfile
#import PCS_tgt_8par as pcs8
#from col_sub import col_sub

import sassie.analyze.pcs.generate_vertices as genfile
import sassie.analyze.pcs.PCS_tgt_8par as pcs8
from sassie.analyze.pcs.col_sub import col_sub

def pcs_simplex_8par(self, x=None, data_PCS=None, coords=None, tolerance = None):
    #--------------------------------------------------------------------------
    # This program performs simplex to minimize the calculated PCS values
    # Varies 8 parameters: the x,y,z coordinates of the spin label and the
    # five values of the X-tensor
    #------------------------------------------------------------------------
    n = np.size(x)
    
    #Fix below code to reflect new options data structure
    tol = 1e-4 #Tolerance in parameter delta
    tol2 = float(tolerance) #Tolerance in function (chi-square) delta
    cnt_Max = 10000 #Maximum number of function evaluations
    
    # Set up a simplex near the initial guess.
    xin = np.copy(np.transpose(np.asarray(x)))# Force xin to be a column vector, AG: probably do not even need this, since Python does not care about vector orientation
    v = np.copy(xin)# Place input guess in the simplex! (credit L.Pfeffer at Stanford)
    fv = pcs8.PCS_tgt_8par(x, data_PCS, coords)

    # Following improvement suggested by L.Pfeffer at Stanford
    usual_delta = 0.05# 5 percent deltas for non-zero terms
    zero_term_delta = 0.00025# Even smaller delta for zero elements of x
    for j in range(n):
        y = np.copy(xin)
        if (y[j] != 0):
            y[j] = (1 + usual_delta) * np.copy(y[j])
        else:
            y[j] = zero_term_delta
        v = np.column_stack((v, y))
        x = np.copy(y)
        f = pcs8.PCS_tgt_8par(x, data_PCS, coords)
        fv = np.copy(np.concatenate([fv, f]))
    j = np.argsort(fv)
    fv = np.sort(fv)
    v = np.asarray(v[:, j])
    cnt = n + 1
    alpha = 1.0
    beta = 1.0 / 2.0
    gamma = 2.0

    n = np.shape(v)[0]
    ot = np.arange(1, n+1, dtype=int) #In python 3, this will turn
    on = np.arange(n, dtype=int)
    
    # Iterate until the diameter of the simplex is less than tol.
    while cnt < cnt_Max:
        if np.max(abs(col_sub(v[:, ot], v[:, 0]))) <= tol and np.max(abs(fv[ot] - fv[0])/fv[ot]) <= tol2: break
        # One step of the Nelder-Mead simplex algorithm, review the below section
        vbar = (np.sum(np.copy(v[:,on]),1) / float(n)) 
        vr = (1 + alpha) * np.copy(vbar) - alpha * np.copy(v[:, n])
        x[:] = np.copy(vr)
        fr = pcs8.PCS_tgt_8par(x, data_PCS, coords)
        cnt = cnt + 1
        vk = np.copy(vr)
        fk = np.copy(fr)
        how = 'reflect '

        if fr < fv[n-1]:
            if fr < fv[0]:
                ve = gamma * np.copy(vr) + (1 - gamma) * vbar
                x[:] = np.copy(ve)
                fe = pcs8.PCS_tgt_8par(x, data_PCS, coords)
                cnt = cnt + 1
                if fe < fv[0]:
                    vk = np.copy(ve)
                    fk = np.copy(fe)

                    how = 'expand  '
        else:
            vt = np.copy(v[:, n])
            ft = np.copy(fv[n])

            if fr < ft:
                vt = np.copy(vr)
                ft = np.copy(fr)

            
            vc = beta * np.copy(vt) + (1.0 - beta) * vbar
            x[:] = np.copy(vc)
            fc = pcs8.PCS_tgt_8par(x, data_PCS, coords)
            cnt = cnt + 1
            if fc < fv[n-1]:
                vk = np.copy(vc)
                fk = np.copy(fc)

                how = 'contract'
            else:
                for j in range(1, n):
                    v[:, j] = (np.copy(v[:, 0]) + np.copy(v[:, j])) / 2.0
                    x[:] = np.copy(v[:, j])
                    fv[j] = pcs8.PCS_tgt_8par(x, data_PCS, coords)
                
                cnt = cnt + n - 1
                vk = (np.copy(v[:, 0]) + np.copy(v[:, n])) / 2.0
                x[:] = np.copy(vk)
                fk = pcs8.PCS_tgt_8par(x, data_PCS, coords)
                cnt = cnt + 1
                how = 'shrink  '
            
        v[:, n] = np.copy(vk)
        fv[n] = np.copy(fk)
        j = np.argsort(fv)
        fv = np.sort(fv)
        v = v[:, j]
    
    v[:, 0]
    FVmin = fv[0]#return function value at the minimum

    if cnt == cnt_Max:
        print("Warning: Maximum number of iterations (", int(cnt_Max), ") has been exceeded")
        print("increase Maximum Function Evaluations")
    else:
        print("Done")

    return [x, FVmin]
    #============================================================
