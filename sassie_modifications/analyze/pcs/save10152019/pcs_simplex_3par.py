import numpy as np

#import generate_vertices as genfile
#import PCS_tgt_3par as pcs3
#from col_sub import col_sub
#import vertindex as vi

import sassie.analyze.pcs.generate_vertices as genfile
import sassie.analyze.pcs.PCS_tgt_3par as pcs3
from sassie.analyze.pcs.col_sub import col_sub

def pcs_simplex_3par(self, x=None, data_PCS=None, coords=None, tolerance=None):
    #--------------------------------------------------------------------------
    # Performs simplex optimization by varying 3 parameters (x,y,z coordinates)
    # Optimization starts at an initial guess generated by generate_vertices
    # Uses the simplex algorithm to find a relative local minimum 
    #------------------------------------------------------------------------
    n = np.size(x)
    
    #Two tolerances are defined in case we wish to have separate values for the two
    tol = 1e-4 #Tolerance in parameter delta
    tol2 = float(tolerance) #Tolerance in function (chi-square) delta
    cnt_Max = 10000 #Maximum number of function evaluations
    
    
    # Set up a simplex near the initial guess.
    xin = np.copy(np.transpose(np.asarray(x)))# Force xin to be a column vector, AG: probably do not even need this, since Python does not care about vector orientation
    v = np.copy(xin)# Place input guess in the simplex! (credit L.Pfeffer at Stanford)
    fv = pcs3.PCS_tgt_3par(x, data_PCS, coords)

    # Following improvement suggested by L.Pfeffer at Stanford
    usual_delta = 0.05# 5 percent deltas for non-zero terms
    zero_term_delta = 0.00025# Even smaller delta for zero elements of x
    for j in range(n):
        y = np.copy(xin)
        if (y[j] != 0):
            y[j] = (1 + usual_delta) * (y[j])
        else:
            y[j] = zero_term_delta
        v = np.column_stack((v, y))
        x = np.copy(y)
        f = pcs3.PCS_tgt_3par(x, data_PCS, coords) #Calculate chi-square
        fv = np.copy(np.concatenate([fv, f])) #Append to list of chi-square
    j = np.argsort(fv)
    fv = np.sort(fv)
    v = np.asarray(v[:, j])
    
    cnt = n + 1
    alpha = 1.0
    beta = 1.0 / 2.0
    gamma = 2.0

    n = np.shape(v)[0]
    ot = np.arange(1, n+1, dtype=int)
    on = np.arange(n, dtype=int)

    # Iterate until the diameter of the simplex is less than tol.
    while cnt < cnt_Max:
        if np.max(abs(col_sub(v[:, ot], v[:, 0]))) <= abs(tol) and np.max(abs((fv[ot] - fv[0])/fv[ot])) <= abs(tol2): break
        # One step of the Nelder-Mead simplex algorithm, review the below section
        vbar = (np.sum(np.copy(v[:,on]), 1) / float(n)) #Sum each row
        vr = (1 + alpha) * np.copy(vbar) - alpha * np.copy(v[:, n])
        x[:] = np.copy(vr)
        fr = pcs3.PCS_tgt_3par(x, data_PCS, coords)
        cnt = cnt + 1
        vk = np.copy(vr)
        fk = np.copy(fr)
        how = 'reflect '

        if fr < fv[n-1]:
            if fr < fv[0]:
                ve = gamma * np.copy(vr) + (1 - gamma) * np.copy(vbar)
                x[:] = np.copy(ve)
                fe = pcs3.PCS_tgt_3par(x, data_PCS, coords)
                cnt = cnt + 1
                if fe < fv[0]:
                    vk = np.copy(ve)
                    fk = np.copy(fe)

                    how = 'expand'
        else:
            vt = np.copy(v[:, n])
            ft = np.copy(fv[n])

            if fr < ft:
                vt = np.copy(vr)
                ft = np.copy(fr)

            
            vc = beta * np.copy(vt) + (1 - beta) * np.copy(vbar)
            x[:] = np.copy(vc)
            fc = pcs3.PCS_tgt_3par(x, data_PCS, coords)
            cnt = cnt + 1
            if fc < fv[n-1]:
                vk = np.copy(vc)
                fk = np.copy(fc)

                how = 'contract'
            else:
                for j in range(1, n):
                    v[:, j] = (np.copy(v[:, 0]) + np.copy(v[:, j])) / 2.0
                    x[:] = np.copy(v[:, j])
                    fv[j] = pcs3.PCS_tgt_3par(x, data_PCS, coords)
                
                cnt = cnt + n - 1
                vk = (np.copy(v[:, 0]) + np.copy(v[:, n])) / 2.0
                x[:] = np.copy(vk) 
                fk = pcs3.PCS_tgt_3par(x, data_PCS, coords)
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
        print("increase Maximum function evaluations")
    else:
        print("Done")

    return [x, FVmin]
    #============================================================

'''
f = np.loadtxt("../testset_Tm/coordinates_H.dat")
coordH = f[:,1:4]
coord = np.loadtxt("../testset_Tm/coordinates_all.dat")
g = np.loadtxt("../testset_Tm/PCS_data.dat")
pcs = g[:,1]
guess = np.array([66.1, -95.1450, -25.2660])
#guess=np.array([63.7638, -91.6308, -11.426])
vertices = genfile.generate_vertices(coord)
output = pcs_simplex_3par(guess, pcs, coordH)
print (output)
'''
