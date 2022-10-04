import numpy as np
import math

from pdb2nhcoor import pdb2nhcoor
from SL_add2pdb import SL_add2pdb
from SLfit_ss import SLfit_ss
import scipy.optimize

#from config import init

def slfit(ratio, T2dia, Htime, freq, TAUc, pdb_filename, out_filename = [], guess = np.array([]), reslist = np.array([]), scaling = 1, model = 1, chainID = []):
    if scaling == None:
        scaling = 1

    omega = freq * 2 * math.pi * 1e6   
    d2 = 1.23e-44
    d2 = d2 * (1e10 ** 6)
    tauC = TAUc * 1e-9
    R2dia = 1/float(T2dia)
    factor = d2 * (4 * tauC + 3 * tauC / float(1 + (omega * tauC) ** 2))

    position  = np.full((3,), float('nan'))
    ratio_dist = []
    Chi2 = float('nan')

    ratio[:, 1] = ratio[:, 1] * scaling
    rlist_ratio = ratio[:, 0]
    nres_ratio = len(rlist_ratio)

    NH_coord = pdb2nhcoor(pdb_filename, [], model, 0, chainID)
    NH_coord[:, 0] = np.subtract(NH_coord[:, 0], 1)
    nNH = NH_coord[:, 0].shape[0]
    if guess.size == 0:
        guess = np.mean(NH_coord[~np.isnan(NH_coord[:, 1]), 1:4], axis = 0)
        print('guess = COM')
        
    ind_ratio = np.full((nres_ratio, 1), float('nan'))
    ind_NH = np.full((nNH, 1), float('nan'))
    for i in range(nres_ratio):
        if reslist.size > 0:
            indB = np.nonzero(reslist[:] == rlist_ratio[i])[0]
        else:
            indB = np.array([1])
        
        if ~np.isnan(ratio[i, 1]) and indB.size > 0:
            indA = np.nonzero(NH_coord[:, 0] == rlist_ratio[i])[0]
            if indA.size > 0: 
                if ~(np.isnan(NH_coord[indA, 1::]).all()):
                    ind_NH[indA] = indA
                    ind_ratio[i] = i
                
    ind_ratio = ind_ratio[~np.isnan(ind_ratio)].astype(int)
    ind_NH = ind_NH[~np.isnan(ind_NH)].astype(int)
    
    fit_coord = NH_coord[ind_NH][:, [0, 4, 5, 6]]
    fit_ratio = ratio[ind_ratio, :]

    
    [position, Chi2, iter, funcalls, warnflag] = scipy.optimize.fmin(func = SLfit_ss, x0 = guess, args = (fit_ratio, fit_coord, factor, R2dia, Htime), xtol = 1e-5, ftol = 1e-5, maxiter = 100000, maxfun = 1000000, full_output = True, disp = False)

    print(' ')
    print(['Chi**2 = ', str(Chi2)])
    print(' ')
    print('Position of the Spin Label in the PDB coordinate frame: ')
    print([' X = ', str(position[0]), ' A'])
    print([' Y = ', str(position[1]), ' A'])
    print([' Z = ', str(position[2]), ' A'])
    print(' ')

    if len(out_filename) > 0:
        [s, w] = SL_add2pdb(pdb_filename, out_filename, position)

    nres = NH_coord.shape[0]
    ratio_dist = np.full((nres, 4), float('nan'))
    fit_ind = []
    not_fit_ind = []
    for i in range(nres): 
        Hvect = NH_coord[i, 4:7] - position
        dist = np.sqrt(np.dot(Hvect, np.conj(Hvect.reshape(Hvect.size, 1))))
        R2para = factor/(dist ** 6)
        ratio_sim = R2dia * math.exp(-R2para * Htime) / float(R2para + R2dia)
        ratio_dist[i, 0] = NH_coord[i, 0]
        ratio_dist[i, 1] = ratio_sim
        ratio_dist[i, 2] = dist
        indnn = np.nonzero(fit_ratio[:, 0] == ratio_dist[i, 0])[0]
        if indnn.size > 0:
            ratio_dist[i, 3] = fit_ratio[indnn, 1]
            fit_ind.append(i)
        else:
            not_fit_ind.append(i)

    import matplotlib.pyplot as plt

    fit_ind = np.asarray(fit_ind).flatten()
    not_fit_ind = np.asarray(not_fit_ind).flatten()

    coords1X = ratio_dist[:, 0]
    coords1Y = ratio_dist[:, 1]
    coords2X = ratio_dist[fit_ind, 0]
    coords2Y = ratio_dist[fit_ind, 1]
    coords3X = ratio_dist[not_fit_ind, 0]
    coords3Y = ratio_dist[not_fit_ind, 1]

    combinationX = np.concatenate([coords1X, coords2X])
    combinationY = np.concatenate([coords1Y, coords2Y])

    combinationY = combinationY[np.argsort(combinationX)]
    combinationX = np.sort(combinationX)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(ratio[:, 0], ratio[:, 1], zorder = 1)
    ax.plot(combinationX, combinationY, linewidth = 0.25, color = "red", zorder = 2)
    ax.scatter(coords2X, coords2Y, edgecolors = "r", zorder = 3)
    ax.scatter(coords3X, coords3Y, facecolors = "none", edgecolors = "r", zorder = 4)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Ratio")
    ax.set_title("bars = experiment, o = prediction, * = actual fit")
    plt.savefig('fig1.png')


    dist_plot = np.arange(min(10, dist), max(40, dist) + 0.1, 0.1)
    R2para_plot = np.divide(factor, dist_plot ** 6)
    ratio_plot = R2dia * np.divide(np.exp(-R2para_plot * Htime), (R2para_plot + R2dia))

    fig = plt.figure('figure2')
    ax = fig.add_subplot(111)
    ax.plot(dist_plot, ratio_plot, color = "black", zorder = 1)
    ax.scatter(ratio_dist[fit_ind, 2], fit_ratio[:, 1], color = "red", zorder = 2)
    ax.set_xlabel("Distance, A")
    ax.set_ylabel("Ratio")
    ax.set_title("o = experiment, line = prediction")
    plt.savefig('fig2.png')

    a = sum((ratio_dist[fit_ind, 3] - ratio_dist[fit_ind, 1]) ** 2) / 2.0
    b = (ratio_dist[fit_ind, 3] - np.mean(ratio_dist[fit_ind, 3])) ** 2
    c = a / float(sum(b))

    qR_fact = round(np.sqrt(c), 5)

    R = np.corrcoef(ratio_dist[fit_ind, 3], ratio_dist[fit_ind, 1])
    r = round(R[0, 1], 5)

    fig = plt.figure('figure4')
    ax = fig.add_subplot(111)
    ax.plot([0, 1], [0, 1], color = "red", zorder = 1)
    ax.scatter(ratio_dist[fit_ind, 3], ratio_dist[fit_ind, 1], color = "blue", zorder = 2)
    if (scaling == 1):
        ax.set_xlabel("Experiment")
    else:
        ax.set_xlabel("Experiment")

    ax.set_ylabel("Prediction")

    ax.set_title("qR-factor = {}; CoorCoeff = {}".format(qR_fact, r))

    plt.savefig('fig3.png')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(coords2X, ratio_dist[fit_ind, 2])
    ax.set_xlabel("Residue")
    ax.set_ylabel("Distance, A")
    plt.savefig('fig4.png')

    return [position, ratio_dist, Chi2]
