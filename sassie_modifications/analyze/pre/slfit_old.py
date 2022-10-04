import numpy as np
import math, sys, os
import scipy.optimize

from bokeh.layouts import layout, row, column
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import HoverTool, Span, Range1d, ColumnDataSource,FactorRange , tickers, SingleIntervalTicker, LinearAxis, WheelZoomTool, LabelSet, Label, OpenURL, TapTool, Div, LayoutDOM
from bokeh.transform import factor_cmap
from bokeh.embed import components

from sassie.analyze.pre.pdb2nhcoor import pdb2nhcoor
from sassie.analyze.pre.SL_add2pdb import SL_add2pdb
from sassie.analyze.pre.SLfit_ss import SLfit_ss

# For local test
#sys.path.append('./')
#from pdb2nhcoor import pdb2nhcoor
#from SL_add2pdb import SL_add2pdb
#from SLfit_ss import SLfit_ss

#from config import init

def slfit(self, app, ratio, T2dia, Htime, freq, TAUc, pdb_filename, scaling = 1.0, model = 1, reslist = np.array([]), chainID = [], out_filename = [], guess = np.array([])):
    mvars = self.mvars
    pgui = self.run_utils.print_gui

    omega = freq * 2.0 * math.pi * 1.0e6
    d2 = 1.23e-44
    d2 = d2 * (1.0e10 ** 6)
    tauC = TAUc * 1.0e-9
    R2dia = 1.0/T2dia
    factor = d2 * (4.0 * tauC + 3.0 * tauC / (1.0 + (omega * tauC) ** 2))

    position  = np.full((3,), float('nan'))
    ratio_dist = []
    Chi2 = float('nan')

    ratio[:, 1] = ratio[:, 1] * scaling
    rlist_ratio = ratio[:, 0]
    nres_ratio = len(rlist_ratio)

    NH_coord = pdb2nhcoor(self,pdb_filename, [], model, 0, chainID)
    NH_coord[:, 0] = np.subtract(NH_coord[:, 0], 1)
    nNH = NH_coord[:, 0].shape[0]
    if guess.size == 0:
        guess = np.mean(NH_coord[~np.isnan(NH_coord[:, 1]), 1:4], axis = 0)
        pgui('Center of mass was taken as an intial guess\n')

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

    pgui('Chi^2 = ' +  "{:.4f}".format(Chi2) + "\n")
    pgui('Position of the Spin Label in the PDB coordinate frame: \n')
    pgui(' X = ' + "{:.4f}".format(position[0]) + ' A\n')
    pgui(' Y = ' + "{:.4f}".format(position[1]) + ' A\n')
    pgui(' Z = ' + "{:.4f}".format(position[2])+ ' A\n\n')

    if len(out_filename) > 0:
        [s, w] = SL_add2pdb(pdb_filename, out_filename, position)

#   For Plotting
#   ToDo: Refactorize later

    Tools = 'pan,box_zoom,wheel_zoom,save,reset,crosshair'

#   Place holder for later
    html_dir = os.path.join(mvars.runname, app, 'html')

#    wtitle = "PRE results ( Model nubmber = " + str(frame + 1) + " )"
    wtitle = "PRE results"
    output_dir = os.path.join(mvars.runname, app)
    output_html = os.path.join(output_dir, mvars.runname + '_single_model.html')
    output_file(output_html, title=wtitle )
    colors = ["blue", "green", "cyon", "black", "yellow", "violet", "red"]

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

    source1 = ColumnDataSource(data=dict(
        resid = ratio[:,0],
        ratio = ratio[:,1]
        ))
    source1_predict = ColumnDataSource(data=dict(
        resid = coords2X,
        ratio = coords2Y
        ))
    source1_actualfit= ColumnDataSource(data=dict(
        resid = coords3X,
        ratio = coords3Y
        ))

    hover1 = HoverTool(
            tooltips = [
            ("Resid", "@resid"),
            ("Ratio", "@ratio"),
            ]
            )

    height = 400
    width = 600
    s1_title = "Ratio vs. Residue Number"
    s1 = figure(plot_height=height, plot_width=width, title=s1_title, tools=[hover1,Tools])

    s1.vbar(x='resid', top = 'ratio', width = 0.5, color = "blue", line_color = "black", legend = "Experiment", alpha = 0.7, source=source1)
    s1.circle(x='resid', y = 'ratio', size = 8, line_width = 2, color = "red", fill_color = "white",fill_alpha = 0, legend = "Prediction", source=source1_predict)
    s1.line(x='resid', y = 'ratio', color = "red", line_width = 2, source=source1_predict)
#    s1.asterisk(x='resid', y = 'ratio', size = 10, line_width = 2, color = "black", fill_color = "white",fill_alpha = 0, legend = "Actual fit", source=source1_predict)
#    s1.line(x='resid', y = 'ratio', color = "black", line_width = 2, source=source1_predict)

    y_range = 0, np.max(ratio[:,1])*1.2
    x_range = np.min(ratio[:,0]) - 2, np.max(ratio[:,0]) + 2
    s1.x_range = Range1d(*x_range)
    s1.y_range = Range1d(*y_range)

    s1.legend.orientation = "horizontal"
    s1.legend.location = "top_left"
    s1.legend.click_policy="hide"
    s1.yaxis.axis_label = "Ratio"
    s1.xaxis.axis_label = "Residue Number"

    dist_plot = np.arange(min(10, dist), max(40, dist) + 0.1, 0.1)
    R2para_plot = np.divide(factor, dist_plot ** 6)
    ratio_plot = R2dia * np.divide(np.exp(-R2para_plot * Htime), (R2para_plot + R2dia))

    source2 = ColumnDataSource(data=dict(
        distance = dist_plot,
        ratio = ratio_plot,
        resid = ratio[:,0]
        ))
    source2_predict = ColumnDataSource(data=dict(
        distance = ratio_dist[fit_ind, 2],
        ratio = fit_ratio[:, 1],
        resid = ratio[:,0]
        ))

    hover2 = HoverTool(
            tooltips = [
            ("Distance", "@distance"),
            ("Ratio", "@ratio"),
            ("Resid", '@resid'),
            ]
        )

    s2_title = "Ratio vs. Distance"
    s2 = figure(plot_height=height, plot_width=width, title=s2_title, tools=[hover2,Tools])

    s2.line(x='distance', y = 'ratio', color = "black", line_width = 2, legend = "Prediction", source=source2)
    s2.circle(x='distance', y = 'ratio', size = 8, line_width = 2, color = "red", fill_color = "red", fill_alpha = 0, legend = "Experiment", source=source2_predict)
#    s1.asterisk(x='resid', y = 'ratio', size = 10, line_width = 2, color = "black", fill_color = "white",fill_alpha = 0, legend = "Actual fit", source=source1_predict)
#    s1.line(x='resid', y = 'ratio', color = "black", line_width = 2, source=source1_predict)

    y_range = -0.05, 1.15
    x_range = np.min(dist_plot) - 2.0, np.max(dist_plot) + 2.0
    s2.x_range = Range1d(*x_range)
    s2.y_range = Range1d(*y_range)

    s2.legend.orientation = "horizontal"
    s2.legend.location = "top_left"
    s2.legend.click_policy="hide"
    s2.yaxis.axis_label = "Ratio"
    s2.xaxis.axis_label = "Distance [A]"

    a = sum((ratio_dist[fit_ind, 3] - ratio_dist[fit_ind, 1]) ** 2) / 2.0
    b = (ratio_dist[fit_ind, 3] - np.mean(ratio_dist[fit_ind, 3])) ** 2
    c = a / float(sum(b))

    qR_fact = round(np.sqrt(c), 5)

    R = np.corrcoef(ratio_dist[fit_ind, 3], ratio_dist[fit_ind, 1])
    r = round(R[0, 1], 5)

    source3 = ColumnDataSource(data=dict(
        experiment = ratio_dist[fit_ind, 3],
        prediction = ratio_dist[fit_ind, 1],
        resid = ratio[:,0]
        ))

    hover3 = HoverTool(
            tooltips = [
            ("Experiment", "@experiment"),
            ("Prediction", "@prediction"),
            ("Resid", '@resid'),
            ]
        )

    s3_title = "Experiment vs. Prediction (qR-factor = {}; CorrCoeff = {})".format(qR_fact, r)
    s3 = figure(plot_height=height, plot_width=width, title=s3_title, tools=[hover3,Tools])

    s3.line(x=[0,1], y = [0,1], color = "black", line_width = 2)
    s3.circle(x='experiment', y = 'prediction', size = 8, line_width = 2, color = "red", fill_color = "red", fill_alpha = 0, source=source3)

    x_range = -0.05, 1.15
    y_range = -0.05, 1.15
    s3.x_range = Range1d(*x_range)
    s3.y_range = Range1d(*y_range)

    s3.legend.orientation = "horizontal"
    s3.legend.location = "top_left"
    s3.legend.click_policy="hide"
    s3.yaxis.axis_label = "Prediction"
    s3.xaxis.axis_label = "Experiment"

    source4 = ColumnDataSource(data=dict(
        resid = coords2X,
        distance = ratio_dist[fit_ind, 2]
        ))

    hover4 = HoverTool(
            tooltips = [
            ("Residue Number", "@resid"),
            ("Distance", "@distance"),
            ("Resid", '@resid'),
            ]
        )

    s4_title = "Distance vs. Residue Number"
    s4 = figure(plot_height=height, plot_width=width, title=s4_title, tools=[hover4,Tools])

    s4.vbar(x='resid', top = 'distance', width = 0.5, color = "blue", line_color = "black", legend = "Distance", alpha = 0.7, source=source4)

    y_range = 0, np.max(ratio_dist[fit_ind, 2])*1.2
    x_range = np.min(ratio[:,0]) - 2, np.max(ratio[:,0]) + 2
    s4.x_range = Range1d(*x_range)
    s4.y_range = Range1d(*y_range)

    s4.legend.orientation = "horizontal"
    s4.legend.location = "top_left"
    s4.legend.click_policy="hide"
    s4.xaxis.axis_label = "Residue Number"
    s4.yaxis.axis_label = "Distance [A]"


    result = layout([[s1,s2],[s3,s4]])
    save(result)
    script,div = components(result)

    return script, div, position, ratio_dist, Chi2, qR_fact, r
