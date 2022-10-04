from bokeh.layouts import layout, row, column
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import HoverTool, Span, Range1d, ColumnDataSource,FactorRange , tickers, SingleIntervalTicker, LinearAxis, WheelZoomTool, LabelSet, Label, OpenURL, TapTool, Div, LayoutDOM
from bokeh.transform import factor_cmap
from bokeh.embed import components

import string, locale, re, os, json
import numpy as np

import sasmol.sasmol as sasmol

Tools = 'pan,box_zoom,wheel_zoom,save,reset,crosshair'
Tools_mc = 'pan,box_zoom,wheel_zoom,save,reset,crosshair,hover'
Tools_summary = 'pan,box_zoom,wheel_zoom,save,reset,crosshair, tap'

def multi_model_plot(self, height, width, app):

    mvars = self.mvars

    for nf in range(mvars.number_of_frames + 1):
       single_model_link(self, nf, 400, 500, app) 
  
    legends, rdata = read_out_summary(self,app)

    output_dir = os.path.join(mvars.runname, app)
    output_html = os.path.join(output_dir, mvars.runname + '_results_per_model.html')

    wtitle = "ALTENS results ( Total number of models = " + str(mvars.number_of_frames ) + " )"
    output_file(output_html, title=wtitle )

    s1_title = "R-factor vs. Correlation Cofficient (Total Number of Models = " + str(mvars.number_of_frames) + " )"
    s2_title = "Alignment Tensor Eigenvalues vs. Model Number"
    s3_title = "Euler Angles vs. Model Number"

    source1 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[1],
            y = rdata[2],
        ))

    source2_3 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[0],
            y = rdata[3],
        ))
    source2_4 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[0],
            y = rdata[4],
        ))
    source2_5 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[0],
            y = rdata[5],
        ))

    source3_6 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[0],
            y = rdata[6],
        ))
    source3_7 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[0],
            y = rdata[7],
        ))

    source3_8 = ColumnDataSource(data=dict(
            modelid = rdata[0],
            x = rdata[0],
            y = rdata[8],
        ))

    hover1 = HoverTool(
            tooltips = [
            ("Model", "@modelid"),
            ("Corr. Coeff", "@x"),
            ("R-factor", "@y{(0.0000)}"),
            ]
            )

    hover2 = HoverTool(
            tooltips = [
            ("Model", "@modelid"),
            ("Eigenvalue", "@y"),
            ]
            )
    hover3 = HoverTool(
            tooltips = [
            ("Model", "@modelid"),
            ("Angle", "@y"),
            ]
            )

    s1 = figure(plot_height=height, plot_width=width, title=s1_title, tools=[hover1,Tools_summary])
    s2 = figure(plot_height=height, plot_width=width, title=s2_title, tools=[hover2,Tools_summary]) 
    s3 = figure(plot_height=height, plot_width=width, title=s3_title, tools=[hover3,Tools_summary])

    xmin = min(rdata[1])
    xmax = max(rdata[1])
    ymin = min(rdata[2])
    ymax = max(rdata[2])

    xave = 0.0
    yave = 0.0

    for i in range(len(rdata[0])):
        xave += rdata[1][i]
        yave += rdata[2][i]

    xave = xave/float(len(rdata[0]))
    yave = yave/float(len(rdata[0]))

    source1_ave = ColumnDataSource(data=dict(
            names = ["Average"],
            modelid = ["Average"],
            x = [xave],
            y = [yave],
        ))

    x_range = 0.7, 1.05 
    y_range = -0.05, 0.35
    s1.x_range = Range1d(*x_range)
    s1.y_range = Range1d(*y_range)
    s1.circle(x='x', y='y', size = 8, color = "red", alpha = 0.8, source=source1)

    s1.cross(x='x', y='y', size = 16, color = "blue", line_width = 6.0, source=source1_ave)

    labels = LabelSet(x='x', y='y', text='names', level='glyph', 
              x_offset=-80, y_offset=5, source=source1_ave, render_mode='canvas') 

    s1.add_layout(labels)

    s1.legend.orientation = "horizontal"
    s1.legend.location = "top_left"
    s1.legend.click_policy="hide"
    s1.yaxis.axis_label = "R-factor"
    s1.xaxis.axis_label = "Correlation coefficient"

    x_range = min(rdata[0])*0.8,max(rdata[0])*1.1
    ymin = min(rdata[3]+rdata[4]+rdata[5])
    ymax = max(rdata[3]+rdata[4]+rdata[5])
    if (ymin > 0 and ymax > 0):
        y_range = ymin*0.9, ymax*1.2
    elif ( ymin < 0 and ymax > 0):
        y_range = ymin*1.1, ymax*1.8
    else:
        y_range = ymin*1.1, ymax*0.9

    s2.x_range = Range1d(*x_range)
    s2.y_range = Range1d(*y_range)
    legends[3] = "Axx"
    legends[4] = "Ayy"
    legends[5] = "Azz"
    s2.line(x='x', y='y', color = "blue", line_width = 2, legend = legends[3], source=source2_3)
    s2.line(x='x', y='y', color = "red", line_width = 2, legend = legends[4],source=source2_4)
    s2.line(x='x', y='y', color = "green", line_width = 2, legend = legends[5],source=source2_5)

    s2.circle(x='x', y='y', size = 10, fill_color = "white", color = "blue", source=source2_3)
    s2.circle(x='x', y='y', size = 10, fill_color = "white", color = "red", source=source2_4)
    s2.circle(x='x', y='y', size = 10, fill_color = "white", color = "green", source=source2_5)
 
    s2.legend.orientation = "horizontal"
    s2.legend.location = "top_left"
    s2.legend.click_policy="hide"
    s2.xaxis.axis_label = "Model number"
    s2.yaxis.axis_label = "Eigenvalues [*10^-3]"

    x_range = min(rdata[0])*0.8,max(rdata[0])*1.1
    ymin = min(rdata[6]+rdata[7]+rdata[8])
    ymax = max(rdata[6]+rdata[7]+rdata[8])
    if (ymin > 0 and ymax > 0):
        y_range = ymin*0.9, ymax*1.2
    elif ( ymin < 0 and ymax > 0):
        y_range = ymin*1.1, ymax*1.8
    else:
        y_range = ymin*1.1, ymax*0.9
    s3.x_range = Range1d(*x_range)
    s3.y_range = Range1d(*y_range)

#   Override legends for Euler angles  
    legends[6] = "Alpha"
    legends[7] = "Beta"
    legends[8] = "Gamma"

    s3.line(x='x', y='y', color = "blue", line_width = 2, legend = legends[6],source=source3_6)
    s3.line(x='x', y='y', color = "red", line_width = 2, legend = legends[7],source=source3_7)
    s3.line(x='x', y='y', color = "green", line_width = 2, legend = legends[8],source=source3_8)

    s3.circle(x='x', y='y', size = 10, fill_color = "white", color = "blue",source=source3_6)
    s3.circle(x='x', y='y', size = 10, fill_color = "white", color = "red", source=source3_7)
    s3.circle(x='x', y='y', size = 10, fill_color = "white", color = "green", source=source3_8)

    s3.legend.orientation = "horizontal"
    s3.legend.location = "top_left"
    s3.legend.click_policy="hide"
    s3.xaxis.axis_label = "Model number"
    s3.yaxis.axis_label = "Angles [degrees]"

    path_array = mvars.pdbfile.split('/')
    path_index =  path_array.index( "results" )
    path_array_slice = path_array[ path_index: -1]
    new_url =  '/'.join( path_array_slice) + '/' + output_dir + '/html/' + '@modelid' + '.html'

    url = new_url

    taptool1 = s1.select(type=TapTool)
    taptool1.callback = OpenURL(url=url)
    taptool2 = s2.select(type=TapTool)
    taptool2.callback = OpenURL(url=url)
    taptool3 = s3.select(type=TapTool)
    taptool3.callback = OpenURL(url=url)

    result = layout( [[s1, s2, s3]] )

    save(result)
#    show(result)

#   For embedding
    script, div = components(result)
    return script, div
 
def single_model_plot(self, frame, height, width, app):

    mvars = self.mvars
    output_dir = os.path.join(mvars.runname, app)
    if (frame == 0):
        output_html = os.path.join(output_dir, mvars.runname + '_single_model_first.html')
    elif (frame == mvars.number_of_frames - 1):
        output_html = os.path.join(output_dir, mvars.runname + '_single_model_last.html')
#    elif (frame == mvars.number_of_frames):
#        output_html = os.path.join(output_dir, mvars.runname + '_single_model_average.html')
    else:
        output_html = os.path.join(output_dir, mvars.runname + '_single_model.html')    


    wtitle = "PCS results ( Model nubmber = " + str(frame + 1) + " )" 
#    if ( frame == mvars.number_of_frames):
#        wtitle = "ALTENS results averaged over " + str(frame) + " models"

    output_file(output_html, title=wtitle )

    colors = ["blue", "green", "cyon", "black", "yellow", "violet", "red"]

    type_id, chain_all, res_all, rdc_exp_all, rdc_calc_all, delta_rdc_all = read_calc_pcs(self,frame,app)
    corr_trj,rfactor_trj = read_data_per_model(self, app)
    
    ntype = len(type_id)

#  Initializing configuration of plotting

    source1 = ["" for i in range(ntype)]
    source2 = ["" for i in range(ntype)]
    source3 = ["" for i in range(ntype)]
    source4 = ["" for i in range(ntype)]

    for i in range(ntype):
        source1[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = rdc_exp_all[i],
        ))

    hover1 = HoverTool(
            tooltips = [
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("PCS_exp", "@y"),
            ]
            )
    hover2 = HoverTool(
            tooltips = [
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("Delta", "@y"),
            ]
            )
    hover3 = HoverTool(
            tooltips = [
            ("Type", "@title"),
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("PCS_exp", "@x"),
            ("PCS_recalc", "@y"),
            ]
            )

    hover4 = HoverTool(
            tooltips = [
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("SL_distance", "@y"),
            ]
            )

    if ( mvars.number_of_frames > 1 and mvars.number_of_frames == frame ):
        f_name = "( Average over " + str(mvars.number_of_frames) + " Models )"
        corr_of_frame = mvars.ave_corr_coef[0][1]
        r_factor_of_frame = mvars.ave_r_factor
    else:
        f_name = " ( Model number = " + str(frame + 1) + " )"
        corr_of_frame = corr_trj[frame]
        r_factor_of_frame = rfactor_trj[frame]

    s1_title = "Experimental PCS vs. Residue Number " + f_name  
    s2_title = "Experimental PCS - Back-calculated PCS " + f_name
    #s3_title = "Corr.coeff = " + str('%.3f'%(mvars.ave_corr_coef[0][1])) + " , R_factor = " + str('%.3f'%(mvars.ave_r_factor)) + " (taken from normalized PCSs)" + f_name
    s3_title = "Corr.coeff = " + str('%.3f'%(corr_of_frame)) + " , R-factor = " + str('%.3f'%(r_factor_of_frame))  + f_name

    s4_title = "Distances between Hydrogens and Spin Label"

#    s1 = figure(plot_height=height, title=s1_title, tools=Tools, active_scroll="wheel_zoom", active_inspect="hover")
#    s1 = figure(plot_height=height, plot_width=width, title=s1_title, tools=[hover])
#    s2 = figure(plot_height=height, plot_width=width, title=s2_title, tools=Tools, active_scroll="wheel_zoom", active_inspect="hover")
#    s3 = figure(plot_height=height, plot_width=width, title=s3_title, tools=[hover,Tools])

    s1 = figure(plot_height=height, plot_width=width, title=s1_title, tools=[hover1,Tools])
    s2 = figure(plot_height=height, plot_width=width, title=s2_title, tools=[hover2,Tools])
    s3 = figure(plot_height=height, plot_width=width, title=s3_title, tools=[hover3,Tools])
    s4 = figure(plot_height=height, plot_width=width, title=s4_title, tools=[hover4,Tools])

    for i in range(ntype):
        source1[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = rdc_exp_all[i],
        ))
        
        s1.vbar(x='x', top = 'y', width = 0.5, color = colors[i], line_color = "black", legend = type_id[i], alpha = 0.5,source=source1[i])
    s1.legend.orientation = "horizontal"
    s1.legend.location = "top_left"
    s1.legend.click_policy="hide"
    s1.yaxis.axis_label = "PCS [ppm]"
    s1.xaxis.axis_label = "Residue Number"

    for i in range(ntype):
        source2[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = delta_rdc_all[i],
        ))
        s2.vbar(x='x', top = 'y', width = 0.5, color = colors[i], line_color = "black", legend = type_id[i], alpha = 0.5, source=source2[i])

    s2.legend.orientation = "horizontal"
    s2.legend.location = "top_left"
    s2.legend.click_policy="hide"
    s2.yaxis.axis_label = "PCS_exp - PCS_calc [ppm]"
    s2.xaxis.axis_label = "Residue Number"

    concat_exp = []
    concat_calc = []
    for i in range(ntype):
        concat_exp.extend(rdc_exp_all[i]) 
        concat_calc.extend(rdc_calc_all[i])

    concat_exp_mat = np.vstack([concat_exp, np.ones(len(concat_exp))]).T
    #slope, intercept = np.linalg.lstsq(concat_exp_mat,concat_calc,rcond=None)[0]
    slope, intercept = np.linalg.lstsq(concat_exp_mat,concat_calc,rcond=-1)[0] 

    x_range = min(concat_exp)*1.1,max(concat_exp)*1.1
    y_range = min(concat_calc)*1.1,max(concat_calc)*1.7

    s3.x_range = Range1d(*x_range)
    s3.y_range = Range1d(*y_range)


    for i in range(ntype):
        txt_title = [type_id[i] for va in range(len(chain_all[i]))]
        source3[i] = ColumnDataSource(data=dict(
            title = txt_title,
            chainid = chain_all[i],
            resid = res_all[i],
            x = rdc_exp_all[i],
            y = rdc_calc_all[i],
        ))
        s3.circle(x = 'x', y = 'y', legend = type_id[i], size = 8, color = colors[i], alpha = 0.8,source=source3[i])
#        s3.circle(x = rdc_exp_all[i], y = rdc_calc_all[i], size = 10, color = "red", alpha = 0.8,source=source3[i])
    #prepare linear regression plot

    n_model = 4

    x_model = np.linspace(min(concat_exp)*1.5, max(concat_exp)*1.5, num=n_model)
    y_model = np.add( np.multiply(x_model,slope), intercept)

    txt_blank = ["N/A" for va in range(len(y_model))]
    txt_title = ["Fit line" for va in range(len(y_model))]
    source3_fit = ColumnDataSource(data=dict(
        title = txt_title,
        chainid = txt_blank,
        resid = txt_blank,
        x = x_model,
        y = y_model,
    ))

#    s3.line(x_model, y_model, color="red", line_width = 3)
    s3.line(x = 'x', y='y', color="red", line_width = 3, source=source3_fit)
    s3.legend.orientation = "horizontal"
    s3.legend.location = "top_left"
    s3.legend.click_policy="hide"
    s3.yaxis.axis_label = "Back-calculated PCS [ppm]"
    s3.xaxis.axis_label = "Experimental PCS [ppm]"

    for i in range(ntype):
        source4[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = mvars.sl_distance,
        ))
        s4.vbar(x='x', top = 'y', width = 0.5, color = colors[i], line_color = "black", legend = type_id[i], alpha = 0.5, source=source4[i])

    y_range = 0, np.max(mvars.sl_distance)*1.2
    s4.y_range = Range1d(*y_range)

    s4.legend.orientation = "horizontal"
    s4.legend.location = "top_left"
    s4.legend.click_policy="hide"
    s4.yaxis.axis_label = "Distance from SL (Angstrom)"
    s4.xaxis.axis_label = "Residue Number"

    result = layout([[s1, s2],[s3, s4]] )

    save(result)
#    show(result)
#   For embedding
    script, div = components(result)
#    print ("##### This is Script for single model plot#####")
#    print (script)
#    print ("##### This is New Script for single model plot#####")
#    print (div)
    return script, div

def single_model_link(self, frame, height, width, app):

    mvars = self.mvars

    html_dir = os.path.join(mvars.runname, app, 'html') 

    if not os.path.exists(html_dir ):
        os.makedirs(html_dir)

    if (frame == mvars.number_of_frames):
        html_name = os.path.join(html_dir, 'Average.html')
    else:
        html_name = os.path.join(html_dir, str(frame + 1) + '.html')

    wtitle = "<h2>ALTENS results ( Model nubmber = " + str(frame + 1) + " )</h2>"
    if ( frame == mvars.number_of_frames):
        wtitle = "<h2>ALTENS results averaged over " + str(frame) + " models </h2>"
        if ( mvars.use_monte_carlo_flag ):
            wtitle = "<h2>ALTENS results averaged over " + str(frame) + " models (MC results not shown)</h2>"

#    output_file(html_name, title=wtitle )
    output_file(html_name)
    colors = ["blue", "green", "cyon", "black", "yellow", "violet", "red"]

    type_id, chain_all, res_all, rdc_exp_all, rdc_calc_all, delta_rdc_all = read_calc_pcs(self,frame,app)
    corr_trj,rfactor_trj = read_data_per_model(self, app)

    ntype = len(type_id)

#  Initializing configuration of plotting

    source1 = ["" for i in range(ntype)]
    source2 = ["" for i in range(ntype)]
    source3 = ["" for i in range(ntype)]

    for i in range(ntype):
        source1[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = rdc_exp_all[i],
        ))

    hover1 = HoverTool(
            tooltips = [
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("PCS_exp", "@y"),
            ]
            )
    hover2 = HoverTool(
            tooltips = [
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("Delta", "@y"),
            ]
            )
    hover3 = HoverTool(
            tooltips = [
            ("Type", "@title"),
            ("Chain", "@chainid"),
            ("Resid", '@resid'),
            ("PCS_exp", "@x"),
            ("PCS_recalc", "@y"),
            ]
            )

    if ( mvars.number_of_frames > 1 and mvars.number_of_frames == frame ):
        f_name = "( Average over " + str(mvars.number_of_frames) + " Models )"
        corr_of_frame = mvars.ave_corr_coef[0][1]
        r_factor_of_frame = mvars.ave_r_factor
    else:
        f_name = " ( Model number = " + str(frame + 1) + " )"
        corr_of_frame = corr_trj[frame]
        r_factor_of_frame = rfactor_trj[frame]

    s1_title = "Experimental PCS vs. Residue Number " # + f_name
    s2_title = "Experimental PCS - Back-calculated PCS " # + f_name
    #s3_title = "Corr.coeff = " + str('%.3f'%(mvars.ave_corr_coef[0][1])) + " , R_factor = " + str('%.3f'%(mvars.ave_r_factor)) + " (taken from normalized PCSs)" + f_name
    s3_title = "Corr.coeff = " + str('%.3f'%(corr_of_frame)) + " , R-factor = " + str('%.3f'%(r_factor_of_frame)) # + f_name

    s1 = figure(plot_height=height, plot_width=width, title=s1_title, tools=[hover1,Tools])
    s2 = figure(plot_height=height, plot_width=width, title=s2_title, tools=[hover2,Tools])
    s3 = figure(plot_height=height, plot_width=width, title=s3_title, tools=[hover3,Tools])

    for i in range(ntype):
        source1[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = rdc_exp_all[i],
        ))

        s1.vbar(x='x', top = 'y', width = 0.5, color = colors[i], line_color = "black", legend = type_id[i], alpha = 0.5,source=source1[i])
    s1.legend.orientation = "horizontal"
    s1.legend.location = "top_left"
    s1.legend.click_policy="hide"
    s1.yaxis.axis_label = "PCS [ppm]"
    s1.xaxis.axis_label = "Residue Number"

    for i in range(ntype):
        source2[i] = ColumnDataSource(data=dict(
            chainid = chain_all[i],
            resid = res_all[i],
            x = res_all[i],
            y = delta_rdc_all[i],
        ))
        s2.vbar(x='x', top = 'y', width = 0.5, color = colors[i], line_color = "black", legend = type_id[i], alpha = 0.5, source=source2[i])

    s2.legend.orientation = "horizontal"
    s2.legend.location = "top_left"
    s2.legend.click_policy="hide"
    s2.yaxis.axis_label = "PCS_exp - PCS_calc [ppm]"
    s2.xaxis.axis_label = "Residue Number"

    concat_exp = []
    concat_calc = []
    for i in range(ntype):
        concat_exp.extend(rdc_exp_all[i])
        concat_calc.extend(rdc_calc_all[i])

    concat_exp_mat = np.vstack([concat_exp, np.ones(len(concat_exp))]).T
    #slope, intercept = np.linalg.lstsq(concat_exp_mat,concat_calc,rcond=None)[0]
    slope, intercept = np.linalg.lstsq(concat_exp_mat,concat_calc,rcond=-1)[0]

    x_range = min(concat_exp)*1.1,max(concat_exp)*1.1
    y_range = min(concat_calc)*1.1,max(concat_calc)*1.7

    s3.x_range = Range1d(*x_range)
    s3.y_range = Range1d(*y_range)


    for i in range(ntype):
        txt_title = [type_id[i] for va in range(len(chain_all[i]))]
        source3[i] = ColumnDataSource(data=dict(
            title = txt_title,
            chainid = chain_all[i],
            resid = res_all[i],
            x = rdc_exp_all[i],
            y = rdc_calc_all[i],
        ))
        s3.circle(x = 'x', y = 'y', legend = type_id[i], size = 8, color = colors[i], alpha = 0.8,source=source3[i])
#        s3.circle(x = rdc_exp_all[i], y = rdc_calc_all[i], size = 10, color = "red", alpha = 0.8,source=source3[i])
    #prepare linear regression plot

    n_model = 4

    x_model = np.linspace(min(concat_exp)*1.5, max(concat_exp)*1.5, num=n_model)
    y_model = np.add( np.multiply(x_model,slope), intercept)

    txt_blank = ["N/A" for va in range(len(y_model))]
    txt_title = ["Fit line" for va in range(len(y_model))]
    source3_fit = ColumnDataSource(data=dict(
        title = txt_title,
        chainid = txt_blank,
        resid = txt_blank,
        x = x_model,
        y = y_model,
    ))

#    s3.line(x_model, y_model, color="red", line_width = 3)
    s3.line(x = 'x', y='y', color="red", line_width = 3, source=source3_fit)
    s3.legend.orientation = "horizontal"
    s3.legend.location = "top_left"
    s3.legend.click_policy="hide"
    s3.yaxis.axis_label = "Back-calculated PCS [ppm]"
    s3.xaxis.axis_label = "Experimental PCS [ppm]"

#    result = layout([[s1, s2, s3]])
    result = [[ s1, s2, s3 ]]

    if mvars.use_monte_carlo_flag:
        if ( frame != mvars.number_of_frames ):
            mc_result = generate_mc_histogram_link(self,frame,height, width,app)
            result.extend(mc_result)
             
    save( column(Div(text=wtitle, width=800),layout( result ) ) )

#    show(result)
#   For embedding
#    script, div = components(result)
#    print ("##### This is Script for single model plot#####")
#    print (script)
#    print ("##### This is New Script for single model plot#####")
#    print (div)
    return 

