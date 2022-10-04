from bokeh.layouts import layout, row, column
from bokeh.plotting import figure, output_file, show, save
from bokeh.models import HoverTool, Span, Range1d, ColumnDataSource,FactorRange , tickers, SingleIntervalTicker, LinearAxis, WheelZoomTool, LabelSet, Label, OpenURL, TapTool, Div, LayoutDOM, CustomJS
from bokeh.transform import factor_cmap
from bokeh.embed import components

import string, locale, re, os, json
import numpy as np

import sasmol.sasmol as sasmol

Tools = 'pan,box_zoom,wheel_zoom,save,reset,crosshair'
Tools_mc = 'pan,box_zoom,wheel_zoom,save,reset,crosshair,hover'
Tools_summary = 'pan,box_zoom,wheel_zoom,save,reset,crosshair, tap'

def read_calc_rdc(self,frame,app):
    
#   Only works for the first chain

    mvars = self.mvars
    calcfile=open(mvars.file_calcrdc[frame],'r').readlines()

    nchain_limit = 1
    chain_all = []
    res_all = []
    rdc_exp_all = []
    rdc_calc_all = []
    delta_rdc_all = []

    type_id = []
    prev_type = '' 
    prev_chain = ''

    nchain = 0
    ntype = 0

    for line in calcfile:
        txt = string.split(line)
        pat1 = re.compile('\s+|\h+')
        txttrim = pat1.sub('', ''.join(txt)) 

        if ( re.match(r"\d{1,2}", txttrim) ):
            vtype = txt[2]
            chain = txt[3]

            if ( prev_type != vtype ): 
                prev_type = vtype
                type_id.append(vtype)

                if (ntype == 0):
                    chainname = []  
                    resnumber = []
                    rdc_exp = []
                    rdc_calc = []
                    delta_rdc = []
                    rdc_exp = [] 
                else:
                    chain_all.append(chainname)
                    res_all.append(resnumber)
                    rdc_exp_all.append(rdc_exp)
                    rdc_calc_all.append(rdc_calc)
                    delta_rdc_all.append(delta_rdc)

                    chainname = []
                    resnumber = []
                    rdc_exp = []
                    rdc_calc = []
                    delta_rdc = []
                    rdc_exp = []

                ntype += 1

            if ( prev_chain != chain ):
                prev_chain = chain
                nchain += 1
                if (nchain > nchain_limit + 1):
                    break
            chainname.append(chain)
            resnumber.append(int(locale.atof(txt[4])))
            rdc_exp.append(float(locale.atof(txt[5])))
            rdc_calc.append(float(locale.atof(txt[6])))
            delta_rdc.append(float(locale.atof(txt[7])))

    chain_all.append(chainname)
    res_all.append(resnumber)
    rdc_exp_all.append(rdc_exp)
    rdc_calc_all.append(rdc_calc)
    delta_rdc_all.append(delta_rdc)

    return type_id, chain_all, res_all, rdc_exp_all, rdc_calc_all, delta_rdc_all 

def read_mc_trajectory(self,frame,app):
    
    mvars = self.mvars
    mcfile = open(mvars.file_mc_trajectory[frame],'r').readlines()
        
    nline = 0

    for line in mcfile:
        txt = string.split(line)
        if ( nline == 0 ) : 
            item = [x for x in txt]
            nitem = len(item)
            data = [[] for x in range(nitem)]
        else:
            for i in range(nitem):
                data[i].append(float(locale.atof(txt[i])))
        nline += 1

    return item, data

def read_data_per_model(self,app):

    mvars = self.mvars
    datafile = open(mvars.file_out_summary,'r').readlines()

    nline = 0

    for line in datafile:
        txt = string.split(line)
        if ( nline == 0 ) :
            corr_trj = []
            rfactor_trj = []
        else:
            corr_trj.append(float(locale.atof(txt[1])))
            rfactor_trj.append(float(locale.atof(txt[2])))
        nline += 1

    return corr_trj, rfactor_trj


def generate_mc_histogram(self,frame,height, width,app):

    mvars = self.mvars
    hist_legend, hist_data = read_mc_trajectory(self,frame,app)

    output_dir = os.path.join(mvars.runname, app)
    output_html = os.path.join(output_dir, mvars.runname + '_mc_histogram.html')

    wtitle = "MC analysis ( model = " + str(frame + 1) + " )"
    output_file(output_html, title=wtitle )

    nitems = len(hist_legend)

    ndata = len(hist_data)

    hist = ['' for x in range(nitems)]
    hist_count = [[] for x in range(nitems)]
    hist_edge = [[] for x in range(nitems)]
    hist_bin = [[] for x in range(nitems)]
    nbins = 31


#    hist_title = []
#    hist_title.append("X Component")
#    hist_title.append("Y Component")
#    hist_title.append("Z Component")
#    hist_title.append("Alpha")
#    hist_title.append("Beta")
#    hist_title.append("Gamma")

    for i in range(nitems):
        hist_count[i], hist_edge[i] = np.histogram(hist_data[i],bins=nbins)
        for j in range(nbins):
            hist_bin[i].append( 0.5*(hist_edge[i][j+1] + hist_edge[i][j]) )
        hist[i] = figure(plot_height=height, plot_width=width, tools=Tools_mc, active_scroll="wheel_zoom", active_inspect="hover")
        #hist[i] = figure(plot_height=height, plot_width=width, title = hist_legend[i], tools=Tools, active_scroll="wheel_zoom", active_inspect="hover")
        hist_width = 0.4*(hist_bin[i][1] - hist_bin[i][0])
#        hist[i].vbar(x=hist_bin[i], top=hist_count[i], color = "red", line_color = "red", legend = hist_title[i], width = hist_width)
        hist[i].vbar(x=hist_bin[i], top=hist_count[i], color = "red", line_color = "red", width = hist_width)

#        hist[i].legend.orientation = "horizontal" 
#        hist[i].legend.location = "top_left"
#        hist[i].legend.click_policy = "hide"
        if ( i == 3 or i == 5):
            x_range = -5, 185
            hist[i].x_range = Range1d(*x_range)
        if ( i == 4 ):
            x_range = -5, 185
            hist[i].x_range = Range1d(*x_range)
        hist[i].xaxis.axis_label = hist_legend[i].replace("["," [")
        hist[i].yaxis.axis_label = "Counts"
    hist[0].title.text = " MC Analysis for Alignment Tensor Eigenvalues (Model Number = 1) "
    hist[0].title.text_font_size = "14px"
    hist[3].title.text = " MC Analysis for Euler Angles (Model Number = 1) "
    hist[3].title.text_font_size = "14px"
  
    result = layout([ [ hist[0],hist[1],hist[2] ], [ hist[3], hist[4], hist[5] ] ]  )
    save(result) 

#   For embedding 
    script, div = components(result)
    return script, div

def generate_mc_histogram_link(self,frame,height, width,app):

    mvars = self.mvars
    hist_legend, hist_data = read_mc_trajectory(self,frame,app)

    nitems = len(hist_legend)

    ndata = len(hist_data)

    hist = ['' for x in range(nitems)]
    hist_count = [[] for x in range(nitems)]
    hist_edge = [[] for x in range(nitems)]
    hist_bin = [[] for x in range(nitems)]
    nbins = 31


    for i in range(nitems):
        hist_count[i], hist_edge[i] = np.histogram(hist_data[i],bins=nbins)
        for j in range(nbins):
            hist_bin[i].append( 0.5*(hist_edge[i][j+1] + hist_edge[i][j]) )
        hist[i] = figure(plot_height=height, plot_width=width, tools=Tools_mc, active_scroll="wheel_zoom", active_inspect="hover")
        #hist[i] = figure(plot_height=height, plot_width=width, title = hist_legend[i], tools=Tools, active_scroll="wheel_zoom", active_inspect="hover")
        hist_width = 0.4*(hist_bin[i][1] - hist_bin[i][0])
#        hist[i].vbar(x=hist_bin[i], top=hist_count[i], color = "red", line_color = "red", legend = hist_title[i], width = hist_width)
        hist[i].vbar(x=hist_bin[i], top=hist_count[i], color = "red", line_color = "red", width = hist_width)

#        hist[i].legend.orientation = "horizontal" 
#        hist[i].legend.location = "top_left"
#        hist[i].legend.click_policy = "hide"
        if ( i == 3 or i == 5):
            x_range = -5, 185
            hist[i].x_range = Range1d(*x_range)
        if ( i == 4 ):
            x_range = -5, 185
            hist[i].x_range = Range1d(*x_range)
        hist[i].xaxis.axis_label = hist_legend[i].replace("["," [")
        hist[i].yaxis.axis_label = "Counts"
    hist[0].title.text = " MC Analysis for Alignment Tensor Eigenvalues" 
    hist[0].title.text_font_size = "14px"
    hist[3].title.text = " MC Analysis for Euler Angles"
    hist[3].title.text_font_size = "14px"

    result = [ [ hist[0],hist[1],hist[2] ], [ hist[3], hist[4], hist[5] ] ] 

    return result 

def read_out_summary(self,app):
    '''
   rtype item: list for legends of plot
   rtype data: list of model, correlation coeff, r factor, three eigenvalues, three euler angles  

    '''
    mvars = self.mvars
    datafile = open(mvars.file_out_summary,'r').readlines()

    nline = 0

    for line in datafile:
        txt = string.split(line)
        if ( nline == 0 ) :
            item = [x for x in txt]
            nitem = len(item)
            data = [[] for x in range(nitem)]
        else:
            for i in range(nitem):
                data[i].append(float(locale.atof(txt[i])))
        nline += 1

    return item, data

def multi_model_plot(self, height, width, app):

    mvars = self.mvars

    pgui = self.run_utils.print_gui

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

    margin = 0.05
    x_range = 0.7, 1.05 
    y_range = -0.05, 0.35
    x_range = xmin - (xmax-xmin)*margin, xmax + (xmax-xmin)*margin
    y_range = ymin - (ymax-ymin)*margin, ymax + (ymax-ymin)*margin
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

    new_url =  '/'.join(path_array_slice) + '/' + output_dir + '/html/' + '@modelid' + '.html'
    url = new_url

#
#    taptool1 = s1.select(type=TapTool)
#    taptool1.callback = CustomJS(args=dict(source=source1), code="""
#        var mydiv = document.getElementById("link");
#        mydiv.innerHTML = "<a href='""" + url + """'>link</a>";
#        """)

    taptool1 = s2.select(type=TapTool)
    taptool1.callback = OpenURL(url=url)
    taptool2 = s2.select(type=TapTool)
    taptool2.callback = OpenURL(url=url)
    taptool3 = s3.select(type=TapTool)
    taptool3.callback = OpenURL(url=url)

    result = layout( [[s1, s2, s3]] )
    print ( "#### ALL FINE!")

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
    elif (frame == mvars.number_of_frames):
        output_html = os.path.join(output_dir, mvars.runname + '_single_model_average.html')
    else:
        output_html = os.path.join(output_dir, mvars.runname + '_single_model.html')    


    wtitle = "ALTENS results ( Model nubmber = " + str(frame + 1) + " )" 
    if ( frame == mvars.number_of_frames):
        wtitle = "ALTENS results averaged over " + str(frame) + " models"

    output_file(output_html, title=wtitle )

    colors = ["blue", "green", "cyon", "black", "yellow", "violet", "red"]

    type_id, chain_all, res_all, rdc_exp_all, rdc_calc_all, delta_rdc_all = read_calc_rdc(self,frame,app)
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
            ("RDC_exp", "@y"),
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
            ("RDC_exp", "@x"),
            ("RDC_recalc", "@y"),
            ]
            )

    if ( mvars.number_of_frames > 1 and mvars.number_of_frames == frame ):
        f_name = "( from average Back-calculated RDC over " + str(mvars.number_of_frames) + " Models )"
        corr_of_frame = mvars.ave_corr_coef[0][1]
        r_factor_of_frame = mvars.ave_r_factor
    else:
        f_name = " ( Model number = " + str(frame + 1) + " )"
        corr_of_frame = corr_trj[frame]
        r_factor_of_frame = rfactor_trj[frame]

    s1_title = "Experimental RDC vs. Residue Number " + f_name  
    s2_title = "Experimental RDC - Back-calculated RDC " + f_name
    #s3_title = "Corr.coeff = " + str('%.3f'%(mvars.ave_corr_coef[0][1])) + " , R_factor = " + str('%.3f'%(mvars.ave_r_factor)) + " (taken from normalized RDCs)" + f_name
    s3_title = "Corr.coeff = " + str('%.3f'%(corr_of_frame)) + " , R-factor = " + str('%.3f'%(r_factor_of_frame))  + f_name

#    s1 = figure(plot_height=height, title=s1_title, tools=Tools, active_scroll="wheel_zoom", active_inspect="hover")
#    s1 = figure(plot_height=height, plot_width=width, title=s1_title, tools=[hover])
#    s2 = figure(plot_height=height, plot_width=width, title=s2_title, tools=Tools, active_scroll="wheel_zoom", active_inspect="hover")
#    s3 = figure(plot_height=height, plot_width=width, title=s3_title, tools=[hover,Tools])

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
    s1.yaxis.axis_label = "RDC [Hz]"
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
    s2.yaxis.axis_label = "RDC_exp - RDC_calc [Hz]"
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
#    s3.line(x = 'x', y='y', color="red", line_width = 3, source=source3_fit)
    s3.line(x = 'x', y='x', color="red", line_width = 3, source=source3_fit)

    s3.legend.orientation = "horizontal"
    s3.legend.location = "top_left"
    s3.legend.click_policy="hide"
    s3.yaxis.axis_label = "Back-calculated RDC [Hz]"
    s3.xaxis.axis_label = "Experimental RDC [Hz]"

    result = layout([[s1, s2, s3]] )

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
        wtitle = "<h2>ALTENS results averaged over " + str(frame) + " models</h2>"

#    output_file(html_name, title=wtitle )
    output_file(html_name)
    colors = ["blue", "green", "cyon", "black", "yellow", "violet", "red"]

    type_id, chain_all, res_all, rdc_exp_all, rdc_calc_all, delta_rdc_all = read_calc_rdc(self,frame,app)
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
            ("RDC_exp", "@y"),
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
            ("RDC_exp", "@x"),
            ("RDC_recalc", "@y"),
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

    s1_title = "Experimental RDC vs. Residue Number " # + f_name
    s2_title = "Experimental RDC - Back-calculated RDC " # + f_name
    #s3_title = "Corr.coeff = " + str('%.3f'%(mvars.ave_corr_coef[0][1])) + " , R_factor = " + str('%.3f'%(mvars.ave_r_factor)) + " (taken from normalized RDCs)" + f_name
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
    s1.yaxis.axis_label = "RDC [Hz]"
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
    s2.yaxis.axis_label = "RDC_exp - RDC_calc [Hz]"
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
    s3.yaxis.axis_label = "Back-calculated RDC [Hz]"
    s3.xaxis.axis_label = "Experimental RDC [Hz]"

#    result = layout([[s1, s2, s3]])
    result = [[ s1, s2, s3 ]]

    if mvars.use_monte_carlo_flag:
        if ( frame != mvars.number_of_frames ):
            mc_result = generate_mc_histogram_link(self,frame,height, width,app)
            result.extend(mc_result)
             
    save( column(Div(text=wtitle, width=600),layout( result ) ) )

#    show(result)
#   For embedding
#    script, div = components(result)
#    print ("##### This is Script for single model plot#####")
#    print (script)
#    print ("##### This is New Script for single model plot#####")
#    print (div)
    return 

