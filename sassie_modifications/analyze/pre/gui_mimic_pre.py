'''
Driver method to run the Altens module
'''

import sys

sys.path.append('./')
import pre

#import sassie.analyze.altens.altens as altens
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####

#### test multi-frame 
runname = 'run_0'
ratiofile = "./ratio_4_SLfit_test_Q49.txt"
pdbfile='./1D3Z_f.pdb'
dcdfile='./1D3Z_f.pdb'

T2dia = "0.05"
Htime = "0.005"
freq = "600.13"
TAUc = "4.5"

#### end user input ####
#### end user input ####
#### end user input ####

svariables={}

svariables['runname'] = (runname,'string')
svariables['ratiofile'] = (ratiofile,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['dcdfile'] = (dcdfile,'string')
svariables['T2dia'] = (T2dia,'float')
svariables['Htime'] = (Htime,'float')
svariables['freq'] = (freq,'float')
svariables['TAUc'] = (TAUc,'float')

error, variables = input_filter.type_check_and_convert(svariables)

if len(error) > 0:
    print 'error = ', error
    sys.exit()
else:
    pass

import time; start = time.time()
txtQueue = multiprocessing.JoinableQueue()

plotQueues = dict()

plotQueues['bokeh_plot_1'] = multiprocessing.JoinableQueue()
plotQueues['bokeh_plot_2'] = multiprocessing.JoinableQueue()
plotQueues['bokeh_plot_3'] = multiprocessing.JoinableQueue()
plotQueues['bokeh_plot_4'] = multiprocessing.JoinableQueue()
#plotQueues['bokeh_plot_5'] = multiprocessing.JoinableQueue()

pre = pre.pre()
pre.main(variables, txtQueue, plotQueues)

this_text = txtQueue.get(True, timeout=0.1)

print ("time used: ",time.time()-start)

#print 'in GUI and txtOutput = ', this_text, '\n'


