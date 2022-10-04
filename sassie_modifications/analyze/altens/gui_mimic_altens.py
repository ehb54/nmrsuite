'''
Driver method to run the Altens module
'''

import sys

sys.path.append('./')
import altens

#import sassie.analyze.altens.altens as altens
import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####

#### test multi-frame 
runname = 'run_0'
rdc_input_file = "./data/input_rdc/twotypes/three_column/correct/1.txt"
pdbfile='./data/1D3Z_1frame.pdb'
dcdfile='./data/1D3Z_10frames.pdb'

residue_list_file_flag=False
residue_list_file = "reslist_Ub.txt"
use_monte_carlo_flag=True
number_of_monte_carlo_steps = '500' 
seed = '1,123'

#### end user input ####
#### end user input ####
#### end user input ####

svariables={}

svariables['runname'] = (runname,'string')
svariables['rdc_input_file'] = (rdc_input_file,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['dcdfile'] = (dcdfile,'string')
svariables['residue_list_file_flag'] =(residue_list_file_flag,'boolean')
svariables['residue_list_file'] = (residue_list_file,'string')
svariables['use_monte_carlo_flag'] =(use_monte_carlo_flag,'boolean')
svariables['number_of_monte_carlo_steps'] = (number_of_monte_carlo_steps,'int')
svariables['seed'] = (seed, 'int_array')

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
plotQueues['bokeh_plot_5'] = multiprocessing.JoinableQueue()

altens = altens.altens()
altens.main(variables, txtQueue, plotQueues)

this_text = txtQueue.get(True, timeout=0.1)

print ("time used: ",time.time()-start)

#print 'in GUI and txtOutput = ', this_text, '\n'


