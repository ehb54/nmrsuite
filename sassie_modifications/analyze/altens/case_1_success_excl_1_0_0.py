'''
Driver method to run the Altens module
'''

# Test case
# without exclusion list

import sys

sys.path.append('./')
import altens

#import sassie.analyze.altens.altens as altens
import sassie.interface.input_filter as input_filter
#import sassie.interface.altens_filter as altens_filter
import altens_filter as altens_filter
import multiprocessing

runname_base = "case_1_success_excl_1" 
runname_base = "test_ave_coef"
rdc_dir = "./data/input_rdc/three_column/correct/"
exclusion_dir = "./data/input_exclusion/correct/"

svariables = {}

#### user input ####

runname = runname_base + "_" + str(0) + "_" + str(0)
rdc_input_file = rdc_dir + str(0) + ".txt"
pdbfile='./data/1D3Z_mod1.pdb'
#dcdfile='./data/1D3Z_mod1.pdb'
dcdfile='./data/1D3Z_10frames.pdb'
use_monte_carlo_flag=True

residue_list_file_flag=True
residue_list_file = exclusion_dir + str(0) + ".txt"

number_of_monte_carlo_steps = '500'
seed = '1,123'

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

error = altens_filter.check_altens(variables)
if(len(error) > 0):
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

altens_main = altens.altens()
altens_main.main(variables, txtQueue, plotQueues)

this_text = txtQueue.get(True, timeout=0.1)

plot_text = plotQueues['bokeh_plot_1'].get(True, timeout=2.0)

print ("time used: ",time.time()-start)

print ("JOB: " + runname + " Done")
#print 'in GUI and txtOutput = ', this_text, '\n'


