'''
Driver method to run the PCS module
'''

import sys

# For local testing
sys.path.append('./')
import pcs as pcs
import pcs_filter as pcs_filter

# For compiled version
#import sassie.analyze.pcs.pcs as pcs
#import sassie.interface.pcs_filter as pcs_filter  

import sassie.interface.input_filter as input_filter
import multiprocessing

svariables = {}

#### user input ####

#### test multi-frame 
runname = 'run_0'
pcs_input_file = "./data/PCS_data.dat"
pdbfile='./data/1D3Z_mod1.pdb'
dcdfile='./data/1D3Z_mod1.pdb'

residue_exclusion_file_flag=False
residue_exclusion_file = "reslist_Ub.txt"

cond_num_cutoff = "33"
tolerance = "1e-4"

#### end user input ####
#### end user input ####
#### end user input ####

svariables={}

svariables['runname'] = (runname,'string')
svariables['pcs_input_file'] = (pcs_input_file,'string')
svariables['pdbfile'] = (pdbfile,'string')
svariables['dcdfile'] = (dcdfile,'string')
svariables['residue_exclusion_file_flag'] =(residue_exclusion_file_flag,'boolean')
svariables['residue_exclusion_file'] = (residue_exclusion_file,'string')
svariables['cond_num_cutoff'] = (cond_num_cutoff,'int')
svariables['tolerance'] = (tolerance,'float')

error, variables = input_filter.type_check_and_convert(svariables)

if len(error) > 0:
    print 'error = ', error
    sys.exit()
else:
    pass

error = pcs_filter.check_pcs(variables)
if(len(error) > 0):
    print 'error = ', error
    sys.exit()
else:
    pass

import time; start = time.time()
txtQueue = multiprocessing.JoinableQueue()

plotQueues = dict()

#plotQueues['bokeh_plot_1'] = multiprocessing.JoinableQueue()
#plotQueues['bokeh_plot_2'] = multiprocessing.JoinableQueue()
#plotQueues['bokeh_plot_3'] = multiprocessing.JoinableQueue()
#plotQueues['bokeh_plot_4'] = multiprocessing.JoinableQueue()
#plotQueues['bokeh_plot_5'] = multiprocessing.JoinableQueue()

pcs = pcs.pcs()
pcs.main(variables, txtQueue, plotQueues)

this_text = txtQueue.get(True, timeout=0.1)

print ("time used: ",time.time()-start)

#print 'in GUI and txtOutput = ', this_text, '\n'


