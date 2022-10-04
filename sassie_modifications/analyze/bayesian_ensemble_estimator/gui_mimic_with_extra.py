#!/share/apps/local/anaconda2/bin/python

import sys,os
sys.path.append(".")
import numpy as np
import multiprocessing
import sassie.interface.input_filter as input_filter
import ensemble_fit
from mpi4py import MPI
svariables = {}

#### MPI Environment ####
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

#### user input ####

#Standard Input
runname = 'withExtra'
sas_data = 'saxs.dat'
theoretical_profiles_zip = './TheoryFiles.zip' 
acceptance_tuner = '0.1' #Larger values = looser selection.  I.e., very large = random walk, very small = gradient minimization
posterior_burn = '1000' #Number of iterations per walker to remove from beginning (lose initial conformation information
max_iterations = '10000' #Maximum number of MC steps per walker
number_of_MCs = '3'
nproc = '7' #Number of processors to run on
d_max = '88.89' #If d_max is non-zero, do Shannon sampling, otherwise do full space chi2

#Advanced Input
use_all_members = 'False' #Do the entire set, do not conduct iterative search for AIC-identified subset
calc_relative_likelihoods = 'True'
second_dimension_data = 'extra.dat'
use_bic = False
ic_tolerance = '0.0'
## end user input ##

#### Parse inputs
svariables['runname'] = (runname, 'string')
svariables['sas_data'] = (sas_data, 'string')
svariables['theoretical_profiles_zip'] = (theoretical_profiles_zip, 'string')
svariables['acceptance_tuner'] = (acceptance_tuner, 'float')
svariables['posterior_burn'] = (posterior_burn, 'int')
svariables['max_iterations'] = (max_iterations, 'int')
svariables['number_of_MCs'] = (number_of_MCs, 'int')
svariables['nproc'] = (nproc, 'int')
svariables['d_max'] = (d_max, 'float')
svariables['second_dimension_data'] = (second_dimension_data, 'string')
svariables['use_all_members'] = (use_all_members, 'boolean')
svariables['calc_relative_likelihoods'] = (calc_relative_likelihoods, 'boolean')
svariables['ic_tolerance'] = (ic_tolerance, 'float')
svariables['use_bic'] = (use_bic, 'boolean')
####

error,variables = input_filter.type_check_and_convert(svariables)

if len(error) > 0:
	print("error = ",error)
	sys.exit()

txtQueue = multiprocessing.JoinableQueue()
reweighting = ensemble_fit.ensemble_routine()
reweighting.main(variables, txtQueue)
this_text = txtQueue.get(True,timeout=1.0)
