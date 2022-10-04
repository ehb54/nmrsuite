'''
    ENSEMBLE_FIT is a module that will read in a list of scattering
    profiles, along with an optional second set of restraints, and it will
    reweight populations in subsets of these scattering profiles to find
    the best fit to the experimental data. Here, "best fit" means a balance
    between a model with a low chi^2 value that simultaneously uses as few
    fitting parameters (number of profiles) as possible, thereby trying to
    avoid overfitting. This module will output information regarding this
    "best mode", and it has the option to also present the relative
    ability of all other models in comparison. If you use this module,
    please be sure to cite the references below.

    REFERENCES:

    Bowerman et al.
    Journal of Chemical Theory and Computation, 13, pp 2418-2429 (2017)

    Bowerman et al.
    Placeholder for Release Paper
'''

import numpy as np
import os,sys
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import zipfile
import itertools
import time
import pickle
import subprocess
from mpi4py import MPI

if sasconfig.__level__ == 'DEBUG':
    DEBUG = True
else:
    DEBUG = False

app = 'bayesian_ensemble_fit'

#Location of the ensemble modeling runtime
executable= '/home/sbowerma/sassie_module/v1/ensemble_fit_parallel_routine.py'
executable='/share/apps/local/anaconda2/lib/python2.7/site-packages/sassie/analyze/bayesian_ensemble_fit/ensemble_fit_parallel_routine.py'
executable='/share/apps/local/anaconda/lib/python2.7/site-packages/sassie/analyze/bayesian_ensemble_fit/ensemble_fit_parallel_routine.py'
#Location of the proper mpiexec command
mpiexec   = '/share/apps/local/anaconda2/bin/mpirun'
mpiexec   = '/share/apps/local/anaconda/bin/mpirun'

class module_variables():
    def __init__(self, parent = None):
        self.app = app

class efvars():
    ''' _E_nsemble _F_it _VAR_iable_S_ '''
    def __init__(self, parent = None):
        pass

######## MPI Environment ########
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()
#################################


# This is the larger object for managing the Ensemble Fit Routine at top-level
'''
    This object builds up all necessary variables and facilitates the forking 
    of each sub-basis proc. After the sub-basis routines are complete, this 
    object also handles likelihood comparisons and printing of the 'best' model
    information.
'''
class ensemble_routine(object):
    def __init__(self, parent = None):
        pass

    def main(self,input_variables,txtOutput):
        '''
        main method to handle iterative Bayesian MC routine
        '''
        self.mvars = module_variables()
        self.efvars = efvars()
        self.run_utils = module_utilities.run_utils(app,txtOutput)
        self.run_utils.setup_logging(self)
        self.log.debug('In main()')
        self.UnpackVariables(input_variables)
        self.run_utils.general_setup(self)
        self.Initialize()
        self.PickleVars()
        self.EnsembleFit()
        self.Epilogue()

    def UnpackVariables(self,variables):
        '''
        extract variables from gui/mimic into system wide class instance
        '''
        log   = self.log
        mvars = self.mvars

        log.debug('in UnpackVariables()')

        #Standard Options
        mvars.runname                       = variables['runname'][0]
        mvars.sas_data                      = variables['sas_data'][0]
        mvars.theoretical_profiles_zip      = variables['theoretical_profiles_zip'][0]
        mvars.acceptance_tuner              = variables['acceptance_tuner'][0]
        mvars.posterior_burn                = variables['posterior_burn'][0]
        mvars.max_iterations                = variables['max_iterations'][0]
        mvars.number_of_MCs                 = variables['number_of_MCs'][0]
        mvars.nproc                         = variables['nproc'][0]
        mvars.d_max                         = variables['d_max'][0]

        #Advanced User Options
        mvars.use_all                       = variables['use_all_members'][0]
        mvars.second_dimension_data         = variables['second_dimension_data'][0]
        mvars.calc_relative_likelihoods     = variables['calc_relative_likelihoods'][0]
        mvars.use_bic                       = variables['use_bic'][0]
        mvars.ic_tolerance                  = variables['ic_tolerance'][0]
        return

    def PickleVars(self):
        log    = self.log
        mvars  = self.mvars
        efvars = self.efvars
        log.debug('In PickleVars.')

        self.mvarspickle = os.path.join(efvars.output_folder,mvars.runname+"_mvars.p")
        pickle.dump(mvars,open(self.mvarspickle,'wb'))
        log.debug('mvars have been written to pickle: '+self.mvarspickle)

        self.efvarspickle = os.path.join(efvars.output_folder,mvars.runname+'_efvars.p')
        pickle.dump(efvars,open(self.efvarspickle,'wb'))
        log.debug('efvars have been written to pickle: '+self.efvarspickle)
        
        return

    def Initialize(self):
        '''
        Prepare efvars
        '''
        log     = self.log
        mvars   = self.mvars
        efvars  = self.efvars
        efvars.debug = DEBUG
        log.debug('In Initialize.')
        
        efvars.output_folder = os.path.join(mvars.runname,app)
        if not os.path.isdir(efvars.output_folder):
            os.mkdir(efvars.output_folder)
        efvars.Q, efvars.I, efvars.ERR      = np.genfromtxt(mvars.sas_data,unpack=True,usecols=(0,1,2))
        efvars.aic_file = os.path.join(efvars.output_folder,mvars.runname\
                                       +'_iteration_tracker.dat')
        trackf = open(efvars.aic_file,'w')
        trackf.write('#SubSize\tAIC\tChi^2\n')
        trackf.close()
        efvars.figfile = os.path.join(efvars.output_folder,mvars.runname\
                                      +'_AIC_vs_iteration.pdf')
        if mvars.d_max > 0.00:
            log.debug('Using Shannon Sampling (chi^2_free).')
            efvars.do_shannon_sampling  = True
            self.BuildShannonSamples()
            log.debug('There are '+str(efvars.number_of_channels)+' Shannon channels.')
            efvars.samples_Q    = efvars.Q[efvars.shannon_samples]
            efvars.samples_I    = efvars.I[efvars.shannon_samples]
            efvars.samples_ERR  = efvars.ERR[efvars.shannon_samples]
            efvars.num_q        = efvars.number_of_channels
        else:
            log.debug('Using standard chi^2 (no shannon sampling).')
            efvars.do_shannon_sampling  = False
            efvars.samples_Q    = efvars.Q
            efvars.samples_I    = efvars.I
            efvars.samples_ERR  = efvars.ERR
            efvars.num_q        = len(efvars.Q)
        
        if mvars.second_dimension_data != '':
            log.debug('Including Secondary Data Set')
            efvars.include_second_dimension         = True
            efvars.extra_data, efvars.extra_error   = np.genfromtxt(mvars.second_dimension_data,usecols=(0,1),unpack=True,dtype=float)
            efvars.num_extra                        = len(efvars.extra_data)
        else:
            log.debug('Not using Secondary Data Set')
            efvars.include_second_dimension         = False
            efvars.num_extra                        = 0

        efvars.num_points   = efvars.num_q + efvars.num_extra
        self.UnpackTheoreticalProfiles()
        self.BuildBasis()
        efvars.aic_tracker = np.array([],dtype=float)
        efvars.chi_tracker = np.array([],dtype=float)
        return

    def BuildShannonSamples(self):
        '''
        Build Shannon Channels for Calculating Chi_Free
        '''
        log     = self.log
        mvars   = self.mvars
        efvars  = self.efvars

        log.debug('In BuildShannonSamples')
        
        shannon_width = np.pi/mvars.d_max
        n_channels = int(np.floor(efvars.Q[-1]/shannon_width))

        shannon_points = np.random.choice(np.where(efvars.Q < shannon_width)[0],3001)
        for idx in range(1,n_channels):
            shannon_points = np.vstack((shannon_points,np.random.choice(np.where(np.logical_and(efvars.Q >= idx * shannon_width, efvars.Q < (idx+1)*shannon_width))[0],3001)))
        
        efvars.shannon_samples      = np.copy(shannon_points)
        efvars.number_of_channels   = n_channels
        
        return

    def UnpackTheoreticalProfiles(self):
        mvars   = self.mvars
        efvars  = self.efvars
        log     = self.log

        log.debug('In UnpackTheoreticalProfiles')
        
        efvars.profiles_directory = os.path.join(efvars.output_folder,'TheoreticalSASProfiles')
        log.debug('efvars.profiles_directory - '+efvars.profiles_directory)
        archive = zipfile.ZipFile(mvars.theoretical_profiles_zip)
        archive.extractall(efvars.profiles_directory)
        archive.close()

        return

    def BuildBasis(self):
        '''
        Stores individual theoretical profiles into single arrays for faster ensemble averaging
        '''
        log     = self.log
        mvars   = self.mvars
        efvars  = self.efvars

        log.debug('In BuildScatterBasis')

        efvars.name_array = np.genfromtxt(os.path.join(efvars.profiles_directory,'filelist.txt'),usecols=0,dtype=str)
        efvars.number_of_profiles = len(efvars.name_array)
        #Below is used to quickly build subsets, even if it seems redundant
        efvars.id_array = np.arange(efvars.number_of_profiles,dtype=int)
        #Dictionary linking subspace id lists to subspace simulation objects
        efvars.subspace_dict = {}
        efvars.full_scattering_basis = np.zeros((efvars.number_of_profiles,len(efvars.Q)),dtype=float)

        if efvars.include_second_dimension:
            try:
                efvars.extra_name_array = np.genfromtxt(os.path.join(efvars.profiles_directory,'filelist.txt'),usecols=1,dtype=str)
            except:
                log.error('In order to use a second dataset for fitting,'\
                          + ' users must include the list of files in the'\
                          + ' \"filelist.txt\" file.') 
            efvars.number_of_extra_profiles = len(efvars.extra_name_array)
            if efvars.number_of_extra_profiles != efvars.number_of_profiles:
                log.error('A different number of extra profiles ('+str(efvars.number_of_extra_profiles)+') have been supplied than the number of scattering profiles ('+str(efvars.number_of_profiles)+').')
            efvars.full_extra_basis = np.zeros((efvars.number_of_profiles,efvars.num_extra),dtype=float)
        else:
            efvars.full_extra_basis = None

        for idx in efvars.id_array:
            scatter_file = os.path.join(efvars.profiles_directory,efvars.name_array[idx])
            log.debug('Scatter File - '+str(scatter_file))
            the_q, the_i = np.genfromtxt(scatter_file,dtype=float,usecols=(0,1),unpack=True)
            if not np.all(the_q==efvars.Q):
                log.error('Q-values of '+scatter_file+' don\'t match the experimental values.\n\tthe_q = '+str(the_q)+'\n\tefvars.Q = '+str(efvars.Q))
            efvars.full_scattering_basis[idx]   = the_i

            if efvars.include_second_dimension:
                extra_file = os.path.join(efvars.profiles_directory,efvars.extra_name_array[idx])
                log.debug('Extra File - '+str(extra_file))
                extra_data = np.genfromtxt(extra_file,usecols=0,dtype=float)
                if len(extra_data) != efvars.num_extra:
                    log.error('The number of points in '+extra_file+' ('+str(len(extra_data))+') doesn\'t match the experimental values ('+str(efvars.num_extra)+').')
                efvars.full_extra_basis[idx]    = extra_data

        log.debug('Scattering Basis - ')
        log.debug(efvars.full_scattering_basis)
        if efvars.include_second_dimension:
            log.debug('Extra Basis - ')
            log.debug(efvars.full_extra_basis)

        return

    def EnsembleFit(self):
        '''
        This is the actual iterative routine
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        log     = self.log

        runcommand = mpiexec + ' -np ' + str(mvars.nproc)\
                   + ' '+ executable + ' -mpick '\
                   + self.mvarspickle + ' -efpick '\
                   + self.efvarspickle
        log.debug('Executable command: '+runcommand)
        p=subprocess.Popen(runcommand,shell=True,executable='/bin/bash')
        sts=os.waitpid(p.pid,0)[1]
        log.debug('Executable PID = '+str(p.pid))
        
        #self.Individuals()
        #self.best_model = self.FindBest()
        #if mvars.use_all == True:
        #    sim_object = simulated_basis(self,efvars.id_array)
        #    sim_object.BayesMC()
        #    self.best_model = self.FindBest()
        #else:
        #    for subsize in range(2,efvars.number_of_profiles + 1):
        #        log.info('Processing subsets of size '+str(subsize))
        #        if rank==0:
        #            sets = list(itertools.combinations(efvars.id_array,subsize))
        #            for thread in range(1,size):
        #                comm.send(sets,dest=thread)
        #                log.debug('Master thread sending subsets to Thread '+str(thread)+'.')
        #        all_sets_as_list = []
        #        Nsets = len(sets)
        #        for tracker in range(rank,Nsets,size):
        #            as_list = np.array([],dtype=int)
        #            for submember in sets[tracker]:
        #                as_list = np.append(as_list,submember)
        #            all_sets_as_list.append(str(as_list))
        #            #log.debug('all_sets_as_list: '+str(all_sets_as_list))
        #            efvars.subspace_dict[str(as_list)] = simulated_basis(self,as_list)
        #            efvars.subspace_dict[str(as_list)].BayesMC()
                
        #        if (rank!=0 and rank < np.min([size,Nsets])): #The second boolean is to protect from acase where there are fewer sets than requested threads.
        #            log.debug('Thread '+str(rank)+' is sending its subspace dictionary to Master.')
        #            comm.send(efvars.subspace_dict,dest=0)
        #            efvars.subspace_dict = comm.recv(source=0)
        #            log.debug('Thread '+str(rank)+' has received updated subspace dictionary from Master.')
        #        elif (rank >= Nsets):
        #            efvars.subspace_dict = comm.recv(source=0)
        #            log.debug('Thread '+str(rank)+' did not participate in this round of BayesMC simulations, but has received the updated subspace dictionary.')
        #        else:
        #            for thread in range(1,np.min([size,Nsets])):
        #                dum_subspace_dict = comm.recv(source=thread)
        #                efvars.subspace_dict.update(dum_subspace_dict)
        #            for thread in range(1,size):
        #                comm.send(efvars.subspace_dict,dest=thread)
        #            

        #        self.FindBest(sublist=all_sets_as_list)
        #        if self.subset_min_aic > self.min_aic:
        #            if rank==0: log.info('Thread '+str(rank)+': Best subset of size '+str(subsize)+' is a poorer fit than '+str(subsize-1)+' according to AIC.  Terminating Iterative Inclusion.') 
        #            break
        #        else:
        #            if rank==0: log.info('Thread '+str(rank)+': Best subset of size '+str(subsize)+' is an improvement over '+str(subsize-1)+' according to AIC.  Expanding subset size.')
        #            self.min_aic    = self.subset_min_aic
        #            self.best_model = self.subset_best_model
        #if rank==0: log.info('Best model is: '+str(self.best_model)+', AIC = '+str(efvars.subspace_dict[str(self.best_model)].aic)+', saxs_chi2 = '+str(efvars.subspace_dict[str(self.best_model)].saxs_chi2)+', extra_chi2 = '+str(efvars.subspace_dict[str(self.best_model)].extra_chi2))
        #if mvars.calc_relative_likelihoods:
        #    self.RelativeLikelihoods()

    def Individuals(self):
        '''
        This will find the fitting capabilities of each individual basis member.
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        log     = self.log

        log.debug('Thread '+str(rank)+': In Individuals.')
        log.debug('Thread '+str(rank)+': There are '+str(efvars.number_of_profiles)+' profiles.')
        for idx in range(efvars.number_of_profiles):
            efvars.subspace_dict[str(idx)]=simulated_basis(self,[idx])
            efvars.subspace_dict[str(idx)].BayesMC()

    def FindBest(self,sublist=np.array([],dtype=str)):
        '''
        Identifies the best model from the subset, and its relative likelihood to the next best fit.
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        log     = self.log
        if rank == 0:
            tic = time.time()
            log.debug('Thread '+str(rank)+': In FindBest.')
            self.min_aic = 999999999.9
            if len(sublist)==0:
                log.debug('No subset provided, comparing ALL sub-basis')
                for key in efvars.subspace_dict:
                    key_AIC = efvars.subspace_dict[key].aic
                    if key_AIC < self.min_aic:
                        self.min_aic    = key_AIC
                        self.best_model = key
                #self.min_aic = minAIC
                #return best_model
            else:
                self.subset_min_aic  = 9999999.9
                for key in sublist:
                    key_AIC = efvars.subspace_dict[str(key)].aic
                    log.debug(str(key_AIC)+','+str(self.subset_min_aic))
                    if key_AIC < self.subset_min_aic:
                        self.subset_min_aic     = key_AIC
                        self.subset_best_model  = key
                log.info('Thread '+str(rank)+': Best sub-basis  is: '+str(self.subset_best_model)+', AIC: '+str(efvars.subspace_dict[str(key)].aic))

    def RelativeLikelihoods(self):
        mvars   = self.mvars
        efvars  = self.efvars
        log     = self.log
        
        if rank==0:
            log.debug('Thread '+str(rank)+': Calculating Relative Likelihoods')

            outfile = open(os.path.join(efvars.output_folder,runname+'_relative_likelihoods.dat'),'w')
            outfile.write('#relative_likelihood,ID1,ID1_pop,ID2,ID2_pop,...,IDn,IDn_pop,normalized chi^2,AIC\n')
            for key in efvars.subspace_dict:
                log.debug('writing for key:'+str(key))
                model = efvars.subspace_dict[key]
                model.relative_likelihood = np.exp((self.min_aic-model.aic)/2.0)
                outfile.write(str(np.around(model.relative_likelihood,decimals=3))+',')
                if model.subset_size > 1:
                    for idx in range(model.subset_size):
                        outfile.write(str(efvars.name_array[model.subset_members[idx]])+','+str(np.around(model.total_mean_weights[idx],decimals=3))+',')
                    outfile.write(str(np.around(model.total_chi2/efvars.num_points,decimals=3))+','+str(np.around(model.aic,decimals=3))+'\n')
                else:
                    outfile.write(str(efvars.name_array[model.subset_members[0]])+','+str(1.000)+','+str(np.around(model.total_chi2/efvars.num_points,decimals=3))+','+str(np.around(model.aic,decimals=3))+'\n')
            outfile.close()


    def Epilogue(self):
        '''
        Saves files and prints info
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        log     = self.log
        log.debug('In Epilogue.')
#        if rank==0:
#            log.debug('Thread '+str(rank)+': In Epilogue')
#
#            print_model = efvars.subspace_dict[self.best_model]
#            outfile = open(os.path.join(efvars.output_folder,runname+'_final_model.dat'),'w')
#            outfile.write('#ScatterFile avg_weight var_weight lone_saxs_chi2 lone_extra_chi2 lone_total_chi2\n')
#            if print_model.subset_size > 1:
#                for idx in range(print_model.subset_size):
#                    outfile.write(efvars.name_array[print_model.subset_members[idx]]+' '+str(np.around(print_model.total_mean_weights[idx],decimals=3))+' '+str(np.around(efvars.subspace_dict[str(idx)].saxs_chi2,decimals=3))+' '+str(np.around(efvars.subspace_dict[str(idx)].extra_chi2,decimals=3))+' '+str(np.around(efvars.subspace_dict[str(idx)].total_chi2,decimals=3))+'\n')
#            else:
#                outfile.write(efvars.name_array[print_model.subset_members[0]]+'\n')
#            outfile.write('SAXS_chi^2  (model): '+str(np.around(print_model.saxs_chi2,decimals=3))+'\n')
#            outfile.write('Extra_chi^2 (model): '+str(np.around(print_model.extra_chi2,decimals=3))+'\n')
#            outfile.write('Total_chi^2 (model): '+str(np.around(print_model.total_chi2,decimals=3))+'\n')
#            outfile.write('AIC         (model): '+str(np.around(print_model.aic,decimals=3))+'\n')
#            outfile.close()
#
#            for thread in range(1,size):
#                log.debug('Master thread telling other threads that task is complete.')
#                comm.send(1,dest=thread)
#
#        else:
#            comm.recv(source=0)
#            log.debug('Thread '+str(rank)+' has been told that task is complete.')

########################  Subset MonteCarlo object class  ##########################################
''' 
    This object is the one that actually does the Monte Carlo for each subset.
    Relevant metrics (Posterior, Convergence, AIC/BIC, etc.) are stored herein.
'''
class simulated_basis(object):

    def __init__(self, parent, subset_ids):
        self.mvars          = parent.mvars
        self.efvars         = parent.efvars

        mvars   = self.mvars
        efvars  = self.efvars

        self.subset_members = subset_ids
        self.subset_basis   = efvars.full_scattering_basis[subset_ids]
        self.subset_size    = len(subset_ids)

        if efvars.include_second_dimension:
            self.subset_extra   = efvars.full_extra_basis[subset_ids]
        else:
            self.extra_chi2 = 0.0

    def BayesMC(self):
        '''
        This is the Monte Carlo routine.
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        self.logfile = os.path.join(efvars.output_folder,mvars.runname+"_"+str(subset_ids)+"basis_BayesMC.log")
        if os.path.isfile(self.logfile):
            os.system('rm '+self.logfile)
        #The posterior array has form [id_1,id_2,...,id_n,SAXSchi^2,Dimension2chi^2,Likelihood]
        if self.subset_size == 1:
            self.ensemble_I = np.copy(self.subset_basis[0])
            if efvars.include_second_dimension:
                self.ensemble_extra = np.copy(self.subset_extra[0])
            self.Chi2()
            self.max_likelihood = self.Likelihood()
            self.AIC()
        else:
            self.posterior_array=np.zeros((mvars.number_of_MCs,mvars.max_iterations,self.subset_size+3),dtype=float)
            for self.run in range(mvars.number_of_MCs):
                os.system('echo \"Beginning run '+str(self.run)+'\" >> '+self.logfile)
                self.weights    = np.random.random(size=self.subset_size)
                self.weights    = self.weights/float(np.sum(self.weights))
                self.ensemble_I = np.dot(self.weights,self.subset_basis)
                if efvars.include_second_dimension:
                    self.ensemble_extra = np.dot(self.weights,self.subset_extra)
                self.Chi2()
                for member in range(self.subset_size):
                    self.posterior_array[self.run,0,member]             = self.weights[member]
                self.posterior_array[self.run,0,self.subset_size]       = self.saxs_chi2
                self.posterior_array[self.run,0,self.subset_size + 1]   = self.extra_chi2
                self.likelihood = self.Likelihood()
                self.posterior_array[self.run,0,self.subset_size + 2]   = self.likelihood

                self.iteration = 1

                while self.iteration < mvars.max_iterations:
                    self.Walk()
                    self.iteration += 1
                    
            self.max_likelihood = np.max(self.posterior_array[:,:,self.subset_size + 2])
            self.AIC()
            self.WeightsFromPosterior()
            if mvars.number_of_MCs > 1:
                self.Convergence()
            # Tracking the full posteriors of all the subsets is very
            # memory-intensive, so we have to delete them as we go.
            del self.posterior_array

    def Walk(self):
        efvars  = self.efvars
        mvars   = self.mvars
        log     = self.log

        proceed = False

        self.prev_weights       = np.copy(self.weights)
        self.prev_saxs_chi2     = np.copy(self.saxs_chi2)
        self.prev_likelihood    = np.copy(self.likelihood)
        self.prev_ensemble_I    = np.copy(self.ensemble_I)
        if efvars.include_second_dimension:
            self.prev_extra_chi2    = np.copy(self.extra_chi2)
            self.prev_ensemble_extra= np.copy(self.ensemble_extra)
        
        while proceed==False:
            deltas = np.array([np.random.normal(loc=-0.03,scale=0.005),np.random.normal(loc=0.03,scale=0.005)],dtype=float)
            forward_or_back = np.random.randint(2,size=1)
            delta = deltas[forward_or_back]
            increment_member = np.random.randint(self.subset_size,size=1)
            if ((self.weights[increment_member] + delta) >= 0.0):
                self.weights[increment_member] += delta
                self.weights = self.weights/np.sum(self.weights)
                proceed = True
        
        self.ensemble_I = np.dot(self.weights,self.subset_basis)
        if efvars.include_second_dimension:
            self.ensemble_extra = np.dot(self.weights,self.subset_extra)
        self.Chi2()
        self.likelihood = self.Likelihood()

        acceptance_ratio        = self.likelihood/self.prev_likelihood
        acceptance_criterion    = np.power(acceptance_ratio,1/mvars.acceptance_tuner)
        draw                = np.random.uniform(low=0.0,high=1.0)
        if not ((acceptance_ratio >= 1.0) or (draw < acceptance_criterion)):
            self.weights    = np.copy(self.prev_weights)
            self.saxs_chi2  = np.copy(self.prev_saxs_chi2)
            self.likelihood = np.copy(self.prev_likelihood)
            self.ensemble_I = np.copy(self.prev_ensemble_I)
            if efvars.include_second_dimension:
                self.extra_chi2     = np.copy(self.prev_extra_chi2)
                self.ensemble_extra = np.copy(self.prev_ensemble_extra)
        for member in range(self.subset_size):
            self.posterior_array[self.run,self.iteration,member]          = self.weights[member]
        self.posterior_array[self.run,self.iteration,self.subset_size]    = self.saxs_chi2
        self.posterior_array[self.run,self.iteration,self.subset_size + 1]= self.extra_chi2
        self.posterior_array[self.run,self.iteration,self.subset_size + 2]= self.likelihood


    def Chi2(self):
        mvars   = self.mvars
        efvars  = self.efvars
        if efvars.do_shannon_sampling:
            self.ensemble_I_shan    = self.ensemble_I[efvars.shannon_samples]
            sum1 = np.sum(np.divide(np.multiply(efvars.samples_I,self.ensemble_I_shan),np.power(efvars.samples_ERR,2)),axis=0)
            sum2 = np.sum(np.power(np.divide(self.ensemble_I_shan,efvars.samples_ERR),2),axis=0)
            c_array = np.reshape(np.divide(sum1,sum2),(1,3001))
            self.ensemble_I_shan    = np.multiply(c_array,self.ensemble_I_shan)
            chi_array = np.sum(np.power(np.divide(np.subtract(efvars.samples_I,self.ensemble_I_shan),efvars.samples_ERR),2),axis=0)
            chi = np.median(chi_array)
            chi_idx = np.where(chi_array==chi)[0][0]
            c = c_array[0,chi_idx]
            self.ensemble_I = c*self.ensemble_I
            self.saxs_chi2  = chi
        else:
            sum1    = np.sum(np.divide(np.multiply(efvars.samples_I,self.ensemble_I),np.power(efvars.samples_ERR,2)))
            sum2    = np.sum(np.power(np.divide(self.ensemble_I,efvars.samples_ERR),2))
            c = np.divide(sum1,sum2)
            self.ensemble_I = c*self.ensemble_I
            chi = np.sum(np.power(np.divide(np.subtract(efvars.samples_I,self.ensemble_I),efvars.samples_ERR),2))
            self.saxs_chi2  = chi

        if efvars.include_second_dimension:
            chi2 = np.sum(np.power(np.divide(np.subtract(efvars.extra_data,self.ensemble_extra),efvars.extra_error),2))
            self.extra_chi2 = chi2
            self.total_chi2 = self.saxs_chi2 + self.extra_chi2
        else:
            self.total_chi2 = self.saxs_chi2


    def Likelihood(self):
        
        likelihood = np.exp(-self.total_chi2/2.0)

        return likelihood

    def AIC(self):

        self.aic = 2*(self.subset_size-1) + self.total_chi2

        return

    def BIC(self):
        efvars  = self.efvars
        self.bic = self.total_chi2 + ((self.subset_size - 1) * np.log(efvars.num_points))

        return

    def WeightsFromPosterior(self):
        mvars   = self.mvars
        efvars  = self.efvars
        log     = open(self.logfile,'a')

        log.write('Posterior Array:\n'+str(self.posterior_array)+'\n')
        self.mean_weights = np.average(self.posterior_array[:,mvars.posterior_burn:,:self.subset_size],axis=1)
        self.mean_weights = self.mean_weights/np.sum(self.mean_weights,axis=1,keepdims=True)
        log.write('mean_weights: '+str(self.mean_weights)+'\n')
        self.total_mean_weights = np.average(self.mean_weights,axis=0)
        log.write('total_mean_weights: '+str(self.total_mean_weights)+'\n')
        log.close()

        self.final_ensemble_I = np.dot(self.total_mean_weights,self.subset_basis)
        self.final_ensemble_extra = np.dot(self.total_mean_weights,self.subset_extra)
        return

    def Convergence(self):
        '''Gelman-Rubin Convergence Diagnostic'''
        mvars   = self.mvars
        efvars  = self.efvars
        log     = open(self.logfile,'a')
        log.write('Calculating Convergence Criterion.')

        N = float(mvars.max_iterations - mvars.posterior_burn)

        self.within_variance = np.average(np.var(self.posterior_array[:,mvars.posterior_burn:,:self.subset_size],axis=1),axis=0)
        
        self.between_variance = (N/(mvars.number_of_MCs - 1))*np.sum(np.power(np.subtract(self.mean_weights,self.total_mean_weights),2),axis=0)

        self.pooled_var = ((N-1)/N)*self.within_variance + ((mvars.number_of_MCs + 1)/(N*mvars.number_of_MCs))*self.between_variance

        self.PSRF   = np.divide(self.pooled_var,self.within_variance)
        self.Rc     = np.sqrt(((self.subset_size + 2.0)/self.subset_size)*self.PSRF)

        log.write('PSRF: '+str(self.PSRF)+'\n')
        log.write('Rc: '+str(self.Rc)+'\n')
        log.close()

        return

