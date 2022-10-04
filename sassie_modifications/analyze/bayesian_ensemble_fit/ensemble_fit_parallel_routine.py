'''
    ENSEMBLE_FIT is the module that contains the functions
    that are used to fit an ensemble of structres to 
    experimental scattering data with optional additional data

    REFERENCE:

    Bowerman et al.
    Journal of Chemical Theory and Computation, 13, pp 2418-2429 (2017)

    Bowerman et al.
    Placeholder for Methods Paper...
'''
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import ceil
import os,sys
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import zipfile
import itertools
import time
import math
import pickle
import argparse
from scipy.special import binom as binom
from mpi4py import MPI

app = 'bayesian_ensemble_fit'

######## MPI Environment ########
comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()
#################################

class module_variables():
    def __init__(self, parent = None):
        self.app=app

class efvars():
    ''' _E_nsemble _F_it _VAR_iable_S_ '''
    def __init__(self, parent = None):
        pass

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

    def main(self,mpickle='',efpickle=''):
        '''
        main method to handle iterative Bayesian MC routine
        '''
        self.UnpackVariables(mpickle,efpickle)
        self.EnsembleFit()
        self.Epilogue()
        os.system('echo \"Rank '+str(rank)+' completed all tasks.\" >>'+self.logfile)
        quit()

    def UnpackVariables(self,mpickle,efpickle):
        '''
        extract variables from gui/mimic into system wide class instance
        '''
        self.mvars = pickle.load(open(mpickle,'rb'))
        self.efvars= pickle.load(open(efpickle,'rb'))
        mvars = self.mvars
        efvars= self.efvars
        self.logfile = os.path.join(efvars.output_folder,mvars.runname+'_runtime.log')

    def PickleVars(self):
        mvars  = self.mvars
        efvars = self.efvars
        if rank==0:
            os.system('echo \"In PickleVars.\" >> '+self.logfile)

            mvarspickle = os.path.join(efvars.output_folder,mvars.runname+"_mvars.p")
            pickle.dump(mvars,open(mvarspickle,'wb'))
            os.system('echo \"mvars have been written to pickle: '+mvarspickle+'\" >> '+self.logfile)
        
            efvarspickle = os.path.join(efvars.output_folder,mvars.runname+'_efvars.p')
            pickle.dump(efvars,open(efvarspickle,'wb'))
            os.system('echo \"efvars have been written to pickle: '+efvarspickle+'\" >> '+self.logfile)
            for thread in range(1,size):
                comm.send(1,dest=thread)
        else:
            comm.recv(source=0)
        return

    def EnsembleFit(self):
        '''
        This is the actual iterative routine
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        self.tic= time.time()
        self.min_aic = 99999.9
        self.Individuals()
        self.FindBest()
        efvars.best_saxs_chi2 = efvars.subspace_dict[str(efvars.best_model)].saxs_chi2
        efvars.best_extra_chi2= efvars.subspace_dict[str(efvars.best_model)].extra_chi2
        efvars.best_total_chi2= efvars.subspace_dict[str(efvars.best_model)].total_chi2

        if rank == 0:
            trackf = open(efvars.aic_file,'a')
            trackf.write(str(1)+'\t'\
                         +str(self.min_aic)+'\t'\
                         +str(efvars.best_total_chi2)+'\n')
            trackf.close()


        if mvars.use_all == True:
            efvars.best_model = str(efvars.id_array)
            efvars.subspace_dict[efvars.best_model] = simulated_basis(self,efvars.id_array)
            efvars.subspace_dict[efvars.best_model].BayesMC()
            model_object = efvars.subspace_dict[efvars.best_model]
            trackf = open(efvars.aic_file,'a')
            trackf.write(str(efvars.number_of_profiles)+'\t'\
                         +str(model_object.aic)+'\t'\
                         +str(model_object.total_chi2)+'\n')
            trackf.close()
        else:
            for subsize in range(2,efvars.number_of_profiles + 1):
                self.current_subsize = subsize
                if rank==0:
                    self.Progress()
                    sets = list(itertools.combinations(efvars.id_array,subsize))
                    for thread in range(1,size):
                        comm.send(sets,dest=thread)
                else:
                    sets=comm.recv(source=0)
                all_sets_as_list = []
                self.Nsets = len(np.atleast_1d(sets))
                
                for tracker in range(rank,self.Nsets,size):
                    as_list = np.array([],dtype=int)
                    for submember in sets[tracker]:
                        as_list = np.append(as_list,submember)
                    os.system("echo \"Rank "+str(rank)\
                               +': Working on sub-basis '+str(as_list)\
                               +'\" >> '+self.logfile)
                    all_sets_as_list.append(str(as_list))
                    efvars.subspace_dict[str(as_list)] = simulated_basis(self,as_list)
                    efvars.subspace_dict[str(as_list)].BayesMC()
                
                #The second boolean is to protect from case of sets < threads
                if (rank!=0 and rank < np.min([size,self.Nsets])):
                    os.system('echo \"Rank '+str(rank)\
                               +': Sharing sub-basis information'\
                               +' with head node\" >> '+self.logfile)
                    comm.send(efvars.subspace_dict,dest=0)
                    comm.recv(source=0)
                    os.system('echo \"Rank '+str(rank)\
                               +': Received \'continue\' signal'\
                               +' from head node\" >> '+self.logfile)
                elif (rank >= self.Nsets):
                    os.system('echo \"Rank '+str(rank)\
                               +': Waiting for \'continue\' signal'\
                               +' from head node\" >> '+self.logfile)
                    comm.recv(source=0)
                else:
                    for thread in range(1,np.min([size,self.Nsets])):
                        #dum_subspace_list = comm.recv(source=thread)
                        #Nsubsets = comm.recv(source=thread)
                        #for subset in range(int(Nsubsets)):
                        #    efvars.subspace_dict[dum_subspace_list[subset]] = comm.recv(source=thread)
                        dum_subspace_dict = comm.recv(source=thread)
                        efvars.subspace_dict.update(dum_subspace_dict)
                    for thread in range(1,size):
                        comm.send(1,dest=thread) 
                
                self.FindBest(sublist=all_sets_as_list,subsize=self.current_subsize)
                if self.subset_min_aic > self.min_aic:
                    if rank==0:
                        os.system('echo \"Best subset of size '+str(subsize)\
                                  +' is a poorer fit than '+str(subsize-1)\
                                  +' according to AIC.\"'\
                                  +' >> '+self.logfile) 
                    if not mvars.every:
                        break
                else:
                    if rank==0:
                        os.system('echo \"Best subset of size '+str(subsize)\
                                  +' is an improvement over '+str(subsize-1)\
                                  +' according to AIC.'\
                                  +' Expanding subset size.\"'\
                                  +' >> '+self.logfile)
                    	self.min_aic    = self.subset_min_aic
                    	efvars.best_model = self.subset_best_model
                    	efvars.best_saxs_chi2 = efvars.subspace_dict[str(efvars.best_model)].saxs_chi2
                    	efvars.best_extra_chi2= efvars.subspace_dict[str(efvars.best_model)].extra_chi2
                    	efvars.best_total_chi2= efvars.subspace_dict[str(efvars.best_model)].total_chi2
                        for thread in range(1,size):
                            comm.send(self.min_aic,dest=thread)
                    else:
                        self.min_aic = comm.recv(source=0)
        if rank==0:
            os.system('echo \"Best model is: '+str(efvars.best_model)+', '\
                      +'AIC = '+str(self.min_aic)+', '\
                      +'saxs_chi2 = '+str(efvars.best_saxs_chi2)+', '\
                      +'extra_chi2 = '+str(efvars.best_extra_chi2)+'\"'\
                      +' >> '+self.logfile)
            os.system('echo \"STATUS\t0.9999\" >> '+efvars.status_file)
        if mvars.calc_relative_likelihoods:
            self.RelativeLikelihoods()
        
        return

    def Individuals(self):
        '''
        This will find the fitting capabilities of each individual basis member.
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        os.system('echo \"Finding Individual Fits\" >> '+self.logfile)
        if rank==0:
            for idx in range(efvars.number_of_profiles):
                efvars.subspace_dict[str(idx)]=simulated_basis(self,[idx])
                efvars.subspace_dict[str(idx)].BayesMC()
            for thread in range(1,size):
                comm.send(efvars.subspace_dict,dest=thread)
        else:
            efvars.subspace_dict=comm.recv(source=0)

        return

    def FindBest(self,sublist=np.array([],dtype=str),subsize=None):
        '''
        Identifies the best model from the subset, and its relative likelihood to the next best fit.
        '''
        mvars   = self.mvars
        efvars  = self.efvars
        if rank == 0:
            tic = time.time()
            os.system('echo \"Finding best available model\" >> '\
                     + self.logfile)
        if len(sublist)==0:
            if rank==0:
                for key in efvars.subspace_dict:
                    key_AIC = efvars.subspace_dict[key].aic
                    if key_AIC < self.min_aic:
                        self.min_aic    = key_AIC
                        efvars.best_model = key
                for thread in range(1,size):
                    comm.send(self.min_aic,dest=thread)
                    comm.send(efvars.best_model,dest=thread)
            else:
                if subsize==None:
                    self.min_aic    = comm.recv(source=0)
                    efvars.best_model = comm.recv(source=0)
                else:
                    self.subset_min_aic     = comm.recv(source=0)
                    self.subset_best_model  = comm.recv(source=0)
                    self.subset_best_chi    = comm.recv(source=0)
        else:
            self.subset_min_aic  = 9999999.9
            for key in sublist:
                key_AIC = efvars.subspace_dict[str(key)].aic
                os.system('echo \"'+str(key_AIC)+','\
                          +str(self.subset_min_aic)+'\" >> '\
                          +self.logfile)
                if key_AIC < self.subset_min_aic:
                    self.subset_min_aic     = key_AIC
                    self.subset_best_model  = key
                    self.subset_best_chi    = efvars.subspace_dict[key].total_chi2
            if (rank!=0 and rank < np.min([size,self.Nsets])):
                comm.send(self.subset_min_aic,dest=0)
                comm.send(self.subset_best_model,dest=0)
                comm.send(self.subset_best_chi,dest=0)
            if rank==0:
                for thread in range(1,np.min([size,self.Nsets])):
                    dum_min_aic    = comm.recv(source=thread)
                    dum_best_model = comm.recv(source=thread)
                    dum_best_chi   = comm.recv(source=thread)
                    if dum_min_aic < self.subset_min_aic:
                        self.subset_min_aic    = dum_min_aic
                        self.subset_best_model = dum_best_model
                        self.subset_best_chi   = dum_best_chi
                efvars.aic_tracker = np.append(efvars.aic_tracker,\
                                               self.subset_min_aic)
                efvars.chi_tracker = np.append(efvars.chi_tracker,\
                                               self.subset_best_chi)
                trackf = open(efvars.aic_file,'a')
                trackf.write(str(self.current_subsize)+'\t'\
                             +str(self.subset_min_aic)+'\t'\
                             +str(self.subset_best_chi)+'\n')
                trackf.close()
                for thread in range(1,size):
                    comm.send(self.subset_min_aic,dest=thread)
                    comm.send(self.subset_best_model,dest=thread)
                    comm.send(self.subset_best_chi,dest=thread)
                os.system('echo \"Best sub-basis is: '\
                          + str(self.subset_best_model)\
                          +', AIC: '+str(self.subset_min_aic)+'\"'\
                          +' >> '+self.logfile)
                
            else:
                self.subset_min_aic    = comm.recv(source=0)
                self.subset_best_model = comm.recv(source=0)
                self.subset_best_chi   = comm.recv(source=0)

        return

    def RelativeLikelihoods(self):
        mvars   = self.mvars
        efvars  = self.efvars
        
        if rank==0:
            os.system('echo \"Calculating Relative Likelihoods\" >> '\
                      + self.logfile)

            outfile = open(os.path.join(efvars.output_folder,mvars.runname+'_relative_likelihoods.dat'),'w')
            outfile.write('#relative_likelihood,ID1,ID1_pop,ID1_pop_std,ID2,ID2_pop,ID2_pop_std,...,IDn,IDn_pop,ID2_pop_std,normalized chi^2,AIC\n')
            for key in efvars.subspace_dict:
                model = efvars.subspace_dict[key]
                model.relative_likelihood = math.exp((self.min_aic-model.aic)/2.0)
                outfile.write(str(np.around(model.relative_likelihood,decimals=3))+',')
                if model.subset_size > 1:
                    for idx in range(model.subset_size):
                        outfile.write(str(efvars.name_array[model.subset_members[idx]])+','+str(np.around(model.total_mean_weights[idx],decimals=3))+','+str(np.around(model.weights_std[idx],decimals=3))+',')
                    outfile.write(str(np.around(model.total_chi2/efvars.num_points,decimals=3))+','+str(np.around(model.aic,decimals=3))+'\n')
                else:
                    outfile.write(str(efvars.name_array[model.subset_members[0]])+',1.000,0.000,'+str(np.around(model.total_chi2/efvars.num_points,decimals=3))+','+str(np.around(model.aic,decimals=3))+'\n')
            outfile.close()
            for thread in range(1,size):
                comm.send(1,dest=thread)
        else:
            comm.recv(source=0)
        return


    def Epilogue(self):
        '''
        Saves files and prints info
        '''
        mvars   = self.mvars
        efvars  = self.efvars

        if rank==0:
            print_model = efvars.subspace_dict[efvars.best_model]
            outfile = open(os.path.join(efvars.output_folder,mvars.runname+'_final_model.dat'),'w')
            outfile.write('#ScatterFile avg_weight std_weight lone_saxs_chi2 lone_extra_chi2 lone_total_chi2\n')
            if print_model.subset_size > 1:
                for idx in range(print_model.subset_size):
                    ##  FIXED?  ##
                    subid = print_model.subset_members[idx]
                    outfile.write(\
                            efvars.name_array[print_model.subset_members[idx]]\
                            +' '\
                            +str(np.around(print_model.total_mean_weights[idx],decimals=3))\
                            +' '\
                            +str(np.around(print_model.weights_std[idx],decimals=3))\
                            +' '\
                            +str(np.around(efvars.subspace_dict[str(subid)].saxs_chi2,decimals=3))\
                            +' '\
                            +str(np.around(efvars.subspace_dict[str(subid)].extra_chi2,decimals=3))\
                            +' '\
                            +str(np.around(efvars.subspace_dict[str(subid)].total_chi2,decimals=3))\
                            +'\n')
            else:
                outfile.write(efvars.name_array[print_model.subset_members[0]]+'\n')
            outfile.write('SAXS_chi^2  (normalized): '\
                          +str(np.around(print_model.saxs_chi2,decimals=3))\
                          +' ('+str(np.around(print_model.saxs_chi2/efvars.num_q,decimals=3))+')'\
                          +'\n')
            if efvars.include_second_dimension:
                    outfile.write('Extra_chi^2 (normalized): '\
                                  +str(np.around(print_model.extra_chi2,decimals=3))\
                                  +' ('+str(np.around(print_model.extra_chi2/efvars.num_extra,decimals=3))+')'\
                                  +'\n')
            else:
                    outfile.write('Extra_chi^2 (normalized): 0.000 (0.000)\n')
            outfile.write('Total_chi^2 (normalized): '\
                           +str(np.around(print_model.total_chi2,decimals=3))\
                           +' ('+str(np.around(print_model.total_chi2/efvars.num_points,decimals=3))+')'\
                           +'\n')
            if not mvars.use_bic:
                outfile.write('AIC                     : '\
                          +str(np.around(print_model.aic,decimals=3))\
                          +'\n')
            else:
                outfile.write('BIC                     : '\
                          +str(np.around(print_model.aic,decimals=3))\
                          +'\n')
            outfile.close()
            

            SAS_ensemble = efvars.subspace_dict[efvars.best_model].ensemble_I
            ensemblefile = os.path.join(efvars.output_folder,\
                                        mvars.runname+'_ensemble.int')
            ofile = open(ensemblefile,'w')
            for qidx in range(len(efvars.Q)):
                ofile.write(str(efvars.Q[qidx])+'\t'\
                            +str(SAS_ensemble[qidx])+'\n')
            ofile.close()

            if efvars.include_second_dimension:
                extra_ensemble = efvars.subspace_dict[efvars.best_model].ensemble_extra
                extrafile = os.path.join(efvars.output_folder,mvars.runname+'_ensemble_extra.dat')
                ofile2 = open(extrafile,'w')
                for eidx in range(len(efvars.extra_data)):
                    ofile2.write(str(extra_ensemble[eidx])+'\n')
                ofile2.close()
                
            for thread in range(1,size):
                comm.send(1,dest=thread)
            
            if not mvars.debug:
                outputfiles = os.listdir(efvars.output_folder)
                for item in outputfiles:
                    if item.endswith('BayesMC.log'):
                        os.remove(os.path.join(efvars.output_folder,item))
            self.toc = time.time()
            runtime = self.toc - self.tic
            os.system("echo \"Bayesian Ensemble Fit routine complete in "+str(runtime/3600.)+" hours.\" >> "+self.logfile)
            self.CheckPoint()
            self.createPostPlots()
        else:
            comm.recv(source=0)

    def CheckPoint(self):
        '''
        Pickles the current efvars object.
        '''
        mvars  = self.mvars
        efvars = self.efvars
        efvarspickle = os.path.join(efvars.output_folder,\
                                    mvars.runname+'_efvars_checkpoint.p')
        if os.path.isfile(efvarspickle):
            os.remove(efvarspickle)
        pickle.dump(efvars,open(efvarspickle,'wb'))
        
        return

    def Progress(self):
        efvars = self.efvars
        mvars  = self.mvars

        try:
            Ncomplete = self.cumul_sum[self.current_subsize-1]
            Npercent  = (Ncomplete/self.total_combos)
        except:
            self.cumul_sum = np.zeros(efvars.number_of_profiles,dtype=int)
            for idx in range(1,efvars.number_of_profiles+1):
                count = binom(efvars.number_of_profiles,idx)
                if idx > 1:
                    self.cumul_sum[idx-1] = self.cumul_sum[idx-2] + count
                else:
                    self.cumul_sum[idx-1] = count
            
            self.total_combos = float(self.cumul_sum[-1]) + 4

            Ncomplete = self.cumul_sum[self.current_subsize-1]
            Npercent  = (Ncomplete/self.total_combos)
        
        #Don't let it count backwards if Nperc too low
        if Npercent < 0.0001:
            Npercent = Npercent + 0.0001
        #Convert 100% to something lower to reflect epilogue processes
        if Npercent==1.0:
            Npercent = 0.99
        
        statf = open(efvars.status_file,'w')
        statf.write('STATUS\t%f' % Npercent)
        statf.close()

        return
         

    def createPostPlots(self):
    	'''
	Creates relative probability of model and chi^2 histograms
	'''
	efvars = self.efvars
        mvars  = self.mvars
	#Load relative probability and chi^2 values into dictionary as list, using subset_size as key
	relProb_dict = {}
	chiModel_dict = {}
        chiMin_dict = {}
	chiDim1_dict = {}
	chiDim2_dict = {}
	for key in efvars.subspace_dict:
	    model = efvars.subspace_dict[key]
	    subSize = model.subset_size
            #Chi^2 model
	    if subSize in chiModel_dict:
	        chiModel_dict[subSize].append(model.total_chi2/efvars.num_points)
	    else:
	        chiModel_dict[subSize] = [model.total_chi2/efvars.num_points]
	    #Chi^2_min
	    if mvars.use_bic:
	        chiMinValue = (model.aic - (subSize-1)*np.log(efvars.num_points))/efvars.num_points
            else:
	        chiMinValue = (model.aic - 2*(subSize-1))/efvars.num_points
	    if subSize in chiMin_dict:
	        chiMin_dict[subSize].append(chiMinValue)
	    else:
	        chiMin_dict[subSize] = [chiMinValue]
	    #Individual dimensions
            if efvars.include_second_dimension:
	        if subSize in chiDim1_dict:
		    chiDim1_dict[subSize].append(model.saxs_chi2/efvars.num_q)
		    chiDim2_dict[subSize].append(model.extra_chi2/efvars.num_extra)
		else:
		    chiDim1_dict[subSize] = [model.saxs_chi2/efvars.num_q]
		    chiDim2_dict[subSize] = [model.extra_chi2/efvars.num_extra]
            #Rel. Prob.
	    if not mvars.calc_relative_likelihoods:	#Only plot relative likelihoods if the calculation is requested
	    	continue
	    elif subSize in relProb_dict:
	        relProb_dict[subSize].append(model.relative_likelihood)
	    else:
	        relProb_dict[subSize] = [model.relative_likelihood]
	
	#Convert dictionaries to list
	chiModel_list = []
	relProb_list = []
	chiMin_list = []
	label_list = []
	chiDim1_list = []
	chiDim2_list = []
	for subSize in chiModel_dict:
	    chiModel_list.append(chiModel_dict[subSize])
	    chiMin_list.append(chiMin_dict[subSize])
	    label_list.append('%s Member'%str(subSize))
	    if subSize > 1:
	        label_list[-1]+='s'		#Add the plural if needed
	    if mvars.calc_relative_likelihoods:
	        relProb_list.append(relProb_dict[subSize])
            if efvars.include_second_dimension:
	        chiDim1_list.append(chiDim1_dict[subSize])
		chiDim2_list.append(chiDim2_dict[subSize])

	#Plot each key of relProb and chiModel dictionaries
        cMapName = 'nipy_spectral'
        colorMap = plt.get_cmap(cMapName)
        numOfSizes = len(chiModel_dict)
        colors = colorMap(np.arange(0.0,1.0,1.0/numOfSizes))
        alignment='mid'
        displayWidth = .8
        numberOfRows_legend = 99
        plt.rc('legend', fontsize = 9)

	#Total (reduced) Chi^2 plot
        numOfBins = 40		#Wont be the true number of bins in this case
        bestModelChi = efvars.subspace_dict[efvars.best_model].total_chi2/efvars.num_points
        minimum = bestModelChi
        maximum = minimum + 2.0
        width = (maximum-minimum)/numOfBins
        numOfLeftBins = 4
        bins = np.arange(minimum-(numOfLeftBins+.5)*width,maximum+.5*width,width)
        plt.xlabel('Model $\chi^2$ (reduced)')
        plt.ylabel('Counts')
        plt.hist(chiModel_list,bins=bins,color=colors,stacked=True,label=label_list,rwidth=displayWidth,align=alignment)
        plt.axvline(x=minimum,label='Best model',linestyle='-.',color='black')
        plt.tight_layout()
        numOfColumns = int(len(chiModel_list)/numberOfRows_legend)+1
        ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*.8,box.height])
        plt.legend(ncol=numOfColumns,bbox_to_anchor=(1,.5),loc='center left')
        plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_chiOfModels.pdf'))         
        plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_chiOfModels.png'),dpi=400)
        plt.clf()

        #Minimum (reduced) Chi^2 plot
        plt.xlabel('Minimum $\chi^2$ in MC (reduced)')
        plt.ylabel('Counts')
	if mvars.use_bic:
            bestMinChi = (efvars.subspace_dict[efvars.best_model].aic - (efvars.subspace_dict[efvars.best_model].subset_size - 1)*np.log(efvars.num_points))/efvars.num_points
	else:
            bestMinChi = (efvars.subspace_dict[efvars.best_model].aic - 2*(efvars.subspace_dict[efvars.best_model].subset_size - 1))/efvars.num_points
        numOfBins=40
	minimum = bestMinChi
	maximum = minimum + 2.0
	width = (maximum-minimum)/numOfBins
        numOfLeftBins = 4
	bins = np.arange(minimum-(numOfLeftBins+.5)*width,maximum+.5*width,width)
        plt.hist(chiMin_list,bins=bins,color=colors,stacked=True,label=label_list,rwidth=displayWidth,align=alignment)
	plt.tight_layout()
	plt.axvline(x=bestMinChi,label='Best model',linestyle='-.',color='black')
        numOfColumns = int(len(chiMin_list)/numberOfRows_legend)+1
	ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0,box.y0,box.width*.8,box.height])
        plt.legend(ncol=numOfColumns,bbox_to_anchor=(1,.5),loc='center left')
        plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_chiMin.pdf'))
        plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_chiMin.png'),dpi=400)
	plt.clf()

        #Individual (reduced) Chi^2 plot
	if efvars.include_second_dimension:
            plt.figure(figsize=(7,10))
            plt.subplot(311)
            plt.xlabel('$\chi_{SAXS}^2$ (reduced)')
            plt.ylabel('Counts')
            numOfBins = 40
            minimum = np.amin([np.amin(chiDim1_list),np.amin(chiDim2_list)])
            maximum = np.amax([np.amin(chiDim1_list),np.amin(chiDim2_list)]) + 2.0
            width = (maximum-minimum)/numOfBins
            numOfLeftBins = 4
            bins = np.arange(minimum-(numOfLeftBins+.5)*width,maximum+.5*width,width)
            plt.hist(chiDim1_list,bins=bins,color=colors,stacked=True,label=label_list,rwidth=displayWidth,align=alignment)
            plt.legend(ncol=4,bbox_to_anchor=(0.,1.02,1,.102),loc=3,mode='expand',borderaxespad=0)
            ax = plt.gca()
            box = ax.get_position()
            ax.set_position([box.x0,box.y0*1,box.width*1.05,box.height*.9])
            plt.subplot(312)
            plt.xlabel('$\chi_{secondSet}^2$ (reduced)')
            plt.ylabel('Counts')
            plt.hist(chiDim2_list,bins=bins,color=colors,stacked=True,label=label_list,rwidth=displayWidth,align=alignment)
            ax = plt.gca()
            box = ax.get_position()
            ax.set_position([box.x0,box.y0*1,box.width*1.05,box.height*.9])
            plt.subplot(313)
            plt.xlabel('$\chi_{total}^2$ (reduced)')
            plt.ylabel('Counts')
            plt.hist(chiModel_list,bins=bins,color=colors,stacked=True,label=label_list,rwidth=displayWidth,align=alignment)
            plt.axvline(x=minimum,label='Best model',linestyle='-.',color='black')
            numOfColumns = int(len(chiMin_list)/numberOfRows_legend)+1
            ax = plt.gca()
            box = ax.get_position()
            ax.set_position([box.x0,box.y0,box.width*1.05,box.height*.9])
            
            plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_chiDims.pdf'))
            plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_chiDims.pdf'))
            plt.clf()

        #Relative probability plot
        if mvars.calc_relative_likelihoods:
            plt.figure(figsize=(6.4,4.8))
            numOfBins = 40
            minimum = 0.0
            maximum = 1.0
            binWidth = (maximum-minimum)/numOfBins
            bins = np.arange(minimum,maximum+binWidth,binWidth)
            plt.xlabel('Relative Model Probability')
            plt.ylabel('Counts')
            plt.hist(relProb_list,bins=bins,color=colors,rwidth=displayWidth,align=alignment,stacked=True,label=label_list)
            plt.tight_layout()
            numOfColumns = int(len(relProb_list)/numberOfRows_legend)+1
            ax = plt.gca()
            box = ax.get_position()
            ax.set_position([box.x0,box.y0,box.width*.8,box.height])
            plt.legend(ncol=numOfColumns,bbox_to_anchor=(1,.5),loc='center left')
            #Set yrange to a reasonable amount
            relProb_array = [y for x in relProb_list for y in x] #Puts list of lists into a single list
            for i in range(3,len(bins)):
                count = len(np.where(np.array(relProb_array)>=bins[i])[0])
                if count > 0:
                    break
            plt.ylim([0,count])
            plt.xlim([-.05,1.05])
            plt.xticks(np.arange(0.0,1.2,0.2))
            plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_relativeProbabilities.pdf'))
            plt.savefig(os.path.join(efvars.output_folder,mvars.runname+'_relativeProbabilities.png'),dpi=400)
            plt.clf()

        return 
	    


            
#######################################################################        
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
        self.as_string      = ''
        for ID in self.subset_members:
            self.as_string  = self.as_string+'_'+str(ID)
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
        self.logfile = os.path.join(efvars.output_folder,mvars.runname+'_Rank'+str(rank)+"_BayesMC.log")
        os.system('echo \"Running on sub-basis: '+str(self.subset_members)+'\" >> '+self.logfile)
        #The posterior array has form [id_1,id_2,...,id_n,SAXSchi^2,Dimension2chi^2,Likelihood]
        if self.subset_size == 1:
            self.ensemble_I = np.copy(self.subset_basis[0])
            if efvars.include_second_dimension:
                self.ensemble_extra = np.copy(self.subset_extra[0])
            self.Chi2()
            self.max_likelihood = self.Likelihood()
            if not mvars.use_bic:
                self.AIC()
            else:
                self.BIC()
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
                
            self.WeightsFromPosterior()
            if not mvars.use_bic:
                self.AIC()
            else:
                self.BIC()
            if mvars.number_of_MCs > 1:
                self.Convergence()
            self.Epilogue()

    def Walk(self):
        efvars  = self.efvars
        mvars   = self.mvars

        proceed = False

        self.prev_weights       = np.copy(self.weights)
        self.prev_saxs_chi2     = np.copy(self.saxs_chi2)
        self.prev_likelihood    = np.copy(self.likelihood)
        self.prev_ensemble_I    = np.copy(self.ensemble_I)
        if efvars.include_second_dimension:
            self.prev_extra_chi2    = np.copy(self.extra_chi2)
            self.prev_ensemble_extra= np.copy(self.ensemble_extra)

        if mvars.walk_one:
            while proceed==False:
                delta = np.random.normal(scale=mvars.max_delta)
                increment_member = np.random.randint(self.subset_size,size=1)
                if ((self.weights[increment_member] + delta) >= 0.0):
                    self.weights[increment_member] += delta
                    self.weights = self.weights/np.sum(self.weights)
                    proceed = True
        else:
            deltas = np.random.normal(scale=mvars.max_delta,\
                                      size=self.subset_size)
            self.weights = self.weights + deltas
            self.weights = self.weights/np.sum(self.weights)
        zeroers = np.where(self.weights < mvars.weight_threshold)[0]
        self.weights[zeroers] = 0.0
        if np.sum(self.weights) == 0:
            os.system('echo \"Sum of weights is 0!'\
                     +' Try reducing the \'delta\' or'\
                     +' \'zeroing_threshold\' parameters\"'\
                     +' >> ' + self.logfile)
        self.weights = self.weights/np.sum(self.weights)

        self.ensemble_I = np.dot(self.weights,self.subset_basis)
        if efvars.include_second_dimension:
            self.ensemble_extra = np.dot(self.weights,self.subset_extra)
        self.Chi2()
        self.likelihood = self.Likelihood()
        if self.prev_likelihood == 0.0:
            self.prev_likelihood = 1e-8
        acceptance_ratio        = self.likelihood/self.prev_likelihood
        draw                = np.random.uniform(low=0.0,high=1.0)
        if not ((acceptance_ratio >= 1.0) or (draw < acceptance_ratio)):
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
        
        likelihood = math.exp(-self.total_chi2/2.0)

        return likelihood

    def AIC(self):

        if self.subset_size > 1:
            self.maxlike = np.max(self.posterior_array[:,:,-1])
            self.aic = 2*(self.subset_size-1) \
                       - 2 * np.log(self.maxlike)
        else:
            self.aic = self.total_chi2
        return

    def BIC(self):
        efvars  = self.efvars
        if self.subset_size > 1:
            self.maxlike = np.max(self.posterior_array[:,:,-1])
            self.aic = ((self.subset_size - 1) * np.log(efvars.num_points))\
                       - 2 * np.log(self.maxlike)
        else:
            self.aic = self.total_chi2
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
        self.weights_std = np.zeros(np.shape(self.total_mean_weights),dtype=float)
        for idx in range(len(self.weights_std)):
            self.weights_std[idx] = np.std(self.posterior_array[:,mvars.posterior_burn:,idx])
        self.ensemble_I = np.dot(self.total_mean_weights,self.subset_basis)
        if efvars.include_second_dimension:
            self.ensemble_extra = np.dot(self.total_mean_weights,self.subset_extra)
        self.Chi2()
        #if mvars.print_maxlike:
        #    burnt  = self.posterior_array[:,mvars.posterior_burn:,:]
        #    maxlike= np.max(burnt[:,:,-1])
        #    maxid  = np.where(burnt == maxlike)
        #    minchi = burnt[maxid[0][0],maxid[1][0],-3]
        #    maxlike_weights = burnt[maxid[0][0],maxid[1][0],:self.subset_size]
        #    log.write('Weights of Maximum Observed Likelihood: (Chi^2 = '\
        #              +str(minchi)+')\n')
        #    log.write(str(maxlike_weights)+'\n')
        log.close()

    def Epilogue(self):
        mvars  = self.mvars
        efvars = self.efvars
        del self.subset_basis
        if mvars.d_max > 0.0:
            del self.ensemble_I_shan
        # Tracking the full posteriors of all the subsets is very
        # memory-intensive, so we have to delete them as we go.
        del self.posterior_array
        if efvars.include_second_dimension:
            del self.subset_extra
        del self.efvars
        del self.mvars
        return

    def Convergence(self):
        '''Gelman-Rubin Convergence Diagnostic'''
        mvars   = self.mvars
        efvars  = self.efvars
        log     = open(self.logfile,'a')
        log.write('Calculating Convergence Criterion.\n')

        N = float(mvars.max_iterations - mvars.posterior_burn)

        self.within_variance = np.average(np.var(self.posterior_array[:,mvars.posterior_burn:,:self.subset_size],axis=1),axis=0)
        
        self.between_variance = (N/(mvars.number_of_MCs - 1))*np.sum(np.power(np.subtract(self.mean_weights,self.total_mean_weights),2),axis=0)

        self.pooled_var = ((N-1)/N)*self.within_variance + ((mvars.number_of_MCs + 1)/(N*mvars.number_of_MCs))*self.between_variance

        self.PSRF   = np.divide(self.pooled_var,self.within_variance)
        self.Rc     = np.sqrt(((self.subset_size + 2.0)/self.subset_size)*self.PSRF)

        log.write('PSRF: '+str(self.PSRF)+'\n')
        log.write('Rc: '+str(self.Rc)+'\n')
        log.close()
        del self.PSRF
        del self.Rc
        return
########################################################################################

##### Command Line Parsing #####
'''
    What follows below is the command-line parsing from the main
    bayesian_ensemble_fitting module.  It builds and runs the
    above classes.
'''
parser=argparse.ArgumentParser()
parser.add_argument("-mpick",dest='mpick',type=str,help='Pickled module_variables() object.')
parser.add_argument("-efpick",dest='efpick',type=str,help='Pickled ensemble_fitting_variables() object.')
args=parser.parse_args()

reweighting = ensemble_routine()
reweighting.main(mpickle=args.mpick,efpickle=args.efpick)
