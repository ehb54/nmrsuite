'''
    SASSIE: Copyright (C) 2011-2016 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    ALTENS is the module that calculates alignment tensor extracted from RDCs

    From Professor David Fushman, UMD

    J. Magn. Reson. 138, 334-342 (1999)
    J. Magn. Reson. 168, 336-345 (2004)
    J. Magn. Reson. 201, 25-33 (2009)
    Prog. Nuc. Magn. Res. Spec. 44, 189-214 (2004)
'''

from __future__ import division

import sys, os, random, logging, numpy, string, shutil,time

import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities

# For local testing
#sys.path.append('./')
#import pcs_utils as pcs_utils

# For complied version
import sassie.analyze.pcs.pcs_utils_4plot as pcs_utils

if sasconfig.__level__ == "DEBUG": DEBUG = True

app = 'pcs'

class module_variables():
    def __init__(self, parent = None):
        self.app = app

class pcs_input_variables():

    def __init__(self, parent = None):
        pass

class pcs():

    def __init__(self, parent = None):
        pass

    def main(self, input_variables, txtOutput, plotQueues):

        '''
        main method to manage RDC calculation 
        '''

        self.mvars = module_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.run(plotQueues)

        self.epilogue()

        return

    def unpack_variables(self, variables):
        '''
        method to extract variables into system wise class instance
        '''

        mvars = self.mvars

        self.log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.pdbfile = variables['pdbfile'][0]
        mvars.dcdfile = variables['dcdfile'][0]
        mvars.pcs_input_file = variables['pcs_input_file'][0]

        mvars.residue_exclusion_file_flag = variables['residue_exclusion_file_flag'][0]
        if mvars.residue_exclusion_file_flag:
            mvars.residue_exclusion_file = variables['residue_exclusion_file'][0]

        mvars.cond_num_cutoff = variables['cond_num_cutoff'][0]
        mvars.tolerance = variables['tolerance'][0]        
      
        mvars.user_guess_flag = variables['user_guess_flag'][0]
        if mvars.user_guess_flag:
            mvars.user_guess_num = variables['user_guess_num'][0]
            mvars.user_guess_x = variables['user_guess_x'][0]
            mvars.user_guess_y = variables['user_guess_y'][0]
            mvars.user_guess_z = variables['user_guess_z'][0]
        self.log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to initialize input variables 
        type file_out: list containing file names for standard output per each frame(model), 
                       the last element is for file with average values 
                       ( # of element = number of frame + 1 )
        type file_calcrdc: list containing file names for rdc data(exp, calc, delta) per each frame(model), 
                       the last element is for file with average values 
                       ( # of element = number of frame + 1 )
        type file_mc_trajectory: list containing file names for mc trajectory per each frame
                       ( # of element = number of frame )
        '''
        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui
        mvars = self.mvars

        ## code consistent with SACSCALC, so with TAMC

        ''' directory preparation '''
        output_folder = os.path.join(mvars.runname,app)
        if os.path.isdir(output_folder):
            shutil.rmtree(output_folder)
            os.mkdir(output_folder)

        mol = sasmol.SasMol(0)
#        mol.read_pdb(mvars.pdbfile)

        try:
            if(mvars.dcdfile[-3:] == 'dcd'):
                print ">> input file is a DCD file"
                mol.read_dcd(str(mvars.dcdfile))
                self.intype = 'dcd'
            elif(mvars.dcdfile[-3:] =='pdb'):
                print ">> input file is a PDB file"
                mol.read_pdb(mvars.dcdfile)
                nf = mol.number_of_frames()
                self.intype = 'pdb'
        except:
            message = 'trouble reading your input PDB or DCD file, either your file is too large or \n'
            message += 'input filename is a PDB or DCD file but it must end with ".pdb" or ".dcd" '
            message += ' :  stopping here'
            pgui(message)
            sys.exit(1)
        mvars.number_of_frames = mol.number_of_frames()

        mvars.file_calc_pcs = ['' for x in xrange(mvars.number_of_frames + 1)]
        mvars.file_summary = ['' for x in xrange(mvars.number_of_frames + 1)]
        mvars.file_pymol = ['' for x in xrange(mvars.number_of_frames + 1)]
        return

    def run(self, plotQueues):
        '''
        method to perform Altens calculation
        '''

        log = self.log
        mvars = self.mvars
        pgui = self.run_utils.print_gui

        log.debug('in PCS')

        pgui("\n"+"="*60+" \n")
        pgui("DATA FROM RUN: %s \n\n" %(time.asctime( time.gmtime( time.time() ) ) ))

        for frame in range(0, mvars.number_of_frames):
            pcs_utils.pcs_core(self,app, frame, plotQueues)
            ''' save output and report progress '''
            fraction_done = (frame+1)/float(mvars.number_of_frames)
            time.sleep(0.01)
            pgui('STATUS\t'+str(fraction_done))

        pgui('\nProcessed %d DCD frame(s)\n'%mvars.number_of_frames)

        output_dir = os.path.join(mvars.runname, app)
        pgui('\nData stored in directory: %s\n'%output_dir)
#        pgui('\nPlots (HTML) stored in directory: %s\n'%output_dir)
        return

    def epilogue(self):
        '''
        method to print out computational results and to move output data
        to appropriate places.
        '''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        fraction_done = 1
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        pgui("\n"+"="*60+" \n")
        #pgui('\nALTENS IS DONE')

        time.sleep(2)

        return

