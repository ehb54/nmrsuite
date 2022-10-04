import numpy as np
import math as m

#import PCS_fit_3par_simplex as pcsf3
#import PCS_fit_8par_simplex as pcsf8
#import generate_vertices as genfile
#import select_best_solution as sbs
#import generate_xyz_mtrx_pcs as genxyz
#import vector2tensor as v2t
#import rotmat2euler as r2e
#import euler2all as e2a
#import tensor2eigen_sort as t2es
#from stopwatch import Stopwatch
#import plotagreement as pa

import sassie.analyze.pcs.PCS_fit_3par_simplex as pcsf3
import sassie.analyze.pcs.PCS_fit_8par_simplex as pcsf8
import sassie.analyze.pcs.generate_vertices as genfile
import sassie.analyze.pcs.select_best_solution as sbs
import sassie.analyze.pcs.generate_xyz_mtrx_pcs as genxyz
import sassie.analyze.pcs.vector2tensor as v2t
import sassie.analyze.pcs.rotmat2euler as r2e
import sassie.analyze.pcs.euler2all as e2a
import sassie.analyze.pcs.tensor2eigen_sort as t2es
from sassie.analyze.pcs.stopwatch import Stopwatch

def pcs2tensor(self, coord = None, PCS_data = None, guesses = None, cond_num_cutoff = None, tolerance = None):
    #------------------------------------------------------------------------------------------------------
    # This function calculates a tensor using inputs files and generates statndard outputs, called in analyze_PCSs
    # Inputs:
    #    1) Coord is H-coordinates
    #    2) PCS_data is shift values from PCS input file
    #    3) guesses is a nx3 array representing x,y,z coordinates used as starting point for simplex
    #       - n is number of guesses
    #    4) cond_num_cutoff is user specified value
    #    5) k_denovo_8par is a binary input that determines whether denovo 8 parameter fit should be done
    #       - this whole functionality is not required, and can be removed with Dr. Fushman's permission
    # Calls:
    #    1) PCS_fit_3par_simplex to perform 3-parameter simplex minimization
    #------------------------------------------------------------------------------------------------------
   
    log = self.log
    log.debug('in pcs2tensor')
    pgui = self.run_utils.print_gui
    
    #Stopwatch is used to track progress of optimization calculations
    stopwatch = Stopwatch()
    stopwatch.start()

    #Set default inputs
    if cond_num_cutoff is None:
        cond_num_cutoff = 6

    ng = np.shape(guesses)[0]
    solutions_3par = float('nan') * np.ones((ng, 5)) #Placeholder array
    solutions_3par = np.zeros( (ng,5) )
    #For each guess as input, perform simplex minimization
    for i in range(ng):
        [position, chi2, cond_number] = pcsf3.PCS_fit_3par_simplex(self,PCS_data, coord, guesses[i, :], tolerance)
        solutions_3par[i, :] = np.array([position[0], position[1], position[2], chi2, cond_number])
    #pgui('---> 3-parameter fit trials results:')
    #pgui(' ')    
    
    #Select the best solution from all of the guesses
    bestsolution_3par = sbs.select_best_solution(solutions_3par,cond_num_cutoff)
    if (bestsolution_3par.size == 0):
        pgui('no acceptable solutions found!!!\n')
    else:                            #a solution exists!
       
        pgui("="*60+" \n")
        pgui('3-Parameter Fit Trials Best Solution Report:\n')

        #pgui('Best Solution: ' + str(bestsolution_3par[0:-2]) + '\n') 
        Chi2min = bestsolution_3par[-2]
        #determine the tensor
        #shift coordinates
        coord_s = np.copy(coord)
        for i in range(3):
            coord_s[:, i] -= bestsolution_3par.flat[i] # shift in x, y, z
        
        # now we have the coordinates, let's recover the tensor
        xyz_mtrx = genxyz.generate_xyz_mtrx_pcs(coord_s)
        A_inv = np.linalg.pinv(np.copy(xyz_mtrx))
        cond_num = np.linalg.cond(A_inv)
        pgui('  Chi square minimum = ' + str(round(Chi2min,4)) + '\n')
        pgui("  Condition number = " +  str(round(cond_num,4)) + '\n')
        X_vect = np.dot(A_inv, PCS_data) #solve for the tensor in vector form

        X_tensor = v2t.vector2tensor(X_vect)
        ksort = 1 #sort eigenvalues by absolute value
        [X_eigenval, X_eigenvect, Rotmat] = t2es.tensor2eigen_sort(X_tensor, ksort)
        
        #result_3par = [bestsolution_3par[0:3], X_vect.conj().T, Chi2min,cond_num]
        result_3par = float('nan') * np.ones(8)
        for i in range(8):
            if i<3:
                result_3par[i] = bestsolution_3par.flat[i]
            else:
                j = i-3
                result_3par[i] = np.transpose(X_vect[j])
        pgui('  Position of spin label = ' + ', '.join(map(str,np.round(bestsolution_3par[0:3],3).tolist())) + '\n')
        pgui('  Vector form of susceptibility tensor = [ ' + ', '.join(map(str,np.round(np.transpose(X_vect),4).tolist())) + ' ]\n')
        
        DX = "DX"
        pgui('  Susceptibility tensor eigenvalues = ' + ', '.join(map(str, np.round(X_eigenval,4).tolist())) + '\n')
        Xa = X_eigenval[2] * 3.0/2.0
        Xr = X_eigenval[0] - X_eigenval[1] 
        Rhomb = Xr / Xa
        pgui('  Axial component (Xa) = ' + str(round(Xa,4)) + '\n')
        pgui('  Rhombic component (Xr) = ' + str(round(Xr,4)) + ' ( x10^-32 m^3)' + '\n') 
        pgui('  Rhombicity (Xr/Xa)= ' + str(round(Rhomb,4)) + '\n')   
        
        #prepare the output
        backcalc_PCS = np.matmul(np.copy(xyz_mtrx), np.copy(X_vect))
        
        
        k_euler_all = 0       #show all Euler angles
        r2e.printeuler(Rotmat)
        #plot the agreement between experimental and back-calc PCSs  
        #[R_coeff,Q_factor]=pa.plot_agreement(PCS_data,backcalc_PCS,x_label,y_label)

        stopwatch.get()
        
            
        #use the 3-par solution as input to find a 8-par solution    
        guess_8par = np.empty(8)
        for a in range(8):
            if a<3:
                guess_8par[a] = np.copy(bestsolution_3par[a])
            else:
                guess_8par[a] = np.copy(X_vect[a-3])
        guess_8par = np.asarray(guess_8par)
        solution_8par,Chi2_8par,cond_number_8par = pcsf8.PCS_fit_8par_simplex(self,PCS_data,coord,guess_8par, tolerance)
        
        result_8par = np.copy(np.asarray(solution_8par))
        result_8par = np.append(solution_8par, Chi2_8par)
        result_8par = np.append(solution_8par, cond_number_8par)

        position_vect = solution_8par[0:3]
        x_vect = solution_8par[3::]

        pgui("\n"+"="*60+" \n")
        pgui('8-Parameter Fit Report:' + '\n')
        pgui('  Direct hit: Chi square = ' + str(round(Chi2_8par,4)) + '\n')
        pgui('  Condition Number = ' + str(round(cond_number_8par,4)) + '\n')
        pgui('  Position of spin label = ' + ', '.join(map(str,np.round(position_vect,3).tolist())) + '\n')
        pgui('  Vector form of susceptibility tensor = [ ' + ', '.join(map(str,np.round(np.transpose(x_vect),4).tolist())) + ' ]\n')

        X_tensor=v2t.vector2tensor(solution_8par[3::])
        [X_eigenval, X_eigenvect, Rotmat]=t2es.tensor2eigen_sort(X_tensor,ksort)
        
        Xa = X_eigenval[2] * 3./2.
        Xr = X_eigenval[0] - X_eigenval[1] 
        Rhomb = Xr/Xa

        pgui('  Susceptibility tensor eigenvalues = ' + ', '.join(map(str,np.round(X_eigenval,4).tolist())) + '\n')
        pgui('  Axial component (Xa) = ' + str(round(Xa,4)) + '\n')
        pgui('  Rhombic component (Xr) = ' + str(round(Xr,4)) + ' ( x10^-32 m^3)' + '\n')
        pgui('  Rhombicity (Xr/Xa)= ' + str(round(Rhomb,4)) + '\n\n')
 
        #check automatically if a different solution
        dist_tol=0.01
        tens_tol=0.01

        diff_pos=result_8par[0:3]-result_3par[0:3]
        diff_tens=result_8par[3:8]-result_3par[3:8]
        if (m.sqrt(np.sum(diff_pos ** 2)) < dist_tol) or (m.sqrt(np.sum(diff_tens ** 2)) < tens_tol):
            pgui('\n  The 8-param solution is not different from 3-param solution\n')
            pgui("="*60+" \n")
        else:        #it's a different solution, recalculate the Euler angles
            #print('need to re-calculate the Euler angles here') 
            output = r2e.printeuler(Rotmat)
            print (output)

    stopwatch.get()
    return backcalc_PCS, position_vect, x_vect, str(Chi2_8par), str(cond_number_8par), r2e.printeuler(Rotmat)
