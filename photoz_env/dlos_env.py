'''

Code to investigate environment dependence on the 
line-of-sight displacement

Author(s): ChangHoon Hahn


'''
import numpy as np
from scipy.optimize import curve_fit
import sys
import os.path
import subprocess
import cosmolopy as cosmos
import multiprocessing as mp

# --- Local ---
import fibcol_data as fc_data
import fibcol_dlos as fc_dlos
import fibcol_utility as fc_util
import galaxy_environment as genv
import pyspherematch as pysph 
import mpfit as mpfit 

class DlosEnv: 
    def __init__(self, cat_corr, n_NN=3, **kwargs):
        ''' dLOS with corresponding nth nearest neighbor distances (a tracer of galaxy environment)
        '''
        self.cat_corr = cat_corr
        self.n_NN = n_NN

        self.file_name = self.File(n_NN=3, **kwargs) 
        self.dlos = None
        self.env = None 

    def File(self, **kwargs): 
        ''' Get file name for the combined dLOS + Env values 
        '''
        cat = self.cat_corr['catalog'] 
        if 'n_mock' not in cat.keys():
            cat['n_mock'] = 1 
        corr = self.cat_corr['correction']
        if corr['name'].lower() != 'upweight':
            corr = {'name': 'upweight'} 
                
        f_dlos = fc_dlos.dlos(**{'catalog': cat, 'correction': corr})         
        dlos_dir = '/'.join((f_dlos.file_name).split('/')[:-1])+'/'
    
        try: 
            nmock_str = '_' + str(kwargs['n_mocks'])+'mocks'
        except KeyError: 
            if cat['name'].lower() == 'nseries': 
                nmock_str = '_84mocks'
            else: 
                raise NotImplementedError("Only Nseries Implemented For Now") 
        
        file_name = ''.join([dlos_dir, 
            'DLOSENV_d', str(self.n_NN), 'NN_', cat['name'], nmock_str, '.dat']) 

        return file_name 

    def Read(self, **kwargs): 
        ''' Read file 
        '''
        if not os.path.isfile(self.file_name):
            # if file does not exist then make it 
            self.calculate(**kwargs)
            self.Write(**kwargs)
        try: 
            if kwargs['clobber']: 
                self.calculate(**kwargs)
                self.Write(**kwargs)
        except KeyError: 
            pass 

        self.dlos, self.env = np.loadtxt(self.file_name, skiprows=1, unpack=True, usecols=[0,1])

        return None 

    def calculate(self, **kwargs): 
        ''' Construct combined dLOS distribution with corresponding nth nearest neighbor distances (galaxy environment)

        Parameters
        ----------
        cat_corr : catalog and correction dictionary 
        n : n in nth nearest neighbor distance 
        
        Notes
        -----
        * Currently only implemented for nseries

        '''
        catalog = self.cat_corr['catalog']
        correction = {'name': 'upweight'}   # correction has to be upweight
   
        # Read and combine dLOS data 
        if catalog['name'].lower() in ('qpm', 'nseries'):    # QPM or Nseries
            
            try: 
                n_mocks = kwargs['n_mocks'] 
            except KeyError: 
                if catalog['name'].lower() == 'nseries': 
                    n_mocks = 84
                else:
                    raise NotImplementedError("Not implemented error")

            for i_mock in range(1, n_mocks+1):         # combine dLOS values of n mocks 
                # read dLOS for each mock  
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

                # import DLOS values 
                los_disp_i = fc_dlos.dlos(**i_cat_corr)         

                # compute nth nearest neighbor distance for upweighted galaxy  
                NN_dist = genv.dlos_d_NN(n=self.n_NN, **i_cat_corr ) 
                
                # combine dLOS and dNN values  
                try: 
                    combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos])
                    combined_dNN = np.concatenate([combined_dNN, NN_dist]) 
                except NameError: 
                    combined_dlos = los_disp_i.dlos
                    combined_dNN = NN_dist
        else: 
            raise NotImplementedError('asdfasdfasdf') 
        
        self.dlos = combined_dlos
        self.env = combined_dNN

        return None  

    def Write(self, **kwargs): 
        ''' Write dLOS + Env values to file
        '''
        head_str = '# dlos, dNN'
        print '#####################################'
        print 'Writing ', self.file_name 
        print '#####################################'
        np.savetxt(self.file_name, 
                np.c_[self.dlos, self.env], 
                delimiter='\t', header=head_str)

        return None 

    def fpeak_env(self, **kwargs): 
        ''' Calculate fpeak for bins of galaxy environment (dNN)
        '''
        if self.dlos == None: 
            self.Read(**kwargs)

        combined_dlos = self.dlos
        combined_dNN = self.env

        print 'd'+str(self.n_NN)+'NN, Minimum ', min(combined_dNN), ' Maximum ', max(combined_dNN)
    
        # bin dNN values
        if 'stepsize' in kwargs.keys():
            dnn_step = kwargs['stepsize']
        else: 
            dnn_step = 2

        istart = np.int(np.floor(min(combined_dNN)))
        iend = np.int(np.floor(max(combined_dNN)/np.float(dnn_step)))
        dNN_bins = [ (dnn_step*i, dnn_step*(i+1)) for i in np.arange(istart, iend) ] 
        
        for i_bin, dNN_bin in enumerate(dNN_bins): 
            # dNN bin 
            try: 
                bin_index = np.where(( combined_dNN >= dNN_bin[0] ) & ( combined_dNN < dNN_bin[1] )) 
            except TypeError: 
                bin_index = np.where((combined_dNN >= dNN_bin)) 

            bin_dlos = combined_dlos[bin_index] 
            if len(bin_dlos) < 50: 
                # if the bin contains too few dLOSs then skip
                continue
            
            # calculate the histogram and fit peak of the distribution 
            dlos_hist, mpc_mid, peak_param  = \
                    fc_dlos.dlos_hist_peak_fit(bin_dlos, fit='gauss', peak_range=[-15.0, 15.0])
            
            try:        # save avg_dNN and fpeak values 
                dNN_avg.append( 0.5 * np.float(dNN_bin[0] + dNN_bin[1]) ) 
                fpeaks.append( peak_param['fpeak'] ) 
                sigma.append( peak_param['sigma'] ) 
                ndlos.append( len(bin_dlos) )
            except NameError: 
                dNN_avg = [0.5 * np.float(dNN_bin[0] + dNN_bin[1])]
                fpeaks = [peak_param['fpeak']]
                sigma = [peak_param['sigma']]
                ndlos = [len(bin_dlos)]
            except TypeError: 
                dNN_avg.append( np.float(dNN_bin) ) 
                fpeaks.append( peak_param['fpeak'] ) 
                sigma.append( peak_param['sigma'] ) 
                ndlos.append( len(bin_dlos) )

        # MPfit ----
        p0 = [-0.01, 0.8]            # initial guesses for p0[0] * x + p0[1]
        fa = {'x': np.array(dNN_avg), 'y': np.array(fpeaks)}
        fit_param = mpfit.mpfit(fc_dlos.mpfit_linear, p0, functkw=fa, nprint=0)
        
        fit_slope = fit_param.params[0]
        fit_yint = fit_param.params[1]
        
        fit_head_str = ''.join([' Slope = ', str(fit_slope), ", y-int = ", str(fit_yint), '\n', 
            ' Column : dNN, fpeak, sigma, bestfit fpeak, N_dlos']) 

        bestfit_fpeak = np.array([
            dNN_avg[i] * fit_slope + fit_yint for i in range(len(dNN_avg))])
        print len(dNN_avg), len(bestfit_fpeak)

        dlos_dir =  '/'.join( (self.file_name).split('/')[:-1] )+'/'
        dlos_file = (self.file_name).split('/')[-1]
        fpeak_file = ''.join([dlos_dir, 'fpeak_env_', dlos_file]) 

        np.savetxt(fpeak_file, 
                np.c_[dNN_avg, fpeaks, sigma, bestfit_fpeak, ndlos],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i'], 
                delimiter='\t', header=fit_head_str)

        return [dNN_avg, fpeaks, sigma]

def fpeak_dNN(dNN, cat_corr, n_NN=3, **kwargs): 
    ''' fpeak value given nth nearest neighbor distance

    Notes
    -----
    * Upweight correction automatically assigned
    '''
    cat_corr['correction'] = {'name': 'upweight'}

    comb_dlos = DlosEnv(cat_corr, n_NN=n_NN)
    # fpeak file 
    dlos_dir =  '/'.join( (comb_dlos.file_name).split('/')[:-1] )+'/'
    dlos_file = (comb_dlos.file_name).split('/')[-1]
    fpeak_file = ''.join([dlos_dir, 'fpeak_env_', dlos_file]) 
    
    if 'truefpeak' in kwargs.keys():
        if kwargs['truefpeak']: 
            dNN_avg, fpeaks, sigmas = np.loadtxt(fpeak_file, skiprows=2, unpack=True, usecols=[0,1,2])
    else: 
        dNN_avg, fpeaks, sigmas = np.loadtxt(fpeak_file, skiprows=2, unpack=True, usecols=[0,3,2])
    
    if isinstance(dNN, float):
        dNN_list = [dNN]
    else: 
        dNN_list = dNN
    
    min_indx = [] 
    for dNN_i in dNN_list:
        min_indx.append((np.abs(dNN_avg - dNN_i)).argmin())

    return fpeaks[min_indx]

if __name__=="__main__": 
    cat_corr = {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}} 
    for nnn in [1,2,3,4,5,10]: 
        comb_dlos = DlosEnv(cat_corr, n_NN=nnn)
        comb_dlos.fpeak_env(stepsize=2)
    print fpeak_dNN([1.2, 3.4, 5.7], cat_corr, n_NN=3) 
