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

class DlosEnv: 
    def __init__(self, cat_corr, n_NN=3, **kwargs):
        ''' dLOS with corresponding nth nearest neighbor distances (a tracer of galaxy environment)
        '''
        self.cat_corr = cat_corr
        self.n_NN = n_NN

        self.file_name = self.File(n_NN=3, **kwargs) 

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

def dlos_env_multiprocess(params): 
    ''' Wrapper for dlos_env function 
    '''
    n = params[0]
    cat_corr = params[1]

    dlos_env(n=n, **cat_corr)
    return 

def build_fpeak_env(n=3, **cat_corr): 
    ''' Build best fit for fpeak(d_NN) and save parameters

    n_mock ranges hardcoded so far 
    '''
    catalog = cat_corr['catalog']
    correction = {'name': 'upweight'} 
   
    # Read in dLOS data ----------------------------------------------------------------------
    if catalog['name'].lower() == 'qpm':    # QPM --------------------------------------------

        for i_mock in range(1, 11):         # loop through mocks 
            
            # read dLOS for each file 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
            print i_cat_corr

            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # compute nth nearest neighbor distance for upweighted galaxy  
            NN_dist = genv.dlos_d_NN(n=n, **i_cat_corr ) 
            
            # combine dLOS and dNN from files 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
                combined_dNN = NN_dist
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos])
                combined_dNN = np.concatenate([combined_dNN, NN_dist]) 

    else: 
        raise NotImplementedError('asdfasdfasdf') 

    # loop through different dNN measurements 
    dNN_bins = [
            (0, 5), (5, 10), (10, 25), (25, 50)
            ]       # hardcoded right now 
    
    for i_bin, dNN_bin in enumerate(dNN_bins): 

        bin_index = ( combined_dNN >= dNN_bin[0] ) & ( combined_dNN < dNN_bin[1] ) 

        bin_dlos = combined_dlos[bin_index] 
        bin_perc = np.float( len(bin_dlos) )/np.float( len(combined_dlos) ) * 100.0

        dlos_hist, mpc_mid, peak_param  = \
                dlos_hist_peak_fit(bin_dlos, fit='gauss', peak_range=[-15.0, 15.0])
        
        try:        # save avg_dNN and fpeak values 
            dNN_avg.append( 0.5 * np.float(dNN_bin[0] + dNN_bin[1]) ) 
            fpeaks.append( peak_param['fpeak'] ) 
        except NameError: 
            dNN_avg = [0.5 * np.float(dNN_bin[0] + dNN_bin[1])]
            fpeaks = [peak_param['fpeak']]

    # MPfit dNN_avg vs fpeak ---------
    p0 = [-0.01, 0.8]            # initial guesses for p0[0] * x + p0[1]
    fa = {'x': np.array(dNN_avg), 'y': np.array(fpeaks)}
    fit_param = mpfit.mpfit(mpfit_linear, p0, functkw=fa, nprint=0)
    
    fit_slope = fit_param.params[0]
    fit_yint = fit_param.params[1]
    
    # save fpeak(dNN) for catalogs 
    fpeak_dnn_fit_file = ''.join([
        '/mount/riachuelo1/hahn/data/FiberCollisions/', 
        'fit_param_fpeak_d', str(n), 'NN_', catalog['name'].lower(), '.dat']) 
    f = open(fpeak_dnn_fit_file, 'w')
    f.write(str(fit_slope)+'\t'+str(fit_yint))
    f.close() 

def fpeak_dNN(dNN, n=3, **cat_corr): 
    ''' Calculate fpeak given dNN and catalog information 
    '''
    catalog = cat_corr['catalog']
    
    # read fpeak(dNN) best fit file
    fpeak_dnn_fit_file = ''.join([
        '/mount/riachuelo1/hahn/data/FiberCollisions/', 
        'fit_param_fpeak_d', str(n), 'NN_', catalog['name'].lower(), '.dat']) 

    a, b = np.loadtxt(fpeak_dnn_fit_file) 

    return a * dNN + b

if __name__=="__main__": 
    cat_corr = {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}} 
    pool = mp.Pool(processes=5)
    mapfn = pool.map
    
    arglist = [ [i, cat_corr] for i in [1,2,3,4,5,10]]
    
    mapfn( dlos_env_multiprocess, [arg for arg in arglist])

    pool.close()
    pool.terminate()
    pool.join() 
