'''

Code to investigate environment dependence on the 
line-of-sight displacement

Author(s): ChangHoon Hahn


'''
import numpy as np
import scipy as sp 
import os.path
import cosmolopy as cosmos

# --- Local ---
from dlos import Dlos
from galenv.galenv import d_NN

class DlosEnv(Dlos): 

    def __init__(self, cat_corr, n_NN=3, **kwargs):
        """ Child class of Dlos class that describes the nth nearest neighbor distance of 
        upweighted galaxies in fiber collided pair

        Notes
        -----
        * Very clunky because it has to communicate with dLOS parent class 
        """
        self.n_NN = n_NN

        super(DlosEnv, self).__init__(cat_corr, **kwargs)

        self.dlos = None
        self.env = None 

        self.file_name = self.file()
        self.dlos_file = super(DlosEnv, self).file()

    def file(self): 
        """ Name of dLOS + galaxy environment file
        """

        dlos_filename = super(DlosEnv, self).file()

        nNN_str = ''.join([
            'DLOS_d', 
            str(self.n_NN),
            'NN_'
            ]) 
        
        file_name = nNN_str.join( 
                dlos_filename.split('DLOS_')
                ) 

        return file_name 

    def build(self): 
        """ Calculate the nth nearest neighbor distances for fiber collided pairs 
        """
        self.kwargs.pop('clobber', None)
        if self.dlos == None: 
            super(DlosEnv, self).read()

        # dNN for upweighted galaxy in fiber collided 
        # pairs. Each fiber collided pair corresponds to 
        # a value of dLOS 
        NN_dist = d_NN( 
                self.upw_ra, 
                self.upw_dec, 
                self.upw_z, 
                self.cat_corr, 
                n_NN = self.n_NN 
                ) 
    
        # each value of d_NN corresponds to a dLOS value 
        # in dLOS file 
        print self.file_name
        np.savetxt(self.file_name, 
                np.c_[NN_dist], 
                fmt=['%10.5f'],
                header='Columns : d_NN'
                ) 

        return None 

    def read(self, **kwargs): 
        """ Read both dLOS and dNN values 
        """

        if not os.path.isfile(self.file_name):
            self.build()
        elif 'clobber' in self.kwargs.keys():
            if self.kwargs['clobber']: 
                self.build()
        
        if not os.path.isfile(self.dlos_file):
            super(DlosEnv, self).build()
        
        # read dLOS file from parent class
        super(DlosEnv, self).read()

        self.env = np.loadtxt(
                self.file_name, 
                skiprows=1, 
                unpack=True, 
                usecols=[0]
                )

        return None 

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
        dNN_list = np.array([dNN])
    else: 
        dNN_list = dNN
    
    #min_indx = [] 
    #for dNN_i in dNN_list:
    #    min_indx.append((np.abs(dNN_avg - dNN_i)).argmin())

    #interp_fpeak = sp.interpolate.interp1d(dNN_avg, fpeaks, kind='linear')
        
    return np.interp( dNN, dNN_avg, fpeaks )

if __name__=="__main__": 
    pass
    #cat_corr = {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}} 
    #for nnn in [4,10]: 
    #    comb_dlos = DlosEnv(cat_corr, n_NN=nnn)
    #    comb_dlos.fpeak_env(stepsize=2)
    #print fpeak_dNN([1.2, 3.4, 5.7], cat_corr, n_NN=3) 


"""
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
            if len(bin_dlos) < 5: 
                # if the bin contains too few dLOSs then skip
                continue
            
            # calculate the histogram and fit peak of the distribution 
            dlos_hist, mpc_mid, peak_param  = \
                    fc_dlos.dlos_hist_peak_fit(bin_dlos, fit='gauss', peak_range=[-15.0, 15.0])
            
            if peak_param['fpeak'] == 0.0: 
                continue 
            fpeak_err = np.sqrt(peak_param['fpeak'] / np.float(len(bin_dlos)))
            
            try:        # save avg_dNN and fpeak values 
                dNN_avg.append( 0.5 * np.float(dNN_bin[0] + dNN_bin[1]) ) 
                fpeaks.append( peak_param['fpeak'] ) 
                fpeaks_err.append( fpeak_err )
                sigma.append( peak_param['sigma'] ) 
                ndlos.append( len(bin_dlos) )
            except NameError: 
                dNN_avg = [0.5 * np.float(dNN_bin[0] + dNN_bin[1])]
                fpeaks = [peak_param['fpeak']]
                fpeaks_err = [ fpeak_err ]
                sigma = [peak_param['sigma']]
                ndlos = [len(bin_dlos)]
            except TypeError: 
                dNN_avg.append( np.float(dNN_bin) ) 
                fpeaks.append( peak_param['fpeak'] ) 
                fpeaks_err.append( fpeak_err )
                sigma.append( peak_param['sigma'] ) 
                ndlos.append( len(bin_dlos) )

        # MPfit ----
        p0 = [-0.01, 0.8]            # initial guesses for p0[0] * x + p0[1]
        fa = {'x': np.array(dNN_avg), 'y': np.array(fpeaks), 'err': np.array(fpeaks_err)}
        fit_param = mpfit.mpfit(fc_dlos.mpfit_linear, p0, functkw=fa)
        
        fit_slope = fit_param.params[0]
        fit_yint = fit_param.params[1]
        
        fit_head_str = ''.join([' Slope = ', str(fit_slope), ", y-int = ", str(fit_yint), '\n', 
            ' Column : dNN, fpeak, sigma, bestfit fpeak, N_dlos']) 
        print fit_head_str

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
"""
