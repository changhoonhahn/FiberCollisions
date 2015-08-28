'''

Code for the Photometric z + Galaxy environment 
fiber collision correction method 

'''
import numpy as np
import scipy as sp 
import time 
import random
import os.path
import subprocess
import cosmolopy as cosmos
import warnings 
import matplotlib.pyplot as plt

# --- Local ----
import fibcol_data as fc_data
import fibcol_dlos as fc_dlos
import fibcol_nbar as fc_nbar
import photoz as photoz 
import galaxy_environment as genv
import pyspherematch as pysph
from utility.fitstables import mrdfits
import dlos_env

def build_photoz_env_dlospeak_fibcol(cat_corr, **kwargs): 
    ''' Build mock catalog with fiber collision correction using: 
    * Photometric redshift
    * Galaxy environment
    * dLOS distribution peak

    Notes
    -----
    * Interfaces with fibcol_data module 
    * Indexing gets very complex

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    # n th nearest neighbor 
    if 'n_NN' in kwargs.keys(): 
        n_NN = kwargs[n_NN]
    else: 
        n_NN = 5
    
    if 'doublecheck' in kwargs.keys(): 
        if kwargs['doublecheck']:     # check that the PDF of dLOS peak is properly generated
            dlos_values = [] 

    if correction['name'].lower() != 'photozenvpeakshot': 
        raise NameError("Only accepts photozpeakshot correction")

    # fit functions (only gaussian)
    fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    
    # redshift limits 
    if catalog['name'].lower() in ('nseries'): 
        # set up mock catalogs 
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        n_mocks = 1 
    else: 
        raise NotImplementedError('Catalog not yet included')

    # read in fiber collided mocks with assigned photometric redshift  
    start_time = time.time() 
    fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'photoz'}}
    data = fc_data.galaxy_data('data', **fibcoll_cat_corr) 
    cosmo = data.cosmo      # survey comoslogy 
    print (time.time() - start_time),' Seconds to read photoz galaxy catalog' 
    
    # comoving distance of min and max redshifts    
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
    if 'weight' not in data.columns:
        # resolve nomenclature issue
        data.weight = data.wfc            

    f_peak = correction['fpeak']        # peak fraction 
    
    fcoll = np.where(data.weight == 0)     # fiber collided
    # These must be preserved!
    n_fcoll = len(fcoll[0])                 # Ngal fiber collided
    n_fc_peak = int(f_peak * np.float(n_fcoll))     # Ngal fc peak 
    n_fc_tail = n_fcoll - n_fc_peak                 # Ngal fc tail 
    
    # Comoving photometric redshift line-of-sight displacement
    Dc_upw = cosmos.distance.comoving_distance(
            data.zupw[fcoll], **cosmo) * cosmo['h']
    Dc_zphoto = cosmos.distance.comoving_distance(
            data.z_photo[fcoll], **cosmo) * cosmo['h']
    LOS_d_photo = Dc_zphoto - Dc_upw         # photometric redshit dLOS 
   
    # Calculate galaxy environment of upweighted galaxy 
    start_time = time.time() 
    upw_dNN = genv.d_NN( (data.ra)[data.upw_index[fcoll]], (data.dec)[data.upw_index[fcoll]], (data.z)[data.upw_index[fcoll]], 
            n=n_NN, **fibcoll_cat_corr )
    print (time.time() - start_time),' Seconds to calculate environment of upweighted galaxies'

    # dNN bins 
    min_dNN, max_dNN, dNN_step = np.min(upw_dNN), np.max(upw_dNN), 2.0 
    n_bins = int((max_dNN - min_dNN)/dNN_step)

    dum, dNN_binedges = np.histogram(upw_dNN, bins=n_bins, range=[min_dNN, max_dNN]) 
    dNN_low, dNN_high = dNN_binedges[:-1], dNN_binedges[1:]
    
    # for bins of upw_dNN (galaxy environment of upweighted galaxy) 
    for i_dnn in range(len(dNN_low)):   

        dNN_bin_fcoll = np.where( (upw_dNN >= dNN_low[i_dnn]) & (upw_dNN < dNN_high[i_dnn]) ) 
        dNN_bin = (fcoll[0])[dNN_bin_fcoll]
        n_dNN = len(dNN_bin_fcoll[0])       # number og galaxies in dNN bin 

        # Calculate the corresponding fpeak to galaxy environment bin 
        fpeak_dNN = dlos_env.fpeak_dNN( 0.5*(dNN_low[i_dnn] + dNN_high[i_dnn]), 
                fibcoll_cat_corr, n_NN=n_NN)
        n_peak_dNN = int(np.float(n_dNN) * fpeak_dNN)
        n_tail_dNN = n_dNN - n_peak_dNN
    
        # phometric dLOS for dNN bin 
        dlos_photo_dNN_bin = LOS_d_photo[dNN_bin_fcoll]

        def_tail_dlos_bin = np.where( (dlos_photo_dNN_bin < -175.0) | (dlos_photo_dNN_bin > 175.0) ) 
        def_tail_bin = dNN_bin[def_tail_dlos_bin]
        n_def_tail_dNN = len(def_tail_bin)
    
        print str(dNN_low[i_dnn]), ' < dNN < ', str(dNN_high[i_dnn])
        print 'Ngal fiber collided', n_dNN 
        print 'Ngal fiber collided tail', n_tail_dNN 
        print 'Ngal fiber collided definitely in tail', n_def_tail_dNN
    
        # Then the rest of the collided galaxies are not definitely in the tail 
        not_def_tail_bin = list(
                set(list(dNN_bin)) - set(list(def_tail_bin))
                )
        print 'Ngal fiber collided not definitely in tail', len(not_def_tail_bin)

        # upweighted galaxy indices
        upw_def_tail_bin = (data.upw_index)[def_tail_bin]
        upw_def_not_tail_bin = (data.upw_index)[not_def_tail_bin]
    
        not_tail_fpeak_bin = 1.0 - np.float(n_tail_dNN - n_def_tail_dNN)/np.float(n_dNN - n_def_tail_dNN)
        print 'fpeak of not definitely in tail', not_tail_fpeak_bin 

    '''
    n_peakcorrected = 0     # Ngal peak corrected
    for i_mock in not_def_tail: 
        # go through each fibercollided galaxy not definitely in the tail        
        
        # use new fpeak that excludes definitely tail galaxies
        # to deterimine whether galxay is in peak or not 
        rand_num = np.random.random(1)      # random number

        if rand_num <= not_tail_fpeak:     # sampled in the peak 
            # sample a dLOS from the best-fit Gaussian dLOS PDF then place 
            # the collided galaxy with wfc = 1 dLOS away from the upweighted galaxy 
            # and then downweight the upweighted galaxy 

            data.weight[ (data.upw_index)[i_mock] ] -= 1.0  # downweight UPW galaxy 
            data.weight[i_mock] += 1.0  # upweight collided galaxy 
            if data.weight[i_mock] > 1.0:
                raise NameError('something went wrong') 

            # comoving distance of upweighted galaxy
            comdis_upw = cosmos.distance.comoving_distance(
                    data.z[ (data.upw_index)[i_mock] ], **cosmo) * cosmo['h']
        
            # compute the displacement within the peak using best-fit function 
            if correction['fit'].lower() in ('gauss'):  # Gaussian 

                rand1 = np.random.random(1) # random number
                
                # random dLOS +/- 3-sigma of the distribution 
                rand2 = np.random.random(1)
                rand2 = (-3.0 + 6.0 * rand2) * correction['sigma']
                
                peak_pofr = fit_func(rand2, correction['sigma']) # probability distribution
                
                while peak_pofr <= rand1: 
                    rand1 = np.random.random(1) # random number
                    
                    # random dLOS +/- 3-sigma of the distribution 
                    rand2 = np.random.random(1)
                    rand2 = (-3.0 + 6.0 * rand2) * correction['sigma']
                    
                    peak_pofr = fit_func(rand2, correction['sigma']) # probability distribution
                    
            else: 
                NotImplementedError('Not yet implemented') 
                
            # in case the displaced coliided galaxy falls out of bound (may generate large scale issues)
            # this will skew the dLOS displacement slightly at low and high redshift limits 
            if (comdis_upw + rand2 > survey_comdis_max) or (comdis_upw + rand2 < survey_comdis_min): 
                collided_z = comdis2z(comdis_upw-rand2, **cosmo)
            else:
                collided_z = comdis2z(comdis_upw+rand2, **cosmo)

            if doublecheck: 
                # append sample peak dLOS value to doublecheck 
                dlos_values.append(rand2) 

            #print data.z[ (data.upw_index)[i_mock] ], data.z[i_mock], collided_z[0] 
            data.z[i_mock] = collided_z[0]

            n_peakcorrected += 1
    
    print n_peakcorrected, ' galaxies were peak corrected'

    # write to file based on mock catalog  
    corrected_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() in ('nseries'): 
        np.savetxt(corrected_file, 
                np.c_[data.ra, data.dec, data.z, data.weight, data.comp, 
                        data.zupw, data.upw_index, data.z_photo], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i', '%10.5f'], 
                delimiter='\t') 
    else: 
        raise NotImplementedError('asdfasdf')
    
    if doublecheck: 
        np.savetxt(corrected_file+'.dlosvalues', np.c_[dlos_values], fmt=['%10.5f'], delimiter='\t') 
    '''

if __name__=="__main__": 
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'photozenvpeakshot', 'fit': 'gauss', 'sigma': 4.0, 'fpeak': 0.69}
            }
    build_photoz_env_dlospeak_fibcol(cat_corr)
