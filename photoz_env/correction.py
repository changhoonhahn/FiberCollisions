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
        n_NN = 5    # default n_NN = 5

    if correction['name'].lower() != 'photozenvpeakshot': 
        raise NameError("Only accepts photozpeakshot correction")
    
    if 'doublecheck' in kwargs.keys(): 
        if kwargs['doublecheck']:     # check that the PDF of dLOS peak is properly generated
            dlos_values = [] 

    # fit functions (only gaussian)
    fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    
    # survey redshift limits for mock catalogs
    if catalog['name'].lower() in ('nseries'): 
        survey_zmin, survey_zmax = 0.43, 0.7
    else: 
        raise NotImplementedError('Catalog not yet included')

    start_time = time.time() 
    # read in fiber collided mocks with assigned photometric redshift  
    fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'photoz'}}
    data = fc_data.galaxy_data('data', **fibcoll_cat_corr) 
    if 'weight' not in data.columns:        # resolve nomenclature issue
        data.weight = data.wfc            
    cosmo = data.cosmo      # survey comoslogy 
    print (time.time() - start_time),' Seconds to read photoz galaxy catalog' 
    
    # comoving distance of min and max redshifts    
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']

    f_peak = correction['fpeak']        # peak fraction 
    
    fcoll = np.where(data.weight == 0)      # fiber collided galaxies 
    # These must be preserved!
    n_fcoll = len(fcoll[0])                 # Ngal fiber collided
    n_fc_peak = int(f_peak * np.float(n_fcoll))     # Ngal fc peak 
    n_fc_tail = n_fcoll - n_fc_peak                 # Ngal fc tail 
    
    # Comoving photometric redshift line-of-sight displacement
    Dc_upw = cosmos.distance.comoving_distance( data.zupw[fcoll], **cosmo ) * cosmo['h']
    Dc_zphoto = cosmos.distance.comoving_distance( data.z_photo[fcoll], **cosmo ) * cosmo['h']
    LOS_d_photo = Dc_zphoto - Dc_upw         # photometric redshit dLOS 
   
    start_time = time.time() 
    # Calculate galaxy environment of upweighted galaxy 
    upw_dNN = genv.d_NN( (data.ra)[data.upw_index[fcoll]], (data.dec)[data.upw_index[fcoll]], (data.z)[data.upw_index[fcoll]], 
            n=n_NN, **fibcoll_cat_corr )
    print (time.time() - start_time),' Seconds to calculate environment of upweighted galaxies'

    # dNN bins 
    min_dNN, max_dNN, dNN_step = np.min(upw_dNN), np.max(upw_dNN), 2.0 
    n_bins = int((max_dNN - min_dNN)/dNN_step)
    # dNN histogram
    dNN_hist, dNN_binedges = np.histogram(upw_dNN, bins=n_bins, range=[min_dNN, max_dNN]) 
    dNN_low, dNN_high = dNN_binedges[:-1], dNN_binedges[1:]
    
    # not definitely in tail, but sampled as tail or peak 
    Ntot_peak, notdeftail_tail, notdeftail_peak = 0, [], [] 

    for i_dnn in range(len(dNN_low)):   
        # for bins of upw_dNN (galaxy environment of upweighted galaxy) 
        
        dNN_bin_fcoll = np.where( (upw_dNN >= dNN_low[i_dnn]) & (upw_dNN < dNN_high[i_dnn]) ) 
        dNN_bin = (fcoll[0])[dNN_bin_fcoll]     # original index
        n_dNN = len(dNN_bin_fcoll[0])       # number og galaxies in dNN bin 

        # Calculate the corresponding fpeak to galaxy environment bin 
        fpeak_dNN_arr = dlos_env.fpeak_dNN( 0.5*(dNN_low[i_dnn] + dNN_high[i_dnn]), 
                fibcoll_cat_corr, n_NN=n_NN)
        fpeak_dNN = fpeak_dNN_arr[0]
        n_peak_dNN = int(np.float(n_dNN) * fpeak_dNN)
        n_tail_dNN = n_dNN - n_peak_dNN
    
        # phometric dLOS for dNN bin 
        dlos_photo_dNN_bin = LOS_d_photo[dNN_bin_fcoll]
            
        # definitely in the tail from photometric redshift 
        def_tail_dlos_bin = np.where( (dlos_photo_dNN_bin < -175.0) | (dlos_photo_dNN_bin > 175.0) ) 
        def_tail_bin = dNN_bin[def_tail_dlos_bin]
        n_def_tail_dNN = len(def_tail_bin)
    
        print str(dNN_low[i_dnn]), ' < dNN < ', str(dNN_high[i_dnn])
        print 'Ngal fiber collided', n_dNN, ' fpeak = ', fpeak_dNN
        print 'Ngal fiber collided tail', n_tail_dNN 
        print 'Ngal fiber collided definitely in tail', n_def_tail_dNN
        
        # New fpeak_dNN; excluding galaxies that are definitely in the tail  
        try: 
            fpeak_not_tail_bin = 1.0 - np.float(n_tail_dNN - n_def_tail_dNN)/np.float(n_dNN - n_def_tail_dNN)
        except ZeroDivisionError: 
            continue
        print 'fpeak excluding fibercollided galaxies in tail', fpeak_not_tail_bin 
    
        # Then the rest of the collided galaxies are not definitely in the tail 
        not_def_tail_bin = list( set(list(dNN_bin)) - set(list(def_tail_bin)) )
        print 'Ngal fiber collided not definitely in tail', len(not_def_tail_bin)

        for i_mock in not_def_tail_bin:         # for each galaxy that is not definitely in the tail,  
            # use new fpeak, to randomly classify peak/tail fibercollided galaxies  
            rand_num = np.random.random(1)
            
            if rand_num <= fpeak_not_tail_bin:  # randomly sampled in peak 
                notdeftail_peak.append(i_mock)
                Ntot_peak += 1
            else:                               # randomly sampled in tail 
                notdeftail_tail.append(i_mock)    

    # Preserve Ngal_fc_peak by randomly sampling tail galaxies not definitely in tail 
    n_extra_tail = n_fc_peak - Ntot_peak  
    extra_tail = random.sample(notdeftail_tail, n_extra_tail)

    notdeftail_peak += extra_tail
    if len(notdeftail_peak) != n_fc_peak: 
        raise ValueError("Ngal_fc_peak must be conserved!")
    
    data.weight[ (data.upw_index)[notdeftail_peak] ] -= 1.0
    data.weight[notdeftail_peak] += 1.0
    if data.weight[notdeftail_peak] > 1.0:
        raise ValueError("Fibercollided galaxies after correction should never be greater than 1.0")
                
    # comoving distance of upweighted galaxy 
    comdis_upw = cosmos.distance.comoving_distance(
            data.zupw[notdeftail_peak], **cosmo) * cosmo['h']
    
    for ii, i_mock in enumerate(notdeftail_peak): 
        # compute line-of-sight displacement within the peak using best-fit function 
        if correction['fit'].lower() in ('gauss'): 
            # Gaussian best-fit function

            rand1 = np.random.random(1)
            rand2 = np.random.random(1)
            rand2 = (-3.0 + 6.0 * rand2) * correction['sigma']  # +/- 3 sigmas 
            peak_pofr = fit_func(rand2, correction['sigma'])     # P(r)

            while peak_pofr <= rand1: 
                rand1 = np.random.random(1)
                rand2 = np.random.random(1)
                rand2 = (-3.0 + 6.0 * rand2) * correction['sigma']
                peak_pofr = fit_func(rand2, correction['sigma'])     # P(r)
        else: 
            raise NotImplementedError("Only Gaussian best-fit function implemented for peak correction")

        if (comdis_upw[ii] + rand2 > survey_comdis_max) or (comdis_upw[ii] + rand2 < survey_comdis_max): 
            collided_z = fc_data.comdis2z(comdis_upw[ii] - rand2, **cosmo)      # convert comoving distnace to z 
        else: 
            collided_z = fc_data.comdis2z(comdis_upw[ii] + rand2, **cosmo)
        
        if 'doublecheck' in kwargs.keys(): 
            if kwargs['doublecheck']: 
                dlos_values.append(rand2) 

        data.z[i_mock] = collided_z[0]

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

if __name__=="__main__": 
    cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'photozenvpeakshot', 'fit': 'gauss', 'sigma': 4.0, 'fpeak': 0.69}
            }
    build_photoz_env_dlospeak_fibcol(cat_corr)
