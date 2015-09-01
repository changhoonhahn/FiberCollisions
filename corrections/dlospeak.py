
"""
def build_peakcorrected_fibcol_old(sanitycheck=False, **cat_corr): 
    ''' Build peak corrected fibercollided mock catalogs (using cosmolopy) 

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary
    sanitycheck : testing flag 

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # fit functions (using lambda) 
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5*x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 

    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo -----------------

        n_mocks = 160   # total number of mocks

        survey_zmin, survey_zmax = 0.16, 0.44
        
        # only one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or \
                (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
        
        # read in fiber collided mock  
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 
        cosmo = fibcoll_mock.cosmo           # survey comoslogy 

        survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']

        # set peak fraction         
        if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
            f_peak = 1.0
        else: 
            f_peak = correction['fpeak'] 

        appended_ra, appended_dec, appended_z, appended_weight = [], [], [], []
        upweight_again = []
            
        if sanitycheck == True:     # check that the peak p(r) is generated properly
            pr_test = [] 
    
        for i_mock in range(len(fibcoll_mock.weight)): 
            # go through every galaxy in fibercollided mock catalog

            while fibcoll_mock.weight[i_mock] > 1: 
                # for galaxies with wcp > 1
                fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], 
                        **cosmo)*cosmo['h']
                
                rand_num = np.random.random(1) 
                
                if rand_num <= f_peak:          # if in the peak 
                    # keep ra and dec
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])

                    if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
                        # appended galaxy has weight of 1.0 
                        appended_weight.append(correction['fpeak'])
                    else: 
                        # appended galaxy has weight of 1.0 
                        appended_weight.append(1.0)
                    if correction['fit'].lower() in ('gauss', 'expon'):   
                        # compute the displacement within peak using best-fit --------------
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 
                        #--------------------------------------------------------------------- 

                    elif correction['fit'].lower() == 'true': 
                        # compute the displacement within peak using actual distribution   
                        dlos_comb_peak_file = ((fibcoll_mock.file_name).rsplit('/', 1))[0]+'/DLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
                        dlos_mid, dlos_dist = np.loadtxt(dlos_comb_peak_file, unpack=True, usecols=[0,1])

                        dlos_cdf = dlos_dist.cumsum() 

                        rand1 = np.random.random(1) 
                        
                        cdf_closest_index = min(range(len(dlos_cdf)), key = lambda i: abs(dlos_cdf[i]-rand1[0])) 
                        closest_dlos = dlos_mid[cdf_closest_index] 
                       
                        try: 
                            closest_dloses
                        except NameError:
                            closest_dloses = [closest_dlos]
                        else: 
                            closest_dloses.append(closest_dlos)

                        rand2 = np.array([closest_dlos])
                        #--------------------------------------------------------------------- 

                    # in case the displacement falls out of bound (may general large scale issues)
                    if (comdis_imock+rand2 > survey_comdis_max) or \
                            (comdis_imock+rand2 < survey_comdis_min): 
                        collided_z = comdis2z(comdis_imock-rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock+rand2, **cosmo)

                    appended_z.append(collided_z[0]) 

                    if correction['name'] == 'allpeakshot': 
                        # need to re-upweight the fibercollided galaxies by 1-fpeak 
                        upweight_again.append(i_mock)       # save index to be re-upweighted 

                    if sanitycheck == True: 
                        pr_test.append(rand2) 
                else:                           # if not in the peak 
                    if correction['name'] == 'peaktest': 
                        # do nothing 
                        pass
                    elif correction['name'] in ('peak', 'peaknbar'): 
                        # randomly displace 
                        appended_ra.append(fibcoll_mock.ra[i_mock])
                        appended_dec.append(fibcoll_mock.dec[i_mock])
                        appended_weight.append(1.0) 
                        
                        rand1 = np.random.random(1) 
                        appended_z.append(0.16+rand1[0]*0.28)

                        if sanitycheck == True: 
                            original_dm = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                            corrected_dm = cosmos.distance.comoving_distance(0.16+rand1[0]*0.28, **cosmo)*cosmo['h']
                            pr_test.append(corrected_dm-original_dm)

                    elif (correction['name'] == 'allpeak') or (correction['name'] == 'allpeakshot'): 
                        # it should all be in the peak!
                        raise NameError('should not happen')

                    elif correction['name'] == 'peakshot': 
                        upweight_again.append(i_mock)
        
        if correction['name'] == 'peakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_upweightagain in upweight_again: 
                fibcoll_mock.weight[i_upweightagain] = fibcoll_mock.weight[i_upweightagain]+1.

        elif correction['name'] == 'allpeakshot': 
            # re-upweighting for all peak+shotnoise correction 
            for i_upweightagain in upweight_again: 
                fibcoll_mock.weight[i_upweightagain] = fibcoll_mock.weight[i_upweightagain]+(1.0-correction['fpeak']) 

        print len(appended_ra), ' galaxies were peak corrected'
        fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
        fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
        fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])

    elif catalog['name'].lower() in ('tilingmock', 'qpm', 'patchy', 'nseries'): 
        # Tiling mock, QPM, PATCHY 
        
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits

        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        elif catalog['name'].lower() == 'nseries': 
            n_mocks = 84
        
        # read in mock with fibercollisions imposed
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr)  
        
        cosmo = fibcoll_mock.cosmo  # survey cosmology 

        survey_comdis_min = cosmos.distance.comoving_distance(survey_zmin, **cosmo)*cosmo['h']         # in units of Mpc/h
        survey_comdis_max = cosmos.distance.comoving_distance(survey_zmax, **cosmo)*cosmo['h']         # in units of Mpc/h
    
        if catalog['name'].lower() != 'tilingmock':
            # only use fiber collision weights
            fibcoll_mock.weight = fibcoll_mock.wfc            
                        
        if correction['name'] == 'peaknbar': 
            # read in the true galaxies for tail portion  
            # the tail portion of the peak corrections will be generated 
            # similar to the mksample procedures for assigning the random catalog redshift 
            true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
            true_data = galaxy_data('data', **true_cat_corr)
            
            if catalog['name'].lower() in ('qpm',  'patchy', 'nseries'): 
                true_weight = true_data.wfc
            else: 
                true_weight = true_data.weight
            true_z = true_data.z

            true_weight_cum = true_weight.cumsum()
            true_weight_max = true_weight_cum.max()
            i_true = np.arange(0,len(true_weight_cum))

        # set peak fraction (f_peak in correction direction reflects the correct fpeak, 
        # but for all peak correction methods this is simply used as weight 
        if correction['name'].lower() in ('allpeakshot', 'allpeak'): 
            f_peak = 1.0
        else: 
            f_peak = correction['fpeak'] 
    
        append_ra, append_dec, append_z, append_dlos, append_weight = [], [], [], [], []

        if catalog['name'].lower() in ('qpm', 'nseries'): 
            append_comp = []   # save comp

        reupw = []

        if sanitycheck == True:     # check peak p(r) 
            pr_test = [] 
    
        for i_mock in range(len(fibcoll_mock.weight)):  # go through each galaxy 

            while fibcoll_mock.weight[i_mock] > 1:  # for galaxies with wcp > 1

                fibcoll_mock.weight[i_mock] -= 1.0

                # LOS comoving distance of the galaxy 
                comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                
                rand_num = np.random.random(1)  # determine whether galaxy is in the peak
               
                if rand_num < f_peak:          # in the peak 

                    append_ra.append(fibcoll_mock.ra[i_mock])    # keep ra and dec
                    append_dec.append(fibcoll_mock.dec[i_mock])

                    if catalog['name'].lower() in ('qpm', 'nseries'): 
                        append_comp.append(fibcoll_mock.comp[i_mock])

                    append_weight.append(1.0)       # appended galaxy has weight=1.0 
                    
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak ------------------                        
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 
                        #-------------------------------------------------------------------- 

                    elif correction['fit'].lower() == 'true': 
                        # compute displacements within peak using true distribution ------
                        dlos_comb_peak_file = ((fibcoll_mock.file_name).rsplit('/', 1))[0]+'/DLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
                        dlos_mid, dlos_dist = np.loadtxt(dlos_comb_peak_file, unpack=True, usecols=[0,1])

                        dlos_cdf = dlos_dist.cumsum() 

                        rand1 = np.sum(dlos_dist)*np.random.random(1) 
                        
                        cdf_closest_index = min(range(len(dlos_cdf)), key = lambda i: abs(dlos_cdf[i]-rand1[0])) 
                        if rand1 > dlos_cdf[cdf_closest_index]: 
                            cdf_smaller_index = cdf_closest_index
                            cdf_larger_index = cdf_closest_index+1
                        else: 
                            cdf_smaller_index = cdf_closest_index-1
                            cdf_larger_index = cdf_closest_index

                        # linear interpolations
                        closest_dlos = dlos_mid[cdf_smaller_index] + \
                                (dlos_mid[cdf_larger_index] - dlos_mid[cdf_smaller_index])*(rand1[0]-dlos_cdf[cdf_smaller_index])/(dlos_cdf[cdf_larger_index]-dlos_cdf[cdf_smaller_index])
                       
                        try: 
                            closest_dloses
                        except NameError:
                            closest_dloses = [closest_dlos]
                        else: 
                            closest_dloses.append(closest_dlos)

                        rand2 = np.array([closest_dlos])
                        #--------------------------------------------------------------------- 
                    else: 
                        raise NotImplementedError('does not work') 

                    # in case the displacement falls out of bound 
                    # NOTE: THIS IS NOT CORRECT, BUT IMPLEMENTED FOR SIMPLICITY 
                    if (comdis_imock+rand2 > survey_comdis_max) or \
                            (comdis_imock+rand2 < survey_comdis_min): 
                        collided_z = comdis2z(comdis_imock - rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock + rand2, **cosmo)

                    append_z.append(collided_z[0]) 
                    append_dlos.append(rand2)
                    
                    if correction['name'] == 'allpeakshot': 
                        # for allpeak + shot noise correction 
                        # need to re-upweight the fibercollided galaxies by 1-fpeak 
                        reupw.append(i_mock)       # save index to be re-upweighted 

                    if sanitycheck == True:             
                        # append calculated LOS displacement 
                        pr_test.append(rand2) 

                else:                                           # Not in peak 
                    if correction['name'] == 'peaktest': 
                        # do nothing. 
                        # This effectively discards the galaxies that fall into the tail 
                        pass

                    elif (correction['name'] == 'peak') or (correction['name'] == 'peaknbar'): 
                        '''
                        generate random z based on redshift distribution of galaxies 
                        '''
                        # RA, Dec remain the same 
                        # weight = 1
                        append_ra.append(fibcoll_mock.ra[i_mock])
                        append_dec.append(fibcoll_mock.dec[i_mock])
                        append_weight.append(1.0) 
                        
                        wtarg = np.random.random(1)*true_weight_max
                        zindx = np.floor(np.interp(wtarg, true_weight_cum, i_true)).astype(int)+1
                        
                        # fail safe
                        qqq = np.where(wtarg < true_weight_cum[0])[0]
                        zindx[qqq] = 0 

                        # assign redshift 
                        append_z.append(true_z[zindx]) 

                    elif correction['name'] == 'peakshot': 
                        # peak + shot noise correction 
                        # record galaxy index that need to be upweighted again 
                        reupw.append(i_mock)
                    
                    elif correction['name'] in ('allpeakshot', 'allpeak'): 
                        # for all peak correction methods, all galaxies should be in the peak!
                        raise NameError('should not happen')

        if correction['name'] == 'peakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_reupw in reupw: 
                fibcoll_mock.weight[i_reupw] += 1.0

        elif correction['name'] == 'allpeakshot': 
            # re-upweighting for all peak+shotnoise correction 
            for i_reupw in reupw: 
                fibcoll_mock.weight[i_reupw] += 1.0-correction['fpeak']

        print len(append_ra), ' galaxies were peak corrected'
    
        # append the "append" galaxies to the end of the arrays
        fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, append_ra])
        fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, append_dec])
        fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, append_weight])

        if catalog['name'].lower() in ('qpm', 'nseries'): 
            fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, append_comp])

        fibcoll_mock.z = np.concatenate([fibcoll_mock.z, append_z])
        
        if sanitycheck == True:         # output sanitycheck 
            sanitycheck_file = ''.join([
                'peak_corrected_dlos_values_', correction['name'].lower(), '_', 
                correction['fit'].lower(), '.dat']) 
            np.savetxt(sanitycheck_file,
                    np.c_[append_dlos], fmt=['%10.5f']) 
    else: 
        raise NameError('Error here') 

    peakcorr_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif catalog['name'].lower() == 'patchy': 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, 
                    fibcoll_mock.nbar, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif catalog['name'].lower() == 'tilingmock': 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    if sanitycheck == True: 
        np.savetxt(peakcorr_file+'.sanitycheck', np.c_[pr_test], fmt=['%10.5f'], delimiter='\t') 
"""
def build_peakcorrected_fibcol(doublecheck=False, **cat_corr): 
    ''' Build peak corrected fibercollided mock catalogs (using cosmolopy) 

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary
    doublecheck : save dlos values to test it's working  

    Notes
    -----
    * Currently supported peak correction methods: peakshot 
    * nbar(z) interpolation implemented for CMASS like samples  

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # fit functions (using lambda) ------------------------------------------------------------
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 
    
    # survey specific parameters 
    if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):     # LasDamasGeo 
        n_mocks = 160   # total number of mocks
        survey_zmin, survey_zmax = 0.16, 0.44
    elif catalog['name'].lower() in ('tilingmock', 'qpm', 'patchy', 'nseries'): 
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        elif catalog['name'].lower() == 'nseries': 
            n_mocks = 1 
    elif catalog['name'].lower() in ('bigmd'):             # CMASS, BigMD
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        n_mocks = 1 
    elif 'cmass' in catalog['name'].lower():             # CMASS like catalogs
        if catalog['name'].lower() == 'cmass': 
            survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        elif 'cmasslowz' in catalog['name'].lower(): 
            if '_high' in catalog['name'].lower(): 
                survey_zmin, survey_zmax = 0.5, 0.75    
            elif '_low' in catalog['name'].lower(): 
                survey_zmin, survey_zmax = 0.2, 0.5
        else: 
            raise NotImplementedError('CMASS or CMASSLOWZ combined sample')
        n_mocks = 1 
    else: 
        raise NotImplementedError('Mock Catalog not included')

    # read in fiber collided mock  
    if correction['name'] == 'bigfc_peakshot': 
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'bigfc'}}
    else: 
        fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
    fibcoll_mock = galaxy_data('data', **fibcoll_cat_corr) 
    cosmo = fibcoll_mock.cosmo      # survey comoslogy 
    
    # comoving distnace of z_min and z_max 
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
    # fiber collision weights
    if catalog['name'].lower() not in ('tilingmock', 'lasdamasgeo', 'ldgdownnz'):
        fibcoll_mock.weight = fibcoll_mock.wfc            
                        
    if correction['name'] == 'peaknbar': 
        # read in the true galaxies for tail portion  
        # the tail portion of the peak corrections will be generated 
        # similar to the mksample procedures for assigning the random catalog redshift 
        true_cat_corr = {'catalog':catalog, 'correction': {'name': 'true'}}
        true_data = galaxy_data('data', **true_cat_corr)
        
        if catalog['name'].lower() in ('qpm',  'patchy', 'nseries', 'cmass'): 
            true_weight = true_data.wfc
        else: 
            true_weight = true_data.weight
        true_z = true_data.z

        true_weight_cum = true_weight.cumsum()
        true_weight_max = true_weight_cum.max()
        i_true = np.arange(0,len(true_weight_cum))

    # import nbar(z) file for CMASS-like catalogs 
    if 'cmass' in catalog['name'].lower(): 
        # hardcoded nbar(z) files 
        if catalog['name'].lower() == 'cmass': 
            nbar_file = '/mount/riachuelo1/hahn/data/CMASS/nbar-cmass-dr12v4-N-Reid-om0p31_Pfkp10000.dat'
        elif 'cmasslowz' in catalog['name'].lower(): 
            if 'e2' in catalog['name'].lower(): 
                nbar_file = '/mount/riachuelo1/hahn/data/CMASS/dr12v5/nbar_DR12v5_CMASSLOWZE2_North_om0p31_Pfkp10000.dat'
            elif 'e3' in catalog['name'].lower(): 
                nbar_file = '/mount/riachuelo1/hahn/data/CMASS/dr12v5/nbar_DR12v5_CMASSLOWZE3_North_om0p31_Pfkp10000.dat'
            else: 
                nbar_file = '/mount/riachuelo1/hahn/data/CMASS/dr12v5/nbar_DR12v5_CMASSLOWZ_North_om0p31_Pfkp10000.dat'

        # read in nbar(z) file 
        nbar_z, nbar_nbar = np.loadtxt(nbar_file, skiprows=2, unpack=True, usecols=[0, 3]) 
        # nbar(z) interpolation function
        nbarofz = sp.interpolate.interp1d(nbar_z, nbar_nbar, kind='cubic')       

    # peak fraction         
    try: 
        f_peak = correction['fpeak']
    except KeyError: # all peak and all peak shot corrections
        f_peak = 1.0 

    appended_ra, appended_dec, appended_z, appended_weight = [], [], [], []
    upweight_again = []

    if catalog['name'].lower() in ('qpm', 'nseries'): 
        appended_comp = []   # save comp
    elif 'cmass' in catalog['name'].lower(): 
        appended_comp, appended_wsys, appended_wnoz, appended_nbar = [], [], [], [] 
            
    if doublecheck:     # check that the peak p(r) is generated properly
        dlos_values = [] 
    
    for i_mock in range(len(fibcoll_mock.weight)):  # go through every galaxy in fibercollided mock catalog
        # there is a smarter way to do this; however it's not implemented

        while fibcoll_mock.weight[i_mock] > 1:      # for galaxies with wcp > 1

            fibcoll_mock.weight[i_mock] -= 1.0

            # LOS comoving distance of the galaxy 
            comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo)*cosmo['h']
                
            rand_num = np.random.random(1) 
            if rand_num <= f_peak:          # if in the peak 
                appended_ra.append(fibcoll_mock.ra[i_mock])     # keep ra and dec
                appended_dec.append(fibcoll_mock.dec[i_mock])
                
                # keep extra columns
                if catalog['name'].lower() in ('qpm', 'nseries'):  
                    appended_comp.append(fibcoll_mock.comp[i_mock]) 
                elif 'cmass' in catalog['name'].lower(): 
                    appended_comp.append(fibcoll_mock.comp[i_mock]) 
                    appended_wsys.append(fibcoll_mock.wsys[i_mock]) 
                    appended_wnoz.append(1.0) 

                if correction['name'].lower() in ('allpeak', 'allpeakshot'): 
                    # appended galaxy has weight of peak fraction 
                    appended_weight.append(correction['fpeak'])
                else: 
                    # appended galaxy has weight of 1.0 
                    appended_weight.append(1.0)

                # compute the displacement within peak using best-fit
                if correction['fit'].lower() in ('gauss', 'expon'):   
                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = (-3.0+rand2*6.0)*correction['sigma']
                    peakpofr = fit_func(rand2, correction['sigma']) 
                    
                    while peakpofr <= rand1: 
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 

                elif correction['fit'].lower() == 'true': 
                    # compute the displacement within peak using actual distribution   
                    dlos_comb_peak_file = ''.join([
                        ((fibcoll_mock.file_name).rsplit('/', 1))[0], '/', 
                        'DLOS_norm_peak_dist_', catalog['name'].lower(), '_', str(n_mocks), 'mocks_combined.dat'])
                    dlos_mid, dlos_dist = np.loadtxt(dlos_comb_peak_file, unpack=True, usecols=[0,1])

                    dlos_cdf = dlos_dist.cumsum()/dlos_dist.sum()

                    rand1 = np.random.random(1) 
                    
                    cdf_closest_index = min(range(len(dlos_cdf)), key = lambda i: abs(dlos_cdf[i]-rand1[0])) 
                    closest_dlos = dlos_mid[cdf_closest_index] 
                   
                    try: 
                        closest_dloses
                    except NameError:
                        closest_dloses = [closest_dlos]
                    else: 
                        closest_dloses.append(closest_dlos)

                    rand2 = np.array([closest_dlos])
                else: 
                    raise NotImplementedError('asdfasdf')

                # in case the displacement falls out of bound (may general large scale issues)
                if (comdis_imock + rand2 > survey_comdis_max) or (comdis_imock + rand2 < survey_comdis_min): 
                    collided_z = comdis2z(comdis_imock-rand2, **cosmo)
                else: 
                    collided_z = comdis2z(comdis_imock+rand2, **cosmo)

                appended_z.append(collided_z[0]) 
                if 'cmass' in catalog['name'].lower():
                    appended_nbar.append(nbarofz(collided_z[0]))

                if doublecheck: 
                    dlos_values.append(rand2) 

            else:                           # if not in the peak ----------------------------
                if correction['name'] in ('peak', 'peaknbar'): 
                    '''
                    generate random z based on redshift distribution of galaxies 
                    '''
                    # RA, Dec remain the same 
                    # weight = 1
                    appended_ra.append(fibcoll_mock.ra[i_mock])
                    appended_dec.append(fibcoll_mock.dec[i_mock])
                    appended_weight.append(1.0) 
                    
                    wtarg = np.random.random(1)*true_weight_max
                    zindx = np.floor(np.interp(wtarg, true_weight_cum, i_true)).astype(int)+1
                    
                    # fail safe
                    qqq = np.where(wtarg < true_weight_cum[0])[0]
                    zindx[qqq] = 0 

                    # assign redshift 
                    appended_z.append(true_z[zindx]) 

                elif correction['name'] in ('peakshot', 'bigfc_peakshot'): 
                    upweight_again.append(i_mock)

                else: 
                    raise NotImplementedError('asdfasdf') 
    
    if correction['name'] in ('peakshot', 'bigfc_peakshot'): 
        # re-upweighting for peak+shotnoise correction 
        for i_upweightagain in upweight_again: 
            fibcoll_mock.weight[i_upweightagain] += 1.

    print len(appended_ra), ' galaxies were peak corrected'
    # append artificial galaxies to catalog 
    fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
    fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
    fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
    fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])
        
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, appended_comp])
    elif 'cmass' in catalog['name'].lower():
        fibcoll_mock.comp = np.concatenate([fibcoll_mock.comp, appended_comp])
        fibcoll_mock.wsys = np.concatenate([fibcoll_mock.wsys, appended_wsys])
        fibcoll_mock.wnoz = np.concatenate([fibcoll_mock.wnoz, appended_wnoz])
        fibcoll_mock.nbar = np.concatenate([fibcoll_mock.nbar, appended_nbar])

    peakcorr_file = get_galaxy_data_file('data', **cat_corr) 
    
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        np.savetxt(peakcorr_file, 
                np.c_[
                    fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, 
                    fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 

    elif 'cmass' in catalog['name'].lower():          # CAMSS
        head_str = "# columns : ra, dec, z, nbar, w_systot, w_noz, w_fc, comp"
        np.savetxt(peakcorr_file, 
                np.c_[
                    fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.nbar,
                    fibcoll_mock.wsys, fibcoll_mock.wnoz, fibcoll_mock.weight, fibcoll_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t', header=head_str) 

    elif catalog['name'].lower() in ('tilingmock', 'lasdamasgeo', 'ldgdownnz'): 
        np.savetxt(peakcorr_file, 
                np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    else: 
        raise NotImplementedError('asdfasdf')

    if doublecheck: 
        np.savetxt(peakcorr_file+'.dlosvalues', np.c_[dlos_values], fmt=['%10.5f'], delimiter='\t') 
