
def build_photoz_peakcorrected_fibcol(doublecheck=False, **cat_corr): 
    ''' Build peak corrected fibercollided mock catalogs using assigned photometric redshifts
    that emulated photometric redshift errors (using cosmolopy) 

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary
    doublecheck : save dlos values to test it's working  

    Notes
    -----
    * Addresses a number of inefficiences in the previous peakshot code. 
        * No longer appends new peak corrected galaxies; instead the original galaxy's redshift is replaced with 
        peak corrected redshift. 
        * No longer downweights all upweighted galaxies only to re-upweight them later for collided galaxies that reside in 
        the tail fiber collided galaxies
    * Choise of +/- 3 sigma for peak sampled dLOS may need to be reviewed. 

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if correction['name'].lower() != 'photozpeakshot': 
        raise NameError("Only accepts photozpeakshot correction")

    # fit functions (using lambda) ------------------------------------------------------------
    fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    
    # redshift limits 
    if catalog['name'].lower() in ('nseries'): 
        # set up mock catalogs 
        survey_zmin, survey_zmax = 0.43, 0.7    # survey redshift limits
        n_mocks = 1 
    else: 
        raise NotImplementedError('Catalog not yet included')

    # read in fiber collided mocks with assigned photometric redshift  
    fibcoll_cat_corr = {'catalog':catalog, 'correction': {'name': 'photoz'}}

    data = galaxy_data('data', **fibcoll_cat_corr) 
    cosmo = data.cosmo      # survey comoslogy 
    
    # comoving distance of min and max redshifts    
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
    if 'weight' not in data.columns:
        # resolve nomenclature issue
        data.weight = data.wfc            

    f_peak = correction['fpeak']        # peak fraction 
    
    fcoll = np.where(data.weight == 0)     # fiber collided
    n_fcoll = len(fcoll[0])                 # Ngal fiber collided
    n_fc_peak = int(f_peak * np.float(n_fcoll))     # Ngal fc peak 
    n_fc_tail = n_fcoll - n_fc_peak                 # Ngal fc tail 
    
    # Comoving distance of upweighted galaxy
    Dc_upw = cosmos.distance.comoving_distance(
            data.zupw[fcoll], **cosmo) * cosmo['h']
    # Comoving distance of fibcollided galaxy photoz
    Dc_zphoto = cosmos.distance.comoving_distance(
            data.z_photo[fcoll], **cosmo) * cosmo['h']
   
    # photometric redshit dLOS 
    LOS_d_photo = Dc_zphoto - Dc_upw
    
    # first, determine the collided galaxies that are definitely in the tail based on 
    # photometric redshift 
    def_tail_dlos = np.where( (LOS_d_photo < -175.0) | (LOS_d_photo > 175.0) ) 
    def_tail = (fcoll[0])[def_tail_dlos]
    n_def_tail = len(def_tail) 
    print 'Ngal fiber collided', n_fcoll
    print 'Ngal fiber collided tail', n_fc_tail 
    print 'Ngal fiber collided definitely in tail', n_def_tail 
    
    # Then the rest of the collided galaxies are not definitely in the tail 
    not_def_tail = list(
            set(list(fcoll[0])) - set(list(def_tail))
            )
    print 'Ngal fiber collided not definitely in tail', len(not_def_tail)
    upw_def_tail = (data.upw_index)[def_tail]
    upw_def_not_tail = (data.upw_index)[not_def_tail]
    
    not_tail_fpeak = 1.0 - np.float(n_fc_tail - n_def_tail)/np.float(n_fcoll - n_def_tail)
    print 'fpeak of not definitely in tail', not_tail_fpeak 

    if doublecheck:     # check that the PDF of dLOS peak is properly generated
        dlos_values = [] 

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
                raise NotImplementedError('Not yet implemented') 
                
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
