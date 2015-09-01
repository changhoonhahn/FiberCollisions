
def build_corrected_randoms(sanitycheck=False, **cat_corr): 
    ''' 
    Construct corrected randoms based on corrected data
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # Las Damas Geo ------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        # read series of corrected data in order to minimize single catalog dependent variance?
        # especially for las damas, which is divided into 4 parts 
        print 'Combining Mocks'  
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                i_cat = {'name':catalog['name'], 'n_mock': i_mock, 'letter': letter} 
                i_cat_corr = {'catalog': i_cat, 'correction': correction} 
                corr_data = galaxy_data('data', **i_cat_corr)
                
                print corr_data.file_name
            
                # target weights based on the corrected weights
                if (i_mock == 1) and (letter == 'a'): 
                    targ_weight = corr_data.weight 
                    targ_z = corr_data.z
                else: 
                    targ_weight = np.concatenate((targ_weight, corr_data.weight))
                    targ_z = np.concatenate((targ_z, corr_data.z))
    # Tiling Mocks ----------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'tilingmock': 
        # build corrected randoms for tiling mocsk 
        # read corrected data (Only one mock available) 
        corr_data = galaxy_data('data', **cat_corr)
        print 'Generating random catalog based on : ', corr_data.file_name
            
        # target weights based on the corrected weights
        targ_weight = corr_data.weight 
        targ_z = corr_data.z

    # QPM Mocks --------------------------------------------------------------------------------------------------------
    elif catalog['name'].lower() == 'qpm': 
        # read series of corrected data in order to minimize single catalog dependent variance?
        # especially for las damas, which is divided into 4 parts 
        print 'Combining Mocks'  
        for i_mock in range(1,41):                                     #### CURRENTLY ONLY COMBINING 100 BECAUSE 1000 IS TOO FREAKING MANY!
            i_cat = {'name':catalog['name'], 'n_mock': i_mock} 
            i_cat_corr = {'catalog': i_cat, 'correction': correction} 
            corr_data = galaxy_data('data', **i_cat_corr)
                
            print 'combining ... ', corr_data.file_name
            
            # target weights based on the corrected weights
            if i_mock == 1: 
                targ_weight = corr_data.wfc                         # target weight IGNORES veto weights which the randoms take into account separately 
                targ_z = corr_data.z
            else: 
                targ_weight = np.concatenate((targ_weight, corr_data.wfc))
                targ_z = np.concatenate((targ_z, corr_data.z))
    else: 
        # not yet coded 
        raise NameError('not coded yet') 

    targ_weight_cum = targ_weight.cumsum()            # cumulative weights
    targ_weight_max = targ_weight_cum.max()           # max weight (a.k.a. total weight) 
    i_targ = np.arange(0,len(targ_weight_cum))
    print 'Corrected Galaxy Weight Sum = ', targ_weight_max
    
    if catalog['name'].lower() == 'tilingmock': 
        # for tiling mock import original randoms 
        rand_true = {'catalog': catalog, 'correction': {'name': 'true'}}
        true_random = galaxy_data('random', readdata=False, **rand_true)
        
        true_random_data = np.loadtxt('/mount/riachuelo1/hahn/data/tiling_mocks/randoms-boss5003-icoll012-vetoed.dat', unpack=True, usecols=[0,1]) 
        true_random.ra = true_random_data[0]
        true_random.dec = true_random_data[1]
    else:
        # read true random
        rand_true = {'catalog': catalog, 'correction': {'name': 'true'}}
        true_random = galaxy_data('random', **rand_true)

    Nran = len(true_random.ra)  # number of random galaxies
    print 'N_random = ', Nran

    corr_random_z = np.zeros(Nran) 
    # divide into chunks like make_catalog_z.py
    nchunk = 50 
    for i_chunk in range(nchunk): 
        # chunk indicies
        istart = Nran/(nchunk)*i_chunk
        iend = Nran/(nchunk)*(i_chunk+1)
        if i_chunk == nchunk-1: 
            iend = Nran 
        
        # randomly sample the sum of the target weights 
        wtarg = np.random.random(iend-istart)*targ_weight_max

        # get zindex of the target 
        zindx = np.floor(np.interp(wtarg, targ_weight_cum, i_targ)).astype(int)+1       # identical to make_catalog_z.py
        qqq = np.where(wtarg < targ_weight_cum[0])[0]
        zindx[qqq]  = 0 

        # assign redshift 
        corr_random_z[istart:iend] = targ_z[zindx]

    if len(corr_random_z) != Nran: 
        raise NameError("NOT ENOUGH CORRECTED RANDOMS!") 

    ## write Corrected Random File  
    corr_random = galaxy_data('random', readdata=False, **cat_corr)    
    corr_random_filename = corr_random.file_name 

    # write corrected randoms
    if catalog['name'].lower() == 'qpm': 
        # QPM has different data format so needs to be separate
        np.savetxt(corr_random_filename, np.c_[true_random.ra, true_random.dec, corr_random_z, true_random.nbar, true_random.wveto], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%1.5e', '%10.5f'], delimiter='\t') 
    else: 
        np.savetxt(corr_random_filename, np.c_[true_random.ra, true_random.dec, corr_random_z], 
                fmt=['%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_fc_info_catalog(**cat_corr): 
    '''
    build mock catalogs with fiber-collision information (upweighted neighbor info) using specific idl routines
    '''
    catalog = cat_corr['catalog']
    
    # QPM --------------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'qpm': 
        upw_cat_corr = {'catalog': catalog, 'correction': {'name': 'upweight'}}     # upweight cat_corr dict
        fc_upw = galaxy_data('data', readdata=False, **upw_cat_corr)                # get upweight name 
    
        fc_upw_info_file_name = '.fcinfo.'.join((fc_upw.file_name).rsplit('.', 1)) 
        
        fc_info_cat_cmd = ''.join(['idl -e "', 
            "build_fibcoll_catalog, 'qpm', ", 
            "'"+fc_upw.file_name+"'", ", '"+fc_upw_info_file_name+"'"+'"'])
        
        print fc_info_cat_cmd
        os.system(fc_info_cat_cmd) 
    else: 
        raise NameError('not yet coded') 

    return fc_info_cat_cmd

def build_tailupweight_fibcol(**cat_corr): 
    ''' Only upweight galaxies in the tail but keep everything else true 
    tail definition is where abs(dLOS) > 20 Mpc/h
    '''
    catalog = cat_corr['catalog']

    if catalog['name'].lower() == 'qpm': 
        upw_cat_corr = {'catalog': catalog, 'correction': {'name': 'upweight'}}     # upweight cat_corr dict
        fc_upw = galaxy_data('data', readdata=False, **upw_cat_corr)                # get upweight name 
    
        fc_upw_info_file_name = '.fcinfo.'.join((fc_upw.file_name).rsplit('.', 1)) 
    
        if os.path.isfile(fc_upw_info_file_name) == False:                          # if info file doesn't exist
            build_fc_info_catalog(**cat_corr) 

        # read ra, dec, z, wfkp, wcp, wcomp, d_fc from fiber-collision info catalog
        mock_ra, mock_dec, mock_z, wfkp, wcp, wcomp, dfc, fc_index = np.loadtxt(
                fc_upw_info_file_name, unpack=True, usecols=[0,1,2,3,4,5,6,10] 
                ) 

        # run through each line
        for i_row in range(len(mock_ra)): 
            if dfc[i_row] != 0.0:       # if there is some sort of fiber-collision
                if np.abs(dfc[i_row]) > 20.0:    # if in the tail then keep upweights
                    pass
                else:                            # undo upweighting scheme 
                    if wcp[i_row] != 0.0:       
                        # only down-weighted fc-pairs should have dLOS_fc != 0.0 so catch errors
                        raise NameError('something went wrong') 

                    wcp[i_row] = 1.0 
                    wcp[fc_index[i_row]] = wcp[fc_index[i_row]]-1.0         # subtract from fc weight 
    

        tailupw = galaxy_data('data', readdata=False, **cat_corr)
        np.savetxt(tailupw.file_name, 
                np.c_[
                    mock_ra, mock_dec, mock_z, wfkp, wcp, wcomp
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    else: 
        raise NameError('not yet coded')

def build_peak_fpeak_dNN(NN=3, **cat_corr): 
    ''' Build corrected mock catalogs using peak correction with *env dependent fpeak* 

    - Calculates kth nearest neighbor distances for upweighted galaxies 
    - Use bestfit fpeak(d_NN) to calculate the fpeak
    - Peak correct using fpeak 

    Input: 
    **cat_corr dictionary 
    
    (using cosmolopy) 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # lambda function 
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5*x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 

    # LDG -------------------------------------------------------------------- 
    if catalog['name'].lower() == 'lasdamasgeo': 

        survey_zmin, survey_zmax = 0.16, 0.44      # survey redshift limit
        n_mocks = 160

        # LOS comoving distance for redshift limit
        #comdis_lo = cosmos.distance.comoving_distance(0.16, **cosmo)*cosmo['h']
        #comdis_hi = cosmos.distance.comoving_distance(0.44, **cosmo)*cosmo['h'] 
        
        # in case more than one catalog parameter is defined
        if (isinstance(catalog['n_mock'], int) == False) or \
                (isinstance(catalog['letter'], str) == False): 
            raise NameError('only one mock can be corrected at a time') 
        
        # read in fiber collided galaxies
        fc_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fc_mock = galaxy_data('data', **fibcoll_cat_corr) 
        n_gal = len(fc_mock.weight)        # Ngal 

        cosmo = fc_mock.cosmo                   # import cosmology from galaxy class 
        # survey comoving distnace limits 
        survey_comdis_min = \
                cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = \
                cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']

        upw_index = np.arange(n_gal)[ fc_mock.weight > 1 ]  # index of upweighted galaxy 

        upw_dNN = genv.d_NN(fc_mock.ra[upw_index], fc_mock.dec[upw_index], fc_mock.z[upw_index], 
                n=NN, **fc_cat_corr)       # upweight kth Nearest Neighbor distance
        
        fake_ra, fake_dec, fake_z, fake_weight, re_upw = [], [], [], [], []
        for ii_upw, i_upw in enumerate(upw_index):     # go through upweighted galaxy catalog 

            while fc_mock.weight[i_upw] > 1: # for galaxies with wcp > 1

                fc_mock.weight[i_upw] = fc_mock.weight[i_upw] - 1.0

                # LOS comoving distance of the galaxy 
                comdis_iupw = cosmos.distance.comoving_distance(\
                        fc_mock.z[i_upw], **cosmo) * cosmo['h']

                fpeak_iupw = fc_dlos.fpeak_dNN( upw_dNN[ii_upw], n=NN, **fc_cat_corr ) 
                
                rand_num = np.random.random(1)      # random # sampled for the peak 
                if rand_num <= fpeak_iupw:          # IN THE PEAK  

                    # place fake galaxy a dLOS apart from the upweighted galaxy
                    fake_ra.append(fc_mock.ra[i_upw])
                    fake_dec.append(fc_mock.dec[i_upw])
                    fake_weight.append(1.0)

                    if correction['fit'].lower() in ('gauss', 'expon'):   
                        # compute the displacement within peak using best-fit 
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 

                    # if fake galaxy out of bound (may generate large scale issues)
                    # put it in the opposite sign 
                    # **THINK ABOUT THIS OMRE *** 
                    if (comdis_iupw + rand2 > survey_comdis_min) or \
                            (comdis_iupw + rand2 < survey_comdis_max): 
                        collided_z = comdis2z(comdis_iupw - rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_imock + rand2, **cosmo)

                    fake_z.append(collided_z[0])        # append collided z

                else:                       # NOT IN THE PEAK 
                    if correction['name'] == 'peakshot': 
                        re_upw.append(i_upw)
                    else: 
                        raise NotImplementedError
        
        if correction['name'] == 'peakshot': 
            # re-upweighting for peak+shotnoise correction 
            for i_reupw in re_upw: 
                fc_mock.weight[i_reupw] += 1.

        print len(fake_ra), ' galaxies were peak corrected'
        # add the fake galaxies to the catalog 
        fc_mock.ra = np.concatenate([fc_mock.ra, fake_ra])
        fc_mock.dec = np.concatenate([fc_mock.dec, fake_dec])
        fc_mock.z = np.concatenate([fc_mock.z, fake_z])
        fc_mock.weight = np.concatenate([fc_mock.weight, fake_weight])

    # QPM / Tiling Mock ------------------------------------------------
    elif catalog['name'].lower() in ('tilingmock', 'qpm'): 

        survey_zmin, survey_zmax = 0.43, 0.7        # survey redshift limit 
        
        if catalog['name'].lower() == 'qpm':
            n_mocks = 100
        
        # read in fibercollided mocks
        fc_cat_corr = {'catalog':catalog, 'correction': {'name': 'upweight'}}
        fc_mock = galaxy_data('data', **fc_cat_corr) 

        cosmo = fc_mock.cosmo                       # set comoslogy 
        # redshift limits comoving distance  
        survey_comdis_min = \
                cosmos.distance.comoving_distance(survey_zmin, **cosmo)*cosmo['h']  # in units of Mpc/h
        survey_comdis_max = \
                cosmos.distance.comoving_distance(survey_zmax, **cosmo)*cosmo['h'] 

        if catalog['name'].lower() == 'qpm': 
            # fibercollisions weights are ultimately the weights I will be using here
            fc_mock.weight = fc_mock.wfc            
        n_gal = len(fc_mock.weight)            # Ngal 
        
        upw_index = np.arange(n_gal)[ fc_mock.weight > 1 ]  # index of upweighted galaxy 

        upw_dNN = genv.d_NN(fc_mock.ra[upw_index], fc_mock.dec[upw_index], fc_mock.z[upw_index], 
                n=NN, **fc_cat_corr)       # upweight kth Nearest Neighbor distance
        print 'dNN min= ', np.min(upw_dNN), ', max =', np.max(upw_dNN)
        
        fake_ra, fake_dec, fake_z, fake_weight, re_upw = [], [], [], [], []
        if catalog['name'].lower() == 'qpm': 
            fake_wfkp = [] 
            fake_comp = [] 

        for ii_upw, i_upw in enumerate(upw_index):     # go through upweighted galaxy catalog 

            # for galaxies with wcp > 1
            while fc_mock.weight[i_upw] > 1: 
                fc_mock.weight[i_upw] = fc_mock.weight[i_upw]-1.0

                # LOS comoving distance of the galaxy 
                comdis_iupw = \
                        cosmos.distance.comoving_distance(fc_mock.z[i_upw], **cosmo)*cosmo['h']
                
                # fpeak given dNN of galaxy 
                fpeak_iupw = fc_dlos.fpeak_dNN( upw_dNN[ii_upw], n=NN, **fc_cat_corr )      

                rand_num = np.random.random(1)          
                if rand_num < fpeak_iupw:          # if in the peak 
                
                    # keep RA/Dec w = 1
                    fake_ra.append(fc_mock.ra[i_upw])
                    fake_dec.append(fc_mock.dec[i_upw])
                    fake_weight.append(1.0)

                    if catalog['name'].lower() == 'qpm': 
                        fake_wfkp.append(fc_mock.wfkp[i_upw]) 
                        fake_comp.append(fc_mock.comp[i_upw])
                    
                    if correction['fit'].lower() in ('gauss', 'expon'): 
                        # compute the displacement within peak
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*correction['sigma']
                            peakpofr = fit_func(rand2, correction['sigma']) 

                    # in case the displacement falls out of bound 
                    if (comdis_iupw + rand2 > survey_comdis_min) or \
                            (comdis_iupw + rand2 < survey_comdis_max): 

                        collided_z = comdis2z(comdis_iupw - rand2, **cosmo)
                    else: 
                        collided_z = comdis2z(comdis_iupw + rand2, **cosmo)
                    
                    fake_z.append(collided_z[0]) 
                    
                else:       # Tail ------------------------------------------------
                    if correction['name'] in ('peakshot_dnn'): 
                        # peak + shot noise correction 
                        # record galaxy index that need to be upweighted again 
                        re_upw.append(i_upw)
                    else: 
                        raise NotImplementedError

        if correction['name'] in ('peakshot_dnn'): 
            # re-upweight for peak+shotnoise correction 
            for i_reupw in re_upw: 
                fc_mock.weight[i_reupw] = fc_mock.weight[i_reupw]+1.
        else:
            raise NotImplementedError

        print len(fake_ra), ' galaxies were peak corrected'
        fc_mock.ra = np.concatenate([fc_mock.ra, fake_ra])
        fc_mock.dec = np.concatenate([fc_mock.dec, fake_dec])
        fc_mock.weight = np.concatenate([fc_mock.weight, fake_weight])
        if catalog['name'].lower() == 'qpm': 
            fc_mock.wfkp = np.concatenate([fc_mock.wfkp, fake_wfkp])
            fc_mock.comp = np.concatenate([fc_mock.comp, fake_comp])
        fc_mock.z = np.concatenate([fc_mock.z, fake_z])

    else: 
        raise NotImplementedError('Error here') 

    peakcorr = galaxy_data('data', readdata=False, **cat_corr)
    peakcorr_file = peakcorr.file_name 
    
    if catalog['name'].lower() == 'qpm': 
        # QPM has specific formatting
        np.savetxt(peakcorr_file, 
                np.c_[fc_mock.ra, fc_mock.dec, fc_mock.z, 
                    fc_mock.wfkp, fc_mock.weight, fc_mock.comp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                delimiter='\t') 
    else: 
        np.savetxt(peakcorr_file, np.c_[fc_mock.ra, fc_mock.dec, fc_mock.z, fc_mock.weight], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

