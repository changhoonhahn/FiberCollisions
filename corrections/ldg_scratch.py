
def build_ldg_scratch(**cat_corr): 
    ''' Quick function to test fiber collision correction methods on LasDamasGeo mocks
        

    Parameters
    ----------
    cat_corr : catalog and correction dictionary

    Notes
    -----
    * scratch_peak_ang : Instead of using the

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
        
    omega_m = 0.25  # cosmology for LDG

    cosmo = {} 
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.676
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 

    d_peak = 21.0
        
    # read rdzw file 
    data_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
    orig_file = ''.join([data_dir, 'sdssmock_gamma_lrgFull_zm_oriana', 
        str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.fibcoll.dat']) 

    orig_ra, orig_dec, orig_z, orig_wfc, z_upw, upw_index = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,3,4,5], 
            dtype={'names': ('ra', 'dec', 'z', 'wfc', 'zupw', 'upw_index'), 
                'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.int32)})

    if correction['name'].lower() in ('scratch_peakknown'): 
        now_index = np.where(orig_wfc < 1)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < d_peak)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        
        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 
            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

    elif correction['name'].lower() in ('scratch_peakknown_ang'):

        now_index = np.where(orig_wfc < 1)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < d_peak)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

    elif correction['name'].lower() in ('scratch_peakknown_gauss'): 

        now_index = np.where(orig_wfc < 1)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < d_peak)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -1.0*d_peak + 2.0 * d_peak * rand1
            peakpofr = fit_func(rand2, 6.5) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -1.0*d_peak + 2.0 * d_peak * rand1
                peakpofr = fit_func(rand2, 6.5) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
    elif correction['name'].lower() in ('scratch_allpeak_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 1.0
            orig_wfc[upw_index[i_now]] -= 1.0
    
    elif correction['name'].lower() in ('scratch_peakknown_gauss_divide'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 
    
        append_ra, append_dec, append_z, append_wfc, append_wcomp = [], [], [], [], [] 
        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            for jj in range(1,5): 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -15.0 + 30.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
                
                while peakpofr <= rand1: 
                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = -15.0 + 30.0 * rand1
                    peakpofr = fit_func(rand2, 3.8) 
                
                comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']
    
                if jj == 1: 
                    orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
                    orig_wfc[i_now_peak] = 0.25
                else: 
                    append_ra.append(orig_ra[upw_index[i_now_peak]]) 
                    append_dec.append(orig_dec[upw_index[i_now_peak]]) 
                    append_z.append(comdis2z(comdis_upw + rand2, **cosmo))
                    append_wfc.append(0.25)
                    append_wcomp.append(orig_wcomp[i_now_peak])
                    
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
        orig_ra = np.append( orig_ra, np.array(append_ra) ) 
        orig_dec = np.append( orig_dec, np.array(append_dec) ) 
        orig_z = np.append( orig_z, np.array(append_z) ) 
        orig_wfc = np.append( orig_wfc, np.array(append_wfc) ) 
        orig_wcomp = np.append( orig_wcomp, np.array(append_wcomp) ) 

    elif correction['name'].lower() in ('scratch_dividedw_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 0.7
            orig_wfc[upw_index[i_now]] -= 0.7
        
    # save to file 
    scratch_file = get_galaxy_data_file('data', **cat_corr) 
    np.savetxt(scratch_file, 
            np.c_[orig_ra, orig_dec, orig_z, orig_wfc ], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_ldg_nz_down(DorR, **cat_corr): 
    ''' Downsample LasDamas mocks in order to give the nbar(z) redshift dependence
    similar to CMASS sample 

    Parameters 
    ----------
    * cat_corr : catalog correction dictionary

    
    Notes
    -----
    * simultaneously generates nbar(z) file for new downsampled nbar(z) 
    * rather than using the actual CMASS sample, Nseries mock nbar(z) ratio is used 
    * since the main goal is to mostly to construct a redshif tdependence, 
    the details are spared
    * wcp are NOT assigned

    '''
    catalog = cat_corr['catalog']
    ldg_catalog = catalog.copy()
    ldg_catalog['name'] = 'lasdamasgeo'
    true_cat_corr = {'catalog': ldg_catalog, 'correction': {'name': 'true'}} # only true can go her 

    if catalog['name'].lower() != 'ldgdownnz': 
        raise NameError('kasdkfjlaskjdf') 
    
    ldg_nbar = 0.0000944233         # default LDG nbar(z) 
    
    # read in LDG file 
    ldg_file = get_galaxy_data_file(DorR, **true_cat_corr) 
    if DorR == 'data': 
        ra, dec, z, wcp = np.loadtxt(ldg_file, unpack=True, usecols=[0,1,2,3])
    elif DorR == 'random': 
        ra, dec, cz = np.loadtxt(ldg_file, unpack=True, usecols=[0,1,2])
        z =cz/299800.0
    ngal_tot = len(z) 

    # read in Nseries nbar(z) ratio file
    nseries_cat = {'name': 'nseries'} 
    nseries_cat_corr = {'catalog': nseries_cat, 'correction': {'name': 'true'}} 

    nbar_file = fc_nbar.get_nbar_file(**nseries_cat_corr) # avg nbar file
    ratio_nbar_file = ''.join([ '/'.join( nbar_file.split('/')[:-1] ), '/', 
        'ratio_', (nbar_file.split('/')[-1])]) 
    print ratio_nbar_file
    
    nseries_zmid, nseries_zlow, nseries_zhigh, nseries_fraction = np.loadtxt(
            ratio_nbar_file, unpack=True, usecols=[0,1,2,3]) 

    zmid = nseries_zmid 
    zlow = nseries_zlow 
    zhigh = nseries_zhigh

    nseries_zmid = nseries_zmid - 0.27  # shift nseries redshift down 
    nseries_zlow = nseries_zlow - 0.27  # shift nseries redshift down 
    nseries_zhigh = nseries_zhigh - 0.27  # shift nseries redshift down 
    
    remove_indices = [] 
    ngal_remove = 0
    if DorR == 'random':
        new_nz = [] 

    # loop through redshift bins and downsample 
    for i_z, n_zmid in enumerate(nseries_zmid): 

        z_bin = np.where( (z >= nseries_zlow[i_z]) & (z < nseries_zhigh[i_z]) ) 

        ngal_bin = len(z[z_bin]) # Ngal for redshift bin 
        
        bin_down_frac = nseries_fraction[i_z] 
    
        if DorR == 'random':    # make downsampled nbar(z) file
            if n_zmid > 0.0: 
                new_nz.append(ldg_nbar * bin_down_frac) 
            
        # number of galaxies to remove 
        ngal_bin_remove = np.float(ngal_bin) * (1.0 - bin_down_frac) 
    
        if ngal_bin_remove > 0:  
            # indices to remove
            ngal_remove += int(np.rint(ngal_bin_remove))
            remove_indices += random.sample( z_bin[0], int(np.rint(ngal_bin_remove)) )
    
    if DorR == 'random': 
        new_nz += [0.0 for i in range(len(zmid) - len(new_nz))]

        new_nz_file = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
            'nbar-lasdamasgeo.down_nz.dat']) 
        print new_nz_file
        np.savetxt(new_nz_file,
                np.c_[zmid, zlow, zhigh, new_nz],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e'], delimiter='\t')

    print ngal_remove, ' Galaxies will be removed from mock realization'
    print 100.0 * np.float(ngal_remove)/np.float(ngal_tot), '% of galaxies will be removed'
    
    remain_indices = list(set(range(ngal_tot)) - set(remove_indices)) # remove 'remove indices'

    if len(remain_indices) + len(remove_indices) == ngal_tot: 
        pass
    else: 
        raise TypeError('Remove indices not working') 
    
    ldg_down_nz_file = get_galaxy_data_file(DorR, **cat_corr) 
    if DorR == 'data': 
        np.savetxt(ldg_down_nz_file, 
                np.c_[
                    ra[remain_indices], dec[remain_indices], 
                    z[remain_indices], wcp[remain_indices]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 
    elif DorR == 'random': 
        np.savetxt(ldg_down_nz_file, 
                np.c_[
                    ra[remain_indices], dec[remain_indices], 
                    z[remain_indices]
                    ], 
                fmt=['%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    return 

def build_ldg_bigfc(**cat_corr): 
    ''' Build LasDamas realizations with bigger fiber collision angular scale
    to match fiber collided fraction

    Parameters 
    ----------
    * cat_corr : catalog correction dictionary

    Notes
    -----
    * Calls IDL code ldg_bigfibcoll_wcp_assign.pro 
    * theta_ang : 

    '''
    catalog = cat_corr['catalog']
    
    # IDL command 
    fibcollided_cmd = 'idl -e "ldg_bigfibcoll_wcp_assign,'+str(catalog['n_mock'])+", '"+\
            str(catalog['letter'])+"'"+'"'
    print fibcollided_cmd
    os.system(fibcollided_cmd)  # call IDL code 

    return fibcollided_cmd

def build_ldgdownnz_bigfc(**cat_corr): 
    ''' Build LasDamas realizations with bigger fiber collision angular scale
    to match fiber collided fraction

    Parameters 
    ----------
    * cat_corr : catalog correction dictionary

    Notes
    -----
    * Calls IDL code ldg_bigfibcoll_wcp_assign.pro 
    * theta_ang : 

    '''
    catalog = cat_corr['catalog']
    
    # IDL command 
    fibcollided_cmd = 'idl -e "ldgdownnz_bigfc_wcp_assign,'+str(catalog['n_mock'])+", '"+\
            str(catalog['letter'])+"'"+'"'
    print fibcollided_cmd
    os.system(fibcollided_cmd)  # call IDL code 

    return fibcollided_cmd

def build_nseries_scratch(**cat_corr): 
    ''' Quick function to test fiber collision correction methods on Nseries mocks
        

    Parameters
    ----------
    cat_corr : catalog and correction dictionary

    Notes
    -----
    * scratch_peak_ang : Instead of using the

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
        
    omega_m = 0.31

    cosmo = {} 
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.676
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        
    # read rdzw file 
    data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
    orig_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.rdzwc']) 
    orig_ra, orig_dec, orig_z, orig_wfc, z_upw, upw_index = np.loadtxt(orig_file, unpack=True, usecols=[0,1,2,4,5,6], 
            dtype={'names': ('ra', 'dec', 'z', 'wfc', 'zupw', 'upw_index'), 
                'formats': (np.float64, np.float64, np.float64, np.float64, np.float64, np.int32)})

    # file with completeness
    mask_file = ''.join([data_dir, 'CutskyN', str(catalog['n_mock']), '.mask_info']) 
    orig_wcomp = np.loadtxt(mask_file, unpack=True, usecols=[0]) 

    if correction['name'].lower() in ('scratch_peakknown'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        
        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        #print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 
            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

        scratch_file = get_galaxy_data_file('data', **cat_corr) 
        np.savetxt(scratch_file, 
                np.c_[orig_ra, orig_dec, orig_z, orig_wfc, orig_wcomp], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

    elif correction['name'].lower() in ('scratch_peakknown_ang'):

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc 
        
        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        
        np.random.shuffle(dLOS[peak_index])     # shuffle the dLOS in the peak
        shuffle_dlos = dLOS[peak_index]
        
        #print len(shuffle_dlos), np.float(len(shuffle_dlos))/np.float(len(orig_z[now_index]))

        for j, i_now_peak in enumerate(now_peak_index): 

            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + shuffle_dlos[j], **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0

    elif correction['name'].lower() in ('scratch_peakknown_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -15.0 + 30.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -15.0 + 30.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']

            orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now_peak] = 1.0
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
    elif correction['name'].lower() in ('scratch_allpeak_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 1.0
            orig_wfc[upw_index[i_now]] -= 1.0
    
    elif correction['name'].lower() in ('scratch_peakknown_gauss_divide'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        
        now_Dc = cosmos.distance.comoving_distance(orig_z[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h
        upw_Dc = cosmos.distance.comoving_distance(z_upw[now_index], **cosmo)*cosmo['h']  # in units of Mpc/h

        dLOS = now_Dc - upw_Dc  # dLOS values for the pairs 

        peak_index = np.where(np.abs(dLOS) < 15)        # peak of the dLOS distribution
        now_peak_index = (now_index[0])[peak_index] 

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 
    
        append_ra, append_dec, append_z, append_wfc, append_wcomp = [], [], [], [], [] 
        for j, i_now_peak in enumerate(now_peak_index): 

            orig_ra[i_now_peak] = orig_ra[upw_index[i_now_peak]]
            orig_dec[i_now_peak] = orig_dec[upw_index[i_now_peak]]

            for jj in range(1,5): 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -15.0 + 30.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
                
                while peakpofr <= rand1: 
                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = -15.0 + 30.0 * rand1
                    peakpofr = fit_func(rand2, 3.8) 
                
                comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now_peak], **cosmo)*cosmo['h']
    
                if jj == 1: 
                    orig_z[i_now_peak] = comdis2z(comdis_upw + rand2, **cosmo)
                    orig_wfc[i_now_peak] = 0.25
                else: 
                    append_ra.append(orig_ra[upw_index[i_now_peak]]) 
                    append_dec.append(orig_dec[upw_index[i_now_peak]]) 
                    append_z.append(comdis2z(comdis_upw + rand2, **cosmo))
                    append_wfc.append(0.25)
                    append_wcomp.append(orig_wcomp[i_now_peak])
                    
            orig_wfc[upw_index[i_now_peak]] -= 1.0
        
        orig_ra = np.append( orig_ra, np.array(append_ra) ) 
        orig_dec = np.append( orig_dec, np.array(append_dec) ) 
        orig_z = np.append( orig_z, np.array(append_z) ) 
        orig_wfc = np.append( orig_wfc, np.array(append_wfc) ) 
        orig_wcomp = np.append( orig_wcomp, np.array(append_wcomp) ) 

    elif correction['name'].lower() in ('scratch_dividedw_gauss'): 

        now_index = np.where(orig_wfc == 0)   # galaxies with w_fc = 0 
        now_indices = (now_index[0])

        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)     # gaussian fit to dLOS distribution 

        for j, i_now in enumerate(now_indices): 

            orig_ra[i_now] = orig_ra[upw_index[i_now]]
            orig_dec[i_now] = orig_dec[upw_index[i_now]]

            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = -12.0 + 24.0 * rand1
            peakpofr = fit_func(rand2, 3.8) 
            
            while peakpofr <= rand1: 
                rand1 = np.random.random(1) 
                rand2 = np.random.random(1) 

                rand2 = -12.0 + 24.0 * rand1
                peakpofr = fit_func(rand2, 3.8) 
            
            comdis_upw = cosmos.distance.comoving_distance(z_upw[i_now], **cosmo)*cosmo['h']
            
            if z_upw[i_now] < 0.: 
                print z_upw[i_now], comdis_upw, rand2, comdis_upw + rand2
                continue

            orig_z[i_now] = comdis2z(comdis_upw + rand2, **cosmo)
            orig_wfc[i_now] = 0.7
            orig_wfc[upw_index[i_now]] -= 0.7
        
    # save to file 
    scratch_file = get_galaxy_data_file('data', **cat_corr) 
    np.savetxt(scratch_file, 
            np.c_[orig_ra, orig_dec, orig_z, orig_wfc, orig_wcomp], 
            fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

