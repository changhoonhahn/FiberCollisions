'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution. Code could use overall editing


'''

import cosmolopy as cosmos

# --- Local ---
from util import util

def file(cat_corr, **kwargs): 
    """ Specify correction string
    """
    cat = cat_corr['catalog']
    corr = cat_corr['correction']

    if not all(x in corr.keys() for x in ['fpeak', 'sigma', 'fit']): 
        raise KeyError("Specify fpeak, sigma, and fit in correction dictionary")

    corr_str = ''.join([
        '.', corr['fit'].lower(), '.', corr['name'].lower(), 
        '.sigma', str(corr['sigma']), '.fpeak', str(corr['fpeak'])
        ])
    
    return corr_str

def build(cat_corr, **kwargs): 
    ''' Build peak corrected fibercollided mock catalogs (using cosmolopy). 
    dLOS peak corrected mock catalogs are constructed from fibercollided mock catalogs

    Parameters
    ----------
    cat_corr : Catalog + Correction dictionary

    Notes
    -----
    * Currently supported peak correction methods: peakshot 
    * nbar(z) interpolation implemented for CMASS like samples; however this can easily 
    be extended to other mocks
    * dLOS within peak is sampled +/- 3-sigmas

    '''

    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    f_peak = correction['fpeak']    # peak fraction 

    # fit functions
    if correction['fit'].lower() == 'gauss': 
        fit_func = lambda x, sig: np.exp(-0.5 *x**2/sig**2)
    elif correction['fit'].lower() == 'expon': 
        fit_func = lambda x, sig: np.exp(-1.0*np.abs(x)/sig)
    elif correction['fit'].lower() == 'true': 
        pass 
    else: 
        raise NameError('correction fit has to be specified as gauss or expon') 
   
    # upweight corrected galaxy catalog 
    fc_cat_corr = {'catalog': catalog, 'correction': {'name': 'upweight'}}
    fc_mock = spec_gal.Data('data', fc_cat_corr, **kwargs) 
    fc_mock.Read()
    cosmo = fc_mock.Cosmo()      # cosmoslogy 

    if catalog['name'].lower() not in ('tilingmock', 'lasdamasgeo', 'ldgdownnz'):
        # blah
        fc_mock.weight = fc_mock.wfc            

    # survey redshift limits 
    if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):     
        survey_zmin, survey_zmax = 0.16, 0.44
    elif catalog['name'].lower() in ('tilingmock', 'qpm', 'patchy', 'nseries'): 
        survey_zmin, survey_zmax = 0.43, 0.7 
    elif catalog['name'].lower() in ('bigmd'):             
        survey_zmin, survey_zmax = 0.43, 0.7    
    elif 'cmass' in catalog['name'].lower():             
        if catalog['name'].lower() == 'cmass': 
            survey_zmin, survey_zmax = 0.43, 0.7  
        elif 'cmasslowz' in catalog['name'].lower(): 
            if '_high' in catalog['name'].lower(): 
                survey_zmin, survey_zmax = 0.5, 0.75    
            elif '_low' in catalog['name'].lower(): 
                survey_zmin, survey_zmax = 0.2, 0.5
        else: 
            raise NotImplementedError('CMASS or CMASSLOWZ combined sample')
    else: 
        raise NotImplementedError('Mock Catalog not included')
    
    # comoving distnace of z_min and z_max 
    survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
    survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
    # nbar(z) for interpolated nbar(z) (hardcoded for CMASS-like catalogs)
    if 'cmass' in catalog['name'].lower(): 
        nb_z, nb_nbar = temp_nbarz(cat_corr)
        # cubic spline nbar(z) interpolation
        nbarofz = sp.interpolate.interp1d(nb_z, nb_nbar, kind='cubic')       

    appended_ra, appended_dec, appended_z, appended_weight = [], [], [], []
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        appended_comp = []   
    elif 'cmass' in catalog['name'].lower(): 
        appended_comp, appended_wsys, appended_wnoz, appended_nbar = [], [], [], [] 

    sampled_dlos = [] 
    
    upw = np.where(fc_mock.weight > 1)  # upweighted
    notcoll = np.where(fc_mock.weight > 0)  # not collided

    for i_gal in upw[0]:    # for every upweighted galaxy 

        while fc_mock.weight[i_gal] > 1:

            rand_fpeak = np.random.random(1) 
            if rand_fpeak <= f_peak:          
                # in the peak 
                
                # downweight upweighted galaxy 
                fc_mock.weight[i_gal] -= 1.0

                # LOS comoving distance of the upweighted galaxy 
                comdis_igal = \
                        cosmos.distance.comoving_distance(
                                fc_mock.z[i_gal], **cosmo) * cosmo['h']

                appended_ra.append(fc_mock.ra[i_gal])     
                appended_dec.append(fc_mock.dec[i_gal])
                appended_weight.append(1.0)     # "new" galaxy has wfc=1

                if catalog['name'].lower() in ('qpm', 'nseries'):  
                    appended_comp.append(fc_mock.comp[i_gal]) 
                elif 'cmass' in catalog['name'].lower(): 
                    appended_comp.append(fc_mock.comp[i_gal]) 
                    appended_wsys.append(fc_mock.wsys[i_gal]) 
                    appended_wnoz.append(1.0) 

                if correction['fit'].lower() in ('gauss', 'expon'):   
                    # sample dLOS from peak using best-fit

                    rand1 = np.random.random(1) 
                    rand2 = np.random.random(1) 

                    rand2 = (-3.0 + rand2 * 6.0)*correction['sigma']
                    peakpofr = fit_func(rand2, correction['sigma']) 
                    
                    while peakpofr <= rand1: 
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*correction['sigma']
                        peakpofr = fit_func(rand2, correction['sigma']) 

                elif correction['fit'].lower() == 'true': 
                    raise NotImplementedError("Need to revive")
                    '''
                    # sample dLOS within peak from actual distribution   
                    dlos_comb_peak_file = ''.join([
                        ((fc_mock.file_name).rsplit('/', 1))[0], '/', 
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
                    '''
                else: 
                    raise NotImplementedError('asdfasdf')

                # in case the displacement falls out of bound (may general large scale issues)
                if (comdis_igal + rand2 > survey_comdis_max) or (comdis_igal + rand2 < survey_comdis_min): 
                    rand2 = -1.0 * rand2
                
                # convert comoving distance to redshift 
                collided_z = util.comdis2z(comdis_igal + rand2, **cosmo)
                
                appended_z.append(collided_z[0]) 

                if 'cmass' in catalog['name'].lower():
                    appended_nbar.append(nbarofz(collided_z[0]))

                sampled_dlos.append(rand2) 

    print len(appended_ra), ' galaxies were peak corrected'
    # append artificial galaxies to catalog 
    fc_mock.ra = np.concatenate([fc_mock.ra[notcoll], appended_ra])
    fc_mock.dec = np.concatenate([fc_mock.dec[notcoll], appended_dec])
    fc_mock.weight = np.concatenate([fc_mock.weight[notcoll], appended_weight])
    fc_mock.z = np.concatenate([fc_mock.z[notcoll], appended_z])
        
    if catalog['name'].lower() in ('qpm', 'nseries'): 
        fc_mock.comp = np.concatenate([fc_mock.comp[notcoll], appended_comp])
    elif 'cmass' in catalog['name'].lower():
        fc_mock.comp = np.concatenate([fc_mock.comp[notcoll], appended_comp])
        fc_mock.wsys = np.concatenate([fc_mock.wsys[notcoll], appended_wsys])
        fc_mock.wnoz = np.concatenate([fc_mock.wnoz[notcoll], appended_wnoz])
        fc_mock.nbar = np.concatenate([fc_mock.nbar[notcoll], appended_nbar])
    
    # write peak corrected data to file 
    if catalog['name'].lower() in ('qpm', 'nseries'): 

        header_str = "Column : ra, dec, z, wfc, comp"
        data_list = [fc_mock.ra, fc_mock.dec, fc_mock.z, fc_mock.weight, fc_mock.comp]
        data_fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 

    elif 'cmass' in catalog['name'].lower():          # CAMSS

        header_str = "Columns : ra, dec, z, nbar, w_systot, w_noz, w_fc, comp"
        data_list = [fc_mock.ra, fc_mock.dec, fc_mock.z, fc_mock.nbar,
                fc_mock.wsys, fc_mock.wnoz, fc_mock.weight, fc_mock.comp]
        data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f', '%10.5f', '%10.5f', '%10.5f']

    elif catalog['name'].lower() in ('tilingmock', 'lasdamasgeo', 'ldgdownnz'): 

        header_str = "Columns : ra, dec, z, wfc" 
        data_list = [fc_mock.ra, fc_mock.dec, fc_mock.z, fc_mock.weight]
        data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']

    else: 
        raise NotImplementedError('asdfasdf')

    # write to corrected file 
    corr_class = Corr(cat_corr) 
    output_file = corr_class.file()
    np.savetxt(output_file, (np.vstack(np.array(data_list))).T, 
            fmt=data_fmt, delimiter='\t', header=header_str) 
        
    np.savetxt(output_file+'.dlos', np.c_[sampled_dlos], fmt=['%10.5f'], delimiter='\t') 

def temp_nbarz(cat_corr):
    """ nbar(z) data for given catalog and correction. Temporary function. 
    """
    catalog = cat_corr['catalog']
    
    if 'cmass' in catalog['name'].lower(): 
        if catalog['name'].lower() == 'cmass': 
            nbar_dir = '/mount/riachuelo1/hahn/data/CMASS/'
            nbar_file = ''.join([nbar_dir, 'nbar-cmass-dr12v4-N-Reid-om0p31_Pfkp10000.dat'])

        elif 'cmasslowz' in catalog['name'].lower(): 
            nbar_dir = '/mount/riachuelo1/hahn/data/CMASS/dr12v5/'
            if 'e2' in catalog['name'].lower(): 
                nbar_file = ''.join([nbar_dir, 
                    'nbar_DR12v5_CMASSLOWZE2_North_om0p31_Pfkp10000.dat'])
            elif 'e3' in catalog['name'].lower(): 
                nbar_file = ''.join([nbar_dir, 
                    'nbar_DR12v5_CMASSLOWZE3_North_om0p31_Pfkp10000.dat'])
            else: 
                nbar_file = ''.join([nbar_dir, 
                    'nbar_DR12v5_CMASSLOWZ_North_om0p31_Pfkp10000.dat'])
    else: 
        raise NotImplementedError()

    # read nbar(z) file 
    nb_z, nb_nbar = np.loadtxt(nbar_file, skiprows=2, unpack=True, usecols=[0, 3]) 
    return [nb_z, nb_nbar]
