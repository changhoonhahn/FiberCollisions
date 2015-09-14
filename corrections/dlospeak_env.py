'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution and galaxy environment.


'''
import numpy as np
import scipy as sp 
import cosmolopy as cosmos
import time 

# --- Local ---
from util.direc import direc
from util.catalog import Catalog
from corrections import Corrections
from fibcollided import UpweightCorr 
from galenv.galenv import d_NN_dataclass

class DlospeakEnvCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        """ Child class of Correction class in corrections.py 
        Fiber collisions correction using the peak of the line-of-sight displacement 
        distribution and galaxy environment. Depending on the environment of the
        upweighted galaxy of the fiber collided pair, we sample dLOS from different
        gaussians

        --------------------------------------------------------------------------
        Notes
        --------------------------------------------------------------------------
        * Currently supported peak correction methods: peakshot 
        * nbar(z) interpolation implemented for CMASS like samples; however this can easily 
        be extended to other mocks
        * dLOS within peak is sampled +/- 3-sigmas

        """

        super(DlospeakEnvCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        """ Specify correction string for dlos peak correction 
        """
        corr = (self.cat_corr)['correction']

        if not all(x in corr.keys() for x in ['fit', 'n_NN', 'fpeak', 'sigma']): 
            raise KeyError("Specify fpeak, sigma, fit, and n_NN in correction dictionary")

        self.corr_str = ''.join([
            '.', corr['fit'].lower(), 
            '.', str(corr['n_NN']), 'NN',
            '.', corr['name'].lower(),
            '.sigma', str(corr['sigma']), '.fpeak', str(corr['fpeak'])
            ])
        return self.corr_str

    def build(self): 
        """ Build fiber collision corrected mock catalogs. Using the modeled dLOS peak 
        and galaxy environment of upweighted galaxies. Note that this method reads in 
        upweight corrected mock catalogs and if care is not given, may result in circular
        dependencies. 
        """

        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']

        cosmo = self.cosmo()      # cosmoslogy 

        f_peak = corrdict['fpeak']

        # survey redshift limits  
        survey_zmin, survey_zmax = self.survey_zlimits()  

        survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
        # spline interpolation function hardcoded here 
        # to make it faster
        z_arr = np.arange(0.0, 1.01, 0.01)
        dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo) * cosmo['h']
        comdis2z = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 

        # upweight corrected galaxy catalog is used to 
        # in the correction. This is out of convenience and 
        # also due to the fact that the BOSS observed
        # galaxy catalog outputs an upweigh corrected galaxy 
        # catalog. 
        fc_cat_corr = {
                'catalog': (self.cat_corr)['catalog'], 
                'correction': {'name': 'upweight'}
                }
        fc_mock = UpweightCorr(fc_cat_corr, **self.kwargs) 
        fc_mock_file = fc_mock.file()
        fc_mock_cols = fc_mock.datacolumns()

        fc_data = np.loadtxt(
                fc_mock_file, 
                skiprows=1, 
                unpack=True, 
                usecols=range(len(fc_mock_cols))
                )
        
        for i_col, fc_col in enumerate(fc_mock_cols): 
            setattr(fc_mock, fc_col, fc_data[i_col])

        upw = np.where(fc_mock.wfc > 1)  # upweighted

        notcoll = np.where(fc_mock.wfc > 0)  # not collided
        
        # expected number of galaxies to be placed in the peak 
        # of the dLOS distribution calculated based on the peak 
        # fraction. 
        ngal_peak_exp = int( 
                np.floor(
                    f_peak * (np.sum(fc_data.wfc[upw]) - np.float(len(upw[0])))
                    ) 
                )
        
        # nbar(z) cubic spline interpolation, which is currently hardcoded
        # in the function temp_nbarz(). However, this should 
        # ultimately be edited so that it works for any catalog 
        # (currently only works for CMASS-like catalogs)
        if 'cmass' in catdict['name'].lower(): 

            nb_z, nb_nbar = temp_nbarz(self.cat_corr)

            nbarofz = sp.interpolate.interp1d(nb_z, nb_nbar, kind='cubic')       

        # calculate galaxy environment for upweighted galaxies
        # in fiber collided pairs using d_NN_dataclass function in 
        # galenv module 
        start_time = time.time()
        dNNs = d_NN_dataclass(
                fc_data.ra[upw],
                fc_data.dec[upw],
                fc_data.z[upw], 
                fc_data
                )
        print 'd_NN calculation takes ', (time.time()-start_time)/60.0, ' minutes'

        # fpeak and sigma values for dNN values 
        fpeak_envs = dlosenv_fpeak_env(dNNs, n_NN = corrdict['n_NN'])
        sigma_envs = dlosenv_sigma_env(dNNs, n_NN = corrdict['n_NN'])

        append_i, append_z, append_nbar, sampled_dlos = [], [], [], [] 
        
        for ii_gal, i_gal in enumerate(upw[0]):    # for each upweighted galaxy 
            
            # fiber collision weight is in place because
            # we need to make sure that the dLOS peak is 
            # sampled for every upweight
            wfc_counter = fc_mock.wfc[i_gal]

            while wfc_counter > 1:
            
                wfc_counter -= 1.0 

                fpeak_env = fpeak_envs[ii_gal]
                sigma_env = sigma_envs[ii_gal]

                rand_fpeak = np.random.random(1) 
                
                if rand_fpeak <= fpeak_env:        # in the peak 
                
                    # downweight upweighted galaxy 
                    fc_mock.wfc[i_gal] -= 1.0

                    append_i.append(i_gal)

                    comdis_igal = cosmos.distance.comoving_distance(fc_mock.z[i_gal], **cosmo) * cosmo['h']
                
                    # sampled dLOS from dLOS peak distribution best-fit function 
                    d_samp = sample_dlos_peak(corrdict['fit'], sigma_env)

                    # in case the displacement falls out of bounds
                    # of the survey redshift limits. This crude treatment
                    # may generate large scale issues. But currently 
                    # ignoring this issue 
                    if (comdis_igal + d_samp > survey_comdis_max) or (comdis_igal + d_samp < survey_comdis_min): 
                        d_samp = -1.0 * d_samp 
                    
                    # convert comoving distance to redshift 
                    collided_z = comdis2z(comdis_igal + d_samp)
                    
                    append_z.append(collided_z[0]) 
                    
                    if 'cmass' in catdict['name'].lower():
                        append_nbar.append(nbarofz(collided_z[0]))

                    sampled_dlos.append(d_samp) 

                    ngal_peak_exp -= 1.0 
         
        n_append = len(append_i)
        print n_append, ' Galaxies were peak corrected'
        print ngal_peak_exp, ' extra galaxies need to be in the peak' 
        
        """
        data_cols = self.datacolumns()
        data_fmts = self.datacols_fmt()
        data_hdrs = self.datacols_header()
        data_list = []  

        for data_col in data_cols: 

            if data_col not in ('wfc', 'wnoz', 'z', 'nbar'):
                append_arr = getattr(fc_mock, data_col)[append_i]

            elif data_col == 'z':
                append_arr = append_z 

            elif data_col == 'nbar': 
                append_arr = append_nbar

            else: # new galaxies have wfc and wnoz = 1.0 
                append_arr = np.array([1.0 for i in range(n_append)])

            new_col = np.concatenate([
                getattr(fc_mock, data_col)[notcoll], 
                append_arr
                ])

            data_list.append(new_col)

        # write to corrected data to file 
        output_file = self.file()
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmts, 
                delimiter='\t', 
                header=data_hdrs
                ) 
        
        # write sampled dLOS value to match to actual dLOS peak 
        np.savetxt(
                output_file+'.dlos', 
                np.c_[sampled_dlos], 
                fmt=['%10.5f'], 
                delimiter='\t'
                ) 
        """

def temp_nbarz(cat_corr):
    """ nbar(z) data for given catalog and correction. Temporary function. 
    """

    catdict = cat_corr['catalog']
    catalog_name = catdict['name'].lower()
    
    nbar_dir = direc('data', cat_corr)

    if 'cmass' in catalog_name: 

        if catalog_name == 'cmass': 
            nbar_file = ''.join([
                nbar_dir, 
                'nbar-cmass-dr12v4-N-Reid-om0p31_Pfkp10000.dat'
                ])

        elif 'cmasslowz' in catalog_name: 

            if 'e2' in catalog_name: 
                nbar_file = ''.join([
                    nbar_dir, 
                    'nbar_DR12v5_CMASSLOWZE2_North_om0p31_Pfkp10000.dat'
                    ])
            elif 'e3' in catalog_name: 
                nbar_file = ''.join([
                    nbar_dir, 
                    'nbar_DR12v5_CMASSLOWZE3_North_om0p31_Pfkp10000.dat'
                    ])
            else: 
                nbar_file = ''.join([
                    nbar_dir, 
                    'nbar_DR12v5_CMASSLOWZ_North_om0p31_Pfkp10000.dat'
                    ])
    else: 
        raise NotImplementedError()

    # read nbar(z) file 
    nb_z, nb_nbar = np.loadtxt(nbar_file, skiprows=2, unpack=True, usecols=[0, 3]) 
    return [nb_z, nb_nbar]

def sample_dlos_peak(fit, sigma): 
    """ Randomly sample a line-of-sight displacement value from the peak 
    of the distribution described by either a best-fit guassian, 
    exponential, or the true distribution of the peak itself. 
    """

    if fit in ('gauss', 'expon'):   
        # sample dLOS from peak using best-fit

        if fit == 'gauss': 
            peak_fit_func = peak_fit_gauss
        elif fit == 'expon': 
            peak_fit_func = peak_fit_expon

        rand1 = np.random.random(1) 
        rand2 = np.random.random(1) 

        rand2 = (-3.0 + rand2 * 6.0) * sigma 

        peakpofr = peak_fit_func(rand2, sigma) 
        
        while peakpofr <= rand1: 
            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = (-3.0 + rand2 * 6.0) * sigma
            peakpofr = peak_fit_func(rand2, sigma) 

    else: 
        raise NotImplementedError('asdfasdf')

    return rand2 

def peak_fit_gauss(x, sig): 
    """ Gaussian function  
    """
    return np.exp(-0.5 * x**2/sig**2)

def peak_fit_expon(x, sig): 
    """ Exponential function 
    """
    return np.exp(-1.0*np.abs(x)/sig)

def dlosenv_fpeak_env(dNN, n_NN=3): 
    """ Calculate fpeak(dNN) using the best-fit function from the 
    dlosenv_peakfit_fpeak_env_fit function in dlos/fitting.py moduel.
    """
    
    # read in bestfit slope and y-int from file 
    dlosclass = Dlos(cat_corr)
    dlos_fpeak_envbin_fit_file = ''.join([
        (dlosclass.file_name).rsplit('/', 1)[0], '/', 
        'DLOS_fpeak_env_d', str(n_NN), 'NN_bin_bestfit.dat'
        ])

    bestfit_slope, bestfit_yint = np.loadtxt(
            dlos_fpeak_envbin_fit_file, 
            skiprows = 1, 
            unpack = True,
            usecols = [0, 1]
            ) 

    return bestfit_slope * dNN + bestfit_yint 

def dlosenv_sigma_env(dNN, n_NN=3): 
    """ Calculate fpeak(dNN) using the best-fit function from the 
    dlosenv_peakfit_fpeak_env_fit function in dlos/fitting.py moduel.
    """
    
    # read in bestfit slope and y-int from file 
    dlosclass = Dlos(cat_corr)
    dlos_fpeak_envbin_fit_file = ''.join([
        (dlosclass.file_name).rsplit('/', 1)[0], '/', 
        'DLOS_sigma_env_d', str(n_NN), 'NN_bin_bestfit.dat'
        ])

    bestfit_slope, bestfit_yint = np.loadtxt(
            dlos_fpeak_envbin_fit_file, 
            skiprows = 1, 
            unpack = True,
            usecols = [0, 1]
            ) 

    return bestfit_slope * dNN + bestfit_yint 

"""
                    elif corrdict['fit'].lower() == 'true': 
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
"""
