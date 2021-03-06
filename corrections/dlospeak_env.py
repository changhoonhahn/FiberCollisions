'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution and galaxy environment.


'''
import numpy as np
import scipy as sp 
from scipy.stats import norm
import cosmolopy as cosmos
import time 
import random 

# --- Local ---
from util.direc import direc
from util.catalog import Catalog
from corrections import Corrections
from fibcollided import UpweightCorr 
from galenv.galenv import d_NN_dataclass
from dlospeak import sample_dlos_peak
from dlospeak import peak_fit_gauss 
from dlospeak import peak_fit_expon 
from dlospeak import temp_nbarz 

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
        * Somehow NOT correcting for extra peak corrected galaxies improves the overal 
        power spectrum. (This has been tested)

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
        
        fc_mock.upw_index = fc_mock.upw_index.astype(int)

        upw = np.where(fc_mock.wfc > 1)  # upweighted
        collided = np.where(fc_mock.wfc == 0)[0] # collided 
        notcoll = np.where(fc_mock.wfc > 0)  # not collided

        n_fcpair = len(collided)

        # expected number of galaxies to be placed in the peak 
        # of the dLOS distribution calculated based on the peak 
        # fraction. 
        ngal_peak_exp = int(np.rint(np.float(n_fcpair) * f_peak))
        
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
                fc_mock.ra[fc_mock.upw_index[collided]],
                fc_mock.dec[fc_mock.upw_index[collided]],
                fc_mock.z[fc_mock.upw_index[collided]], 
                fc_mock
                )
        print 'd_NN calculation takes ', (time.time()-start_time)/60.0, ' minutes'

        # fpeak and sigma values for dNN values 
        fpeak_envs = dlosenv_fpeak_env(dNNs, fc_cat_corr, n_NN = corrdict['n_NN'])
        sigma_envs = dlosenv_sigma_env(dNNs, fc_cat_corr, n_NN = corrdict['n_NN'])
        print fpeak_envs
        print sigma_envs
   
        np.random.seed()
        rand = np.random.uniform(size = n_fcpair)

        peakcorr = np.where(fpeak_envs > rand)
        n_peakcorr = len(peakcorr[0]) 
        print '' 
        print 'Ngal_(expected in peak) = ', ngal_peak_exp
        print 'Ngal_(actually peak corrected) = ', n_peakcorr 
        print np.abs(np.float(ngal_peak_exp) - np.float(n_peakcorr))/np.float(ngal_peak_exp) 
        
        fc_mock.wfc[collided[peakcorr]] += 1.0
        fc_mock.wfc[fc_mock.upw_index[collided[peakcorr]]] -= 1.0

        fc_mock.ra[collided[peakcorr]] = fc_mock.ra[fc_mock.upw_index[collided[peakcorr]]]
        fc_mock.dec[collided[peakcorr]] = fc_mock.dec[fc_mock.upw_index[collided[peakcorr]]]
            
        # comoving distances of upweighted galaxies in fc pairs 
        # that are going to be peakcorrected
        comdis_upw_gal = cosmos.distance.comoving_distance(
                fc_mock.zupw[collided[peakcorr]], **cosmo) * cosmo['h']

        dlos_pdf = norm(loc = 0.0, scale = 1.0) 
        d_samples = sigma_envs[peakcorr] * dlos_pdf.rvs(n_peakcorr)
            
        # account for peak correction that places galaxies out of bound. 
        outofbounds = np.where( 
                (comdis_upw_gal + d_samples > survey_comdis_max) |
                (comdis_upw_gal + d_samples < survey_comdis_min)
                )

        if len(outofbounds[0]) > 0: 
            d_samples[outofbounds] *= -1.0
                    
        collided_z = comdis2z(comdis_upw_gal + d_samples)
        
        fc_mock.z[collided[peakcorr]] = collided_z

        data_cols = self.datacolumns()
        data_fmts = self.datacols_fmt()
        data_hdrs = self.datacols_header()
        data_list = []  
        for data_col in data_cols: 
            
            new_col = getattr(fc_mock, data_col)

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
                np.c_[d_samples], 
                fmt=['%10.5f'], 
                delimiter='\t'
                ) 

def dlosenv_fpeak_env(dNN, cat_corr, n_NN=3): 
    """ Calculate fpeak(dNN) using the best-fit function from the 
    dlosenv_peakfit_fpeak_env_fit function in dlos/fitting.py moduel.
    """
    
    # read in bestfit slope and y-int from file 
    dlos_fpeak_envbin_fit_file = ''.join([
        direc('data', cat_corr), 
        'DLOS_fpeak_env_d', str(n_NN), 'NN_bin_bestfit.dat'
        ])

    bestfit_slope, bestfit_yint = np.loadtxt(
            dlos_fpeak_envbin_fit_file, 
            skiprows = 1, 
            unpack = True,
            usecols = [0, 1]
            ) 

    return bestfit_slope * dNN + bestfit_yint 

def dlosenv_sigma_env(dNN, cat_corr, n_NN=3): 
    """ Calculate fpeak(dNN) using the best-fit function from the 
    dlosenv_peakfit_fpeak_env_fit function in dlos/fitting.py moduel.
    """
    
    # read in bestfit slope and y-int from file 
    dlos_sigma_envbin_fit_file = ''.join([
        direc('data', cat_corr), 
        'DLOS_sigma_env_d', str(n_NN), 'NN_bin_bestfit.dat'
        ])

    bestfit_slope, bestfit_yint = np.loadtxt(
            dlos_sigma_envbin_fit_file, 
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
