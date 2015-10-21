'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution. Code could use overall editing


'''
import numpy as np
import scipy as sp 
import cosmolopy as cosmos
from scipy.stats import norm
import time 

# --- Local ---
from util.direc import direc
from util.catalog import Catalog
from corrections import Corrections
from fibcollided import UpweightCorr 

class DlospeakCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        """ Child class of Correction class in corrections.py 
        Fiber collisions correction using the peak of the line-of-sight displacement 
        distribution 

        Notes
        -----
        * Currently supported peak correction methods: peakshot 
        * nbar(z) interpolation implemented for CMASS like samples; however this can easily 
        be extended to other mocks
        * dLOS within peak is sampled +/- 3-sigmas
        """

        super(DlospeakCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        """ Specify correction string for dlos peak correction 
        """
        corr = (self.cat_corr)['correction']

        if not all(x in corr.keys() for x in ['fpeak', 'sigma', 'fit']): 
            raise KeyError("Specify fpeak, sigma, and fit in correction dictionary")

        self.corr_str = ''.join([
            '.', corr['fit'].lower(), '.', corr['name'].lower(), 
            '.sigma', str(corr['sigma']), '.fpeak', str(corr['fpeak'])
            ])
        return self.corr_str

    def build(self): 
        """ Build peak corrected fibercollided mock catalogs (using cosmolopy). 
        Note dLOS peak corrected mock catalogs are constructed from fibercollided 
        mock catalogs. 

        """

        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']

        cosmo = self.cosmo()      # cosmoslogy 

        f_peak = corrdict['fpeak']    # peak fraction 
        
        # survey redshift limits  
        survey_zmin, survey_zmax = self.survey_zlimits()  

        survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
        # spline interpolation function hardcoded here 
        # to make it faster
        z_arr = np.arange(0.0, 1.01, 0.01)
        dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo) * cosmo['h']
        comdis2z = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
        
        # nbar(z) cubic spline interpolation, which is currently hardcoded
        # in the function temp_nbarz(). However, this should 
        # ultimately be edited so that it works for any catalog 
        # (currently only works for CMASS-like catalogs)
        if 'cmass' in catdict['name'].lower(): 

            nb_z, nb_nbar = temp_nbarz(self.cat_corr)

            nbarofz = sp.interpolate.interp1d(nb_z, nb_nbar, kind='cubic')       

        # upweight corrected galaxy catalog 
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
        fc_mock_col_fmts = [] 
        for col in fc_mock_cols: 
            if col == 'upw_index':
                fc_mock_col_fmts.append(np.int)
            else: 
                fc_mock_col_fmts.append(np.float)

        fc_data = np.loadtxt(
                fc_mock_file, 
                skiprows=1, 
                unpack=True, 
                usecols=range(len(fc_mock_cols)),
                dtype = {
                    'names': tuple(fc_mock_cols),
                    'formats': tuple(fc_mock_col_fmts)
                    }
                )
        
        for i_col, fc_col in enumerate(fc_mock_cols): 
            setattr(fc_mock, fc_col, fc_data[i_col])

        fc_mock.upw_index = fc_mock.upw_index.astype(int)

        upw = np.where(fc_mock.wfc > 1)         # upweighted
        collided = np.where(fc_mock.wfc == 0)[0]   # collided
        notcoll = np.where(fc_mock.wfc > 0)     # not collided

        n_fcpair = len(collided)     # number of fiber collision pairs 
        
        # expected number of galaxies to be placed in the peak 
        # of the dLOS distribution based on the peak fraction 
        n_peak_exp = int(
                np.rint(f_peak * np.float(n_fcpair))
                )
    
        np.random.seed()
        np.random.shuffle(collided)

        i_peakcorr = collided[:n_peak_exp]  # indices of peak corrected galaxies

        fc_mock.wfc[i_peakcorr] += 1.0
        fc_mock.wfc[fc_mock.upw_index[i_peakcorr]] -= 1.0
        
        fc_mock.ra[i_peakcorr] = fc_mock.ra[fc_mock.upw_index[i_peakcorr]]
        fc_mock.dec[i_peakcorr] = fc_mock.dec[fc_mock.upw_index[i_peakcorr]]
        
        # comoving distances of upweighted galaxies in fc pairs 
        # that are going to be peakcorrected
        comdis_upw_gal = cosmos.distance.comoving_distance(
                fc_mock.zupw[i_peakcorr], **cosmo) * cosmo['h']
        
        # sampled dLOS values from best-fit functional form of 
        # dLOS peak distribution
        if corrdict['fit'] == 'gauss': 
            dlos_pdf = norm(loc = 0.0, scale = corrdict['sigma'])
        else: 
            raise ValueError
        d_samples = dlos_pdf.rvs(n_peak_exp)
        
        # account for peak correction that places galaxies out of bound. 
        outofbounds = np.where( 
                (comdis_upw_gal + d_samples > survey_comdis_max) |
                (comdis_upw_gal + d_samples < survey_comdis_min)
                )

        if len(outofbounds[0]) > 0: 
            d_samples[outofbounds] *= -1.0
                    
        collided_z = comdis2z(comdis_upw_gal + d_samples)
        
        fc_mock.z[i_peakcorr] = collided_z

        if 'cmass' in catdict['name'].lower():
            fc_mock.nbar[i_peakcorr] = nbarofz(collided_z)
        
        print n_peak_exp, ' Galaxies were peak corrected'

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

    np.random.seed()
    if fit in ('gauss'): 
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
