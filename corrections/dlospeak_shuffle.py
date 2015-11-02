'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution. Code could use overall editing


'''
import numpy as np
import scipy as sp 
import cosmolopy as cosmos
from scipy.stats import norm
from scipy.stats import expon
from scipy.stats import binom
import time 

# --- Local ---
from util.direc import direc
from util.catalog import Catalog
from corrections import Corrections
from fibcollided import UpweightCorr 

class DlospeakShuffleCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        """ Child class of Correction class in corrections.py 
        Fiber collisions correction using the peak of the line-of-sight displacement 
        distribution 

        Notes
        -----
        * Currently supported peak correction methods: peakshot 
        * dLOS within peak is sampled +/- 3-sigmas
        """

        super(DlospeakShuffleCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        """ Specify correction string for dlos peak correction 
        """
        corr = (self.cat_corr)['correction']

        self.corr_str = ''.join([
            '.', corr['name'].lower(), 
            '.sigma', str(corr['sigma'])
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

        # survey redshift limits  
        survey_zmin, survey_zmax = self.survey_zlimits()  

        survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
        # spline interpolation function hardcoded here 
        # to make it faster
        z_arr = np.arange(0.0, 1.01, 0.01)
        dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo) * cosmo['h']
        comdis2z = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
        
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

        for i_col, fc_col in enumerate(fc_mock_cols): 
           setattr(fc_mock, fc_col, fc_data[i_col])

        fc_mock.upw_index = fc_mock.upw_index.astype(int)

        collided = np.where(fc_mock.wfc == 0)[0]   # collided
        n_fcpair = len(collided)     # number of fiber collision pairs 

        # D_comov(z_upw), D_comvo(z_coll)
        comdis_upw_gal = cosmos.distance.comoving_distance(
                fc_mock.zupw[collided], **cosmo) * cosmo['h']
        comdis_coll_gal = cosmos.distance.comoving_distance(
                fc_mock.z[collided], **cosmo) * cosmo['h']

        dlos_actual = comdis_coll_gal - comdis_upw_gal   # line-of-sight displacement

        within_peak = np.where(np.abs(dlos_actual) < 3.0 * corrdict['sigma'])

        i_peakcorr = collided[within_peak]
        d_samples = dlos_actual[within_peak]

        np.random.shuffle(d_samples)        # shuffle them up!
        
        fc_mock.wfc[i_peakcorr] += 1.0
        for i_upw in fc_mock.upw_index[i_peakcorr]: 
            fc_mock.wfc[i_upw] -= 1.0

        #fc_mock.wfc[i_tailcorr] += 1.0
        #fc_mock.wfc[fc_mock.upw_index[i_tailcorr]] -= 1.0
            
        #fc_mock.ra[i_peakcorr] = fc_mock.ra[fc_mock.upw_index[i_peakcorr]]
        #fc_mock.dec[i_peakcorr] = fc_mock.dec[fc_mock.upw_index[i_peakcorr]]
            
        # account for peak correction that places galaxies out of bound. 
        #outofbounds = np.where( 
        #        (comdis_upw_gal[within_peak] + d_samples > survey_comdis_max) |
        #        (comdis_upw_gal[within_peak] + d_samples < survey_comdis_min)
        #        )

        #if len(outofbounds[0]) > 0: 
        #    d_samples[outofbounds] *= -1.0
                    
        collided_z = comdis2z(comdis_upw_gal[within_peak] + d_samples)
        fc_mock.z[i_peakcorr] = collided_z

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
