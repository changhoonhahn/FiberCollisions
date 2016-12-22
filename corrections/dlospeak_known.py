'''

LOS Reconstruction Fiber collision correction, for the
case when correlated/chance alignment pairs are KNOWN.


'''
import numpy as np
import scipy as sp 
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
from dlospeak import temp_nbarz 

class DlospeakKnownCorr(Corrections): 

    def __init__(self, cat_corr, **kwargs): 
        ''' Child class of Correction class in corrections.py 
        Fiber collisions correction using the LOS reconstruction method. For
        this correction in particular, the correlated and chance alignment
        collided pairs are known

        '''
        super(DlospeakKnownCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        ''' Specify correction string for dlos peak correction 
        d_photoz_tail_cut has to be specified
        '''
        corr = (self.cat_corr)['correction']

        if not all(x in corr.keys() for x in ['fit', 'fpeak', 'sigma']): 
            raise KeyError("Specify fpeak, sigma, and fit in correction dictionary")

        self.corr_str = ''.join([
            '.', corr['fit'].lower(), 
            '.', corr['name'].lower(),
            '.sigma', str(corr['sigma'])
            ])
        return self.corr_str

    def build(self): 
        ''' Build fiber collision corrected mock catalogs using LOS distribution. 
        '''
        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']

        if catdict['name'] != 'nseries':
            raise ValueError('only implemented for Nseris') 

        cosmo = self.cosmo()      # cosmoslogy 

        #f_peak = corrdict['fpeak']

        # survey redshift limits  
        survey_zmin, survey_zmax = self.survey_zlimits()  
        survey_comdis_min = cosmos.distance.comoving_distance(survey_zmin, **cosmo) *\
                cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance(survey_zmax, **cosmo) *\
                cosmo['h']
    
        # spline interpolation function hardcoded here 
        # to make it faster
        z_arr = np.arange(0.0, 1.01, 0.01)
        dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo) * cosmo['h']
        comdis2z = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
        
        # read in galxies from the fibcollided mocks
        fc_cat_corr = {
                'catalog': catdict, 
                'correction': {'name': 'upweight'}
                }
        fc_mock = UpweightCorr(fc_cat_corr, **self.kwargs) 
        fc_mock_file = fc_mock.file()
        fc_mock_cols = fc_mock.datacolumns()
        
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

        upw = np.where(fc_mock.wfc > 1)  # upweighted
        collided = np.where(fc_mock.wfc == 0) # collided 
        # number of collided galaxy = number of fiber collided pairs 
        n_fcpair = len(collided[0])
        
        # Calculate the LOS displacement for the collided pairs
        comdis_upw = cosmos.distance.comoving_distance(
                fc_mock.zupw[collided], **cosmo) * cosmo['h']
        comdis_coll = cosmos.distance.comoving_distance(
                fc_mock.z[collided], **cosmo) * cosmo['h']
        dlos = np.abs(comdis_coll - comdis_upw) 
        
        # dLOS within 3 sigmas (peak galaxies)
        notin_tail = (collided[0])[np.where(dlos <= 3.0*corrdict['sigma'])]
        print np.float(len(notin_tail))/np.float(len(collided[0]))

        # expected number of galaxies to be placed in the peak 
        # of the dLOS distribution calculated based on the peak 
        # fraction. 
        #n_peak_exp = int(f_peak * np.float(n_fcpair))
        np.random.shuffle(notin_tail)
        
        # collided galaxies that will be peak corrected
        #i_peakcorr = notin_tail[:n_peak_exp]
        i_peakcorr = notin_tail

        fc_mock.wfc[i_peakcorr] += 1.0
        fc_mock.wfc[fc_mock.upw_index[i_peakcorr]] -= 1.0
        
        # comoving distances of upweighted galaxies in fc pairs 
        # that are going to be peakcorrected
        comdis_upw_gal = cosmos.distance.comoving_distance(
                fc_mock.zupw[i_peakcorr], **cosmo) * cosmo['h']

        # sampled dLOS from dLOS peak distribution best-fit function 
        d_samples = np.array([
            sample_dlos_peak(corrdict['fit'], corrdict['sigma'])[0]
            for i in range(len(i_peakcorr))
            ])

        outofbounds = np.where( 
                (comdis_upw_gal + d_samples > survey_comdis_max) |
                (comdis_upw_gal + d_samples < survey_comdis_min)
                )

        if len(outofbounds[0]) > 0: 
            d_samples[outofbounds] *= -1.0

        collided_z = comdis2z( comdis_upw_gal + d_samples) 

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
