'''

FiberCollision Correction using peak of line-of-sight 
displacement distribution and photometric redshift.


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
from photoz_corr import PhotozCorr 
from galenv.galenv import d_NN_dataclass
from dlospeak import sample_dlos_peak
from dlospeak import peak_fit_gauss 
from dlospeak import peak_fit_expon 
from dlospeak import temp_nbarz 

class DlospeakPhotozCorr(Corrections): 

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

        super(DlospeakPhotozCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        """ Specify correction string for dlos peak correction 
        d_photoz_tail_cut has to be specified
        """
        corr = (self.cat_corr)['correction']

        if not all(x in corr.keys() for x in ['fit', 'fpeak', 'sigma', 'd_photoz_tail_cut']): 
            raise KeyError("Specify fpeak, sigma, and fit in correction dictionary")

        self.corr_str = ''.join([
            '.', corr['fit'].lower(), 
            '.', corr['name'].lower(),
            '.dphotozcut', str(corr['d_photoz_tail_cut']), 
            '.sigma', str(corr['sigma']), '.fpeak', str(corr['fpeak'])
            ])
        return self.corr_str

    def build(self): 
        """ Build fiber collision corrected mock catalogs. Using the modeled dLOS peak 
        and modeled photometric redshfit of upweighted galaxies. Note that this method reads in 
        photoz assigned upweight corrected mock catalogs and if care is not given, may result 
        in circular dependencies. 
        """

        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']

        cosmo = self.cosmo()      # cosmoslogy 

        f_peak = corrdict['fpeak']
        dc_photoz_tailcut = corrdict['d_photoz_tail_cut'] 

        # survey redshift limits  
        survey_zmin, survey_zmax = self.survey_zlimits()  

        survey_comdis_min = cosmos.distance.comoving_distance( survey_zmin, **cosmo ) * cosmo['h']
        survey_comdis_max = cosmos.distance.comoving_distance( survey_zmax, **cosmo ) * cosmo['h']
    
        # spline interpolation function hardcoded here 
        # to make it faster
        z_arr = np.arange(0.0, 1.01, 0.01)
        dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo) * cosmo['h']
        comdis2z = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 

        # photoz assigned upweight corrected galaxy catalog is used 
        # in the correction. This is out of convenience and 
        # also due to the fact that the BOSS observed
        # galaxy catalog outputs an upweigh corrected galaxy 
        # catalog. 
        fc_cat_corr = {
                'catalog': (self.cat_corr)['catalog'], 
                'correction': {'name': 'photoz'}
                }
        fc_mock = PhotozCorr(fc_cat_corr, **self.kwargs) 
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
        
        # nbar(z) cubic spline interpolation, which is currently hardcoded
        # in the function temp_nbarz(). However, this should 
        # ultimately be edited so that it works for any catalog 
        # (currently only works for CMASS-like catalogs)
        if 'cmass' in catdict['name'].lower(): 

            nb_z, nb_nbar = temp_nbarz(self.cat_corr)

            nbarofz = sp.interpolate.interp1d(nb_z, nb_nbar, kind='cubic')       

        upw = np.where(fc_mock.wfc > 1)  # upweighted
        collided = np.where(fc_mock.wfc == 0) # collided 
        
        # number of collided galaxy = number of fiber collided pairs 
        n_fcpair = len(collided[0])

        # Difference in comoving distance of upweight galaxy redshift and the 
        # comoving distance of the photometric redshift. The main idea is that
        # this distance roughly traces the actual dLOS 
        # dLOS_photo = Dc(z_upw) - Dc(z_photo)
        comdis_upw = cosmos.distance.comoving_distance(
                fc_mock.zupw[collided], **cosmo) * cosmo['h']
        comdis_coll = cosmos.distance.comoving_distance(
                fc_mock.z[collided], **cosmo) * cosmo['h']
        comdis_zphoto = cosmos.distance.comoving_distance(
                fc_mock.photoz[collided], **cosmo) * cosmo['h']

        dlos_photoz = np.abs( comdis_zphoto - comdis_upw ) 
        dlos_actual = np.abs( comdis_coll - comdis_upw) 

        notin_tail_photoz = (collided[0])[np.where(dlos_photoz <= dc_photoz_tailcut)]
        in_tail_photoz = np.where(dlos_photoz > dc_photoz_tailcut)
        #notin_tail_photoz = (collided[0])[np.where(dlos_actual <= dc_photoz_tailcut)]
        #in_tail_photoz = np.where(dlos_actual> dc_photoz_tailcut)

        dlos_in_tail_photoz = dlos_actual[in_tail_photoz]
        print 'Contamination rate ', np.float(len(np.where(dlos_in_tail_photoz < 3.0*corrdict['sigma'])[0]))/\
                np.float(len(dlos_in_tail_photoz))*100
         
        # expected number of galaxies to be placed in the peak 
        # of the dLOS distribution calculated based on the peak 
        # fraction. 
        n_peak_exp = int(f_peak * np.float(n_fcpair))

        np.random.seed() 
        np.random.shuffle(notin_tail_photoz)
        
        # collided galaxies that will be peak corrected
        i_peakcorr = notin_tail_photoz[:n_peak_exp]

        fc_mock.wfc[i_peakcorr] += 1.0
        for upw_index in fc_mock.upw_index[i_peakcorr]: 
            fc_mock.wfc[upw_index] -= 1.0

        fc_mock.ra[i_peakcorr] = fc_mock.ra[fc_mock.upw_index[i_peakcorr]]
        fc_mock.dec[i_peakcorr] = fc_mock.dec[fc_mock.upw_index[i_peakcorr]]
        
        # comoving distances of upweighted galaxies in fc pairs 
        # that are going to be peakcorrected
        comdis_upw_gal = cosmos.distance.comoving_distance(
                fc_mock.zupw[i_peakcorr], **cosmo) * cosmo['h']

        # sampled dLOS from dLOS peak distribution best-fit function 
        d_samples = np.array([
            sample_dlos_peak(corrdict['fit'], corrdict['sigma'])[0]
            for i in range(n_peak_exp)
            ])

        outofbounds = np.where( 
                (comdis_upw_gal + d_samples > survey_comdis_max) |
                (comdis_upw_gal + d_samples < survey_comdis_min)
                )

        if len(outofbounds[0]) > 0: 
            d_samples[outofbounds] *= -1.0

        collided_z = comdis2z( comdis_upw_gal + d_samples) 

        fc_mock.z[i_peakcorr] = collided_z

        if 'cmass' in catdict['name'].lower():
            fc_mock.nbar[i_peakcorr] = nbarofz(collided_z)

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

    def datacolumns(self): 
        # Data columns
        cols = super(DlospeakPhotozCorr, self).datacolumns()
        if 'photoz' not in cols: 
            cols += ['photoz']

        return cols 

    def datacols_fmt(self): 
        # Data column formats
        cols_fmt = super(DlospeakPhotozCorr, self).datacols_fmt()
        cols_fmt += ['%10.5f']

        return cols_fmt

    def datacols_header(self): 
        # Data column headers 
        cols_header = super(DlospeakPhotozCorr, self).datacols_header()
        cols_header += ", photoz"

        return cols_header 

if __name__=='__main__': 
        
    fc_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'dlospeak_photoz', 'fit': 'gauss', 'fpeak':0.68, 'sigma':3.9, 'd_photoz_tail_cut':200}}
    test = DlospeakPhotozCorr(fc_cat_corr)
    test.build()
