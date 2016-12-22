'''

Photometric Redshift Correction 


'''
import numpy as np
import os 
# --- Local ---
from util.direc import direc
from corrections import Corrections
from fibcollided import UpweightCorr 
from ChangTools.fitstables import mrdfits

class PhotozCorr(Corrections): 

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

        super(PhotozCorr, self).__init__(cat_corr, **kwargs)
        
        self.corrstr() 
    
    def corrstr(self): 
        """ Specify correction string for dlos peak correction 
        """
        corrdict = (self.cat_corr)['correction']

        self.corr_str = ''.join([
            '.', corrdict['name']
            ])
        
        return self.corr_str

    def build(self): 
        """ Build mock catalogs with simulated photometric redshift that mimic
        the photometric redshifts of CMASS. They match the distribution of 
        delta z / (1+z) vs z of the CMASS BOSS photometric redshift training set.  
        """

        catdict = self.cat_corr['catalog']

        if catdict['name'] == 'cmass': 
            raise NameError("Only available for mock catalogs NOT CMASS") 
        
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
                skiprows = 1, 
                unpack = True, 
                usecols = xrange(len(fc_mock_cols))
                )
        
        for i_col, fc_col in enumerate(fc_mock_cols): 
            setattr(fc_mock, fc_col, fc_data[i_col])
    
        # fiber collided galaxies
        coll = np.where(
                fc_mock.wfc == 0
                )      

        # summary statistics (mu and sigma) of delta z / (1+z)
        z_mid, z_low, z_high, mu_deltaz, sigma_deltaz = cmass_deltaz_zspec_zphoto()

        photoz = np.array([-999. for i in xrange(len(fc_mock.ra))])
        
        deltaz_z = []
        for i_coll in coll[0]: 

            i_zbin_closest = np.argmin( np.abs(z_mid - fc_mock.z[i_coll]) )

            delz = np.random.normal(
                    loc = mu_deltaz[i_zbin_closest], 
                    scale = sigma_deltaz[i_zbin_closest]) * (1. + fc_mock.z[i_coll])

            photoz[i_coll] = fc_mock.z[i_coll] - delz

            deltaz_z.append(delz / (1.+fc_mock.z[i_coll]))

        data_cols = self.datacolumns()
        data_fmts = self.datacols_fmt()
        data_hdrs = self.datacols_header()
        data_list = []  

        for data_col in data_cols: 
            
            if data_col == 'photoz': 
                new_col = photoz
            else: 
                new_col = getattr(fc_mock, data_col)

            data_list.append(new_col)

        # write to corrected data to file 
        output_file = self.file()
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt = data_fmts, 
                delimiter = '\t', 
                header = data_hdrs
                ) 

    def datacolumns(self): 
        # Data columns
        cols = super(PhotozCorr, self).datacolumns()
        if 'photoz' not in cols: 
            cols += ['photoz']

        return cols 

    def datacols_fmt(self): 
        # Data column formats
        cols_fmt = super(PhotozCorr, self).datacols_fmt()
        cols_fmt += ['%10.5f']

        return cols_fmt

    def datacols_header(self): 
        # Data column headers 
        cols_header = super(PhotozCorr, self).datacols_header()
        cols_header += ", photoz"

        return cols_header 
    
def cmass_deltaz_zspec_zphoto():
    """ Summary statistics of delta z/(1+z) as a function of z.  
    Mu and stddev of delta z/(1+z) for bins of redshift. 

    THIS IS COPIED FROM photoz.photoz module
    """

    deltaz_file = ''.join([
        '/mount/riachuelo1/hahn/photoz/',
        'cmass_deltaz_zspec_zphotoz_gauss_stat.dat'
        ]) 

    if not os.path.isfile(deltaz_file):  
        raise NameError('check out photoz.photoz')

    z_mid, z_low, z_high, mu_deltaz, sigma_deltaz = np.loadtxt(
            deltaz_file, 
            skiprows = 1,
            unpack = True, 
            usecols = [0,1,2,3,4]
            ) 

    return [z_mid, z_low, z_high, mu_deltaz, sigma_deltaz]
