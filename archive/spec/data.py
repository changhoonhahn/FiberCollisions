'''

Code to handle galaxy data for calculating power/bispectrum 

Author(s): ChangHoon Hahn 


'''


import numpy as np
import scipy as sp 
import time 
import random
import os.path
import subprocess
import cosmolopy as cosmos

# --- Local ----
from util.catalog import Catalog 

# Classes ------------------------------------------------------------
class Data(object): 

    def __init__(self, cat_corr, **kwargs): 
        """ A class describing galaxy catalog of simulations or BOSS data. 
        Corresponds to an ASCII file with galaxy/random catalog 

        Parameters
        ----------
        cat_corr :  Catalog correction Dictionary 

        """ 

        self.cat_corr = cat_corr.copy()    
        self.kwargs = kwargs    

        self.file_name = None

        # catalog 
        self.catclass = Catalog(self.cat_corr)
        self.corr_str = None

        # galaxy properties
        self.ra = None
        self.dec = None
        self.z = None
        self.wfc = None
        self.weight = None 
        self.comp = None

        self.cosmos = None   # cosmology of catalog 

    def read(self): 
        """ Read galaxy/random catalog data 
        """
        if self.file_name is None: 
            self.file_name = self.file()
        
        data_cols = self.datacolumns()
        self.datacolumns = self.datacolumns()
    
        datah = np.loadtxt(
                self.file_name, 
                skiprows=1, 
                unpack=True, 
                usecols=range(len(data_cols))
                )

        for i_col, col in enumerate(data_cols): 
            setattr(self, col, datah[i_col])
        
        return None
    
    def build(self): 
        ''' Construct the default galaxy catalog. 

        Notes
        -----
        * Tiling Mock : Reads in true mock catalogs and imposes predetermined redshift 
        limits on them and attaches a flag mostly hardcoded since it's a simple procedure 
        * QPM: Handles messy weights 
        * PATCHY: Returns mocks with only necessary columns and w_fc = 1
        * Everything is very hardcoded
        '''

        catalog = (self.cat_corr)['catalog']

        output_file = self.file()
        
        if catalog['name'].lower() == 'nseries':       # N Series
            
            # read in original files from Jeremy and adjust them to make them
            # easier to use for fiber collisions
            # read rdzw file 
            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            orig_file = ''.join([
                data_dir, 
                'CutskyN', str(catalog['n_mock']), '.rdzwc'
                ]) 
            orig_ra, orig_dec, orig_z, orig_wfc, orig_zupw, orig_upw_index = np.loadtxt(
                    orig_file, 
                    unpack = True, 
                    usecols=[0,1,2,4,5,6]
                    )
        
            # file with completeness
            mask_file = ''.join([
                data_dir, 
                'CutskyN', str(catalog['n_mock']), '.mask_info'
                ]) 

            orig_wcomp = np.loadtxt(
                    mask_file, 
                    unpack = True, 
                    usecols = [0]
                    ) 

            # true wfc 
            true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
            
            # write to file 
            data_list = [orig_ra, orig_dec, orig_z, true_wfc, orig_wcomp, orig_zupw, orig_upw_index]

        elif catalog['name'].lower() == 'lasdamasgeo':          
            
            orig_true_file = ''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/', 
                'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.dat']) 
            orig_true_data = np.loadtxt(orig_true_file, unpack=True, usecols=[0,1,2])         # ra, dec, ***CZ***

            true_ra = orig_true_data[0]
            true_dec = orig_true_data[1]
            true_z = orig_true_data[2]/299800.0         # convert cz to z

            true_weight = np.array([1.0 for j in range(len(true_z))])   # no weights for true (all 1) 

            header_str = "Columns : ra, dec, z, weight"
            data_list = [true_ra, true_dec, true_z, true_weight]
            data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']

        elif catalog['name'].lower() == 'patchy':               # PATCHY mocks ------------------------

            # read original mock data 
            orig_file = ''.join(['/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/', 
                'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                str("%04d" % catalog['n_mock']), '.dat']) 

            # ra, dec, z, nbar, wfc, veto 
            orig_ra, orig_dec, orig_z, orig_nbar, orig_veto = np.genfromtxt(orig_file, 
                    unpack=True, usecols=[0, 1, 2, 4, 6]) 
            n_gal = len(orig_ra) 
            
            new_wfc = np.array([1.0 for i in range(n_gal)])     # w_fc = 1.0 for true 

            vetomask = (orig_veto == 1)            # only keep galaxies with w_veto = 1
           
            header_str = "Columns : ra, dec, z, nz, wfc" 
            data_list = [orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                        orig_nbar[vetomask], new_wfc[vetomask]] 
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']
        
        elif 'bigmd' in catalog['name'].lower():                # Big MultiDark ------------

            P0 = 20000.0

            # read rdzw file 
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if catalog['name'].lower() == 'bigmd':
                orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.dat'])  # hardcoded
            elif catalog['name'].lower() == 'bigmd1': 
                orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.dat'])
            elif catalog['name'].lower() == 'bigmd2': 
                orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.dat'])
            elif catalog['name'].lower() == 'bigmd3': 
                orig_file = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.dat']) 
            else: 
                raise NotImplementedError('asdlfkjadf') 

            orig_ra, orig_dec, orig_z, orig_wfkp, orig_veto, orig_wfc = np.loadtxt(orig_file, 
                    unpack=True, usecols=[0,1,2,3,4,5])

            # true wfc = 1 for all galaxies 
            true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
            nbar = (1.0/P0) * (1.0/orig_wfkp - 1.0) 

            vetomask = np.where(orig_veto == 1)     # if veto = 1 then keep; otherwise discard 
            
            # write to file 
            header_str = "Columns : ra, dec, z, nz, wfc"
            data_list = [orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], nbar[vetomask], true_wfc[vetomask]]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']

        else: 
            raise NameError('not yet coded') 
            
        # data columns, format and header are predetermined. Consult
        # corrections.corrections and util.catalog for specifics
        data_fmt = self.datacols_fmt()
        header_str = self.datacols_header()
        
        # write to file 
        np.savetxt(
                self.file(), 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmt, 
                delimiter='\t', 
                header=header_str
                ) 

        return None 

    def file(self): 
        """ Name of ASCII file of Data/Random catalogy
        """

        file_list = self.catclass.file()

        # cosmology string 
        try: 
            if self.kwargs['cosmology'] == 'survey': 
                cosmos_str = '.sureycosmo'
            else: 
                cosmos_str = '.fidcosmo'
        except KeyError: 
            cosmos_str = '.fidcosmo'
        file_list.insert( -1, cosmos_str )

        # need to account for the fact that the default 

        # correction string 
        if self.corr_str is not None: 
            file_list.insert( -1, self.corr_str )

            return ''.join(file_list)
        else: 
            return ''.join(file_list)

    def cosmo(self): 
    
        try: 
            if self.kwargs['cosmology'] == 'survey': 
                # survey cosmology
                self.cosmo = self.catclass.cosmo()

                return self.cosmo
            else: 
                # default fiducial cosmology (hardcoded)
                omega_m = 0.31 

        except KeyError: 
            omega_m = 0.31  # default 

        # survey cosmology 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmos = cosmo 

        return self.cosmos
    
    def survey_zlimits(self): 
        """ Catalog survey limits
        """
        return (self.catclass).survey_zlimits()
    
    def datacolumns(self): 
        ''' 
        Data columns for given catalog and correction
        '''
        return (self.catclass).datacolumns()

    def datacols_fmt(self): 
        ''' 
        Data format of columns of catalog data
        '''
        return (self.catclass).datacols_fmt()

    def datacols_header(self): 
        ''' 
        Header string that describes data columsn
        '''
        return (self.catclass).datacols_header()
